import os
import glob
import pylab as plt
import seaborn as sns

from Model import *
from filelocations import *


################################### Output Files ###################################
# Data storing the best fits found with the hyak cluster for the knietics model
bestfits_kinetics_location = os.path.join(fittingfolder, "bestfits_kinetics_indexed.csv")

# Data for the statistics on the parameters for the kinetics model
kinetics_parameter_summary_location = os.path.join(parameterfolder, 'kinetics_fit_summary.csv')

# Data for plotting the expected fits for the kinetic model
kinetics_fit_location = os.path.join(fittingfolder, 'kinetics_fit.csv')

################################### Get Best Fits ###################################
scripts = os.path.join(datafolder, "HyakFittingScripts_CascadeAnalysis_20160104", "output_fittingscript_steadystate_20160104")
files = glob.glob(os.path.join(hyakfolder, 'HyakFittingScripts_CascadeAnalysis_20160104/output_fittingscript_steadystate_20160104/', "*.csv"))
best = []
df = None
for f in files:
    df = pd.read_csv(f, header=[0,1], index_col=0)
    lastrow = df.iloc[-1]
    best.append(lastrow)
best = pd.DataFrame(best)
best = best.reset_index()
u = best['u']
b = best['b']
best.drop('u', axis=1, level=0, inplace=True)
best.drop('b', axis=1, level=0, inplace=True)
best['u','0'] = u['0']
best['b','0'] = b['0']
best.to_csv(bestfits_kinetics_location)

summary = best.describe()
summary.transpose().to_csv(kinetics_parameter_summary_location)

################################### Import Kinetics Data ###################################
layers = []
for i in range(1,5):
    name = os.path.join(expdatafolder, "KineticsLayer%d.csv" % i)
    layers.append(pd.read_csv(name, header=None))
layers1, layers2, layers3, layers4 = layers
header = ['time', 'au']
layers1.columns = header
layers2.columns = header
layers3.columns = header
layers4.columns = header
layers1['layer'] = 1
layers2['layer'] = 2
layers3['layer'] = 3
layers4['layer'] = 4
layers = pd.concat([layers1, layers2, layers3, layers4])
layers_table = pd.pivot_table(layers, values='au', columns=['layer'], index='time')
layers.sort(columns=['time', 'layer'])
layers['au'] = layers['au']
layers.stack()
layers
exp_kinetics = pd.pivot_table(layers, values='au', index=['layer', 'time'])
layers.to_csv(os.path.join(expdatafolder, 'kinetics_exp_data_all.csv'))
# plt.figure()
# plt.xlabel("Time (min)")
# plt.ylabel("Fluorescence (AU)")
# plt.title("Experimental Kinetics")
# exp_kinetics[1].plot()
# exp_kinetics[2].plot()
# exp_kinetics[3].plot()
# exp_kinetics[4].plot()
# plt.show()


print "Plotting Best Kinetics Fits"
################################### Plot Best ###################################

def plotBest(infile, kout = None):
    params = infile
    params = params.loc['50%']
    a = list(params['a'])
    k = list(params['k'])
    u = float(params['u'])
    b = params['b']
    n = float(params['n'])
    b1 = b2 = b
    if len(b) == 2:
        b1, b2 = b
    else:
        b1 = float(b)
        b2 = b1
    k = kineticsModel(a, k, b1, u, n) * 1600.0
    k = k.set_index(k.index * 25.0)

    # Save files
    if not kout is None:
        k_relabeled = k.copy()
        k_relabeled.columns = ['induction layer', 'Strain 6059', 'Strain 6325','Strain 6326', 'Strain 6327']
        k_relabeled.index.name = 'time (min)'
        k_relabeled.to_csv(kout)
        sns.set_context('poster')
        sns.set_style('whitegrid', {'axes.linewidth': 10.0,})
        sns.set_style('ticks')
        plt.plot(k.index, k[1], marker=None, c='black')
        plt.plot(k.index, k[2], marker=None, c='r')
        plt.plot(k.index, k[3], marker=None, c='blueviolet')
        plt.plot(k.index, k[4], marker=None, c='mediumseagreen')
        e1 = pd.DataFrame(exp_kinetics[1])
        e2 = pd.DataFrame(exp_kinetics[2])
        e3 = pd.DataFrame(exp_kinetics[3])
        e4 = pd.DataFrame(exp_kinetics[4])
        plt.plot(e1.index, e1['au'], marker='o', c='black', linestyle='None', ms=8)
        plt.plot(e2.index, e2['au'], marker='o', c='r', linestyle='None', ms=8)
        plt.plot(e3.index, e3['au'], marker='o', c='blueviolet', linestyle='None', ms=8)
        plt.plot(e4.index, e4['au'], marker='o', c='mediumseagreen', linestyle='None', ms=8)
        plt.xlim(0,2500)
        sns.despine()
#         exp_kinetics[1].plot()
#         exp_kinetics[2].plot()
#         exp_kinetics[3].plot()
#         exp_kinetics[4].plot()
        plotfigpath = os.path.join(os.path.dirname(kout), os.path.basename(kout).split(".")[0])
        plt.savefig(plotfigpath + ".png", format="png")
        plt.savefig(plotfigpath + ".pdf", format="pdf")
        return k


k_in = summary
k_out = os.path.join(fittingfolder, 'kinetics_fit.csv')
plotBest(k_in, kout=k_out)