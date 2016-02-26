import pandas as pd
from Model import *
from filelocations import *
import pylab as plt
import seaborn as sns

params = pd.read_csv(os.path.join(parameterfolder, 'kinetics_fit_summary.csv'), index_col=[0,1], header=0)

params = params.transpose()

a = np.array(params['a'].loc['mean',:])
k = np.array(params['k'].loc['mean',:])
b = float(params['b'].loc['mean',:])
u = float(params['u'].loc['mean',:])
n = float(params['n'].loc['mean',:])

e = evaluateSystem(a, k, b, u, n, kmodel=True, dosemodel=False)
experimental_time_to_half_max = e['ttss'] * 25.0/60.0
experimental_time_to_half_max.columns = ['induction layer', 'Strain 6059', 'Strain 6325', 'Strain 6326', 'Strain 6327']

bestfit = pd.read_csv(os.path.join(fittingfolder, 'kinetics_fit.csv'))

def function(x, column_label, matrix):
    index = np.array(bestfit.index)
    col = index
    if not column_label == 'index':
        col = bestfit[column_label]
    i = np.argmin(np.abs(col - x))
    return matrix.loc[i,:]


############
layers = pd.read_csv(os.path.join(expdatafolder, 'kinetics_exp_data_all.csv'))
for row in range(len(layers)):
    t_exp = layers.loc[row, 'time']
    au_exp = layers.loc[row, 'au']
    layer = layers.loc[row, 'layer']
    expected = function(t_exp, 'time (min)', bestfit)
    expected.columns = ['time', 0, 1, 2, 3, 4]
    au_expected = expected[layer+1]
    residual = au_exp - au_expected
    layers.loc[row, 'residual'] = residual

layers.to_csv(os.path.join(fittingfolder, 'kinetics_residuals.csv'))
std_residual = []
for g in layers.groupby('layer'):
    plt.clf()
    sns.set_context('poster')
    sns.set_style('whitegrid', {'axes.linewidth': 10.0,})
    sns.set_style('ticks')
    plt.plot(g[1]['time'], g[1]['residual'], marker='o', c='black', linestyle='None', ms=8)
    sns.despine()
    plotfigpath = os.path.join(fittingfolder, 'kinetics_residual_layer{}'.format(g[0]))
    plt.savefig(plotfigpath + ".pdf", format="pdf")
    plt.clf()

    std = np.std(g[1]['residual'])
    std_residual.append(std)


bestfit.columns = ['time', 0, 1, 2, 3, 4]
t_half_array = np.array(experimental_time_to_half_max)[0][1:] * 60.
upper_error_array = []
lower_error_array = []
for i, t_half in enumerate(t_half_array):
    r = function(t_half, 'time', bestfit)
    au_experimental = np.array(r)[2:]
    au_p = au_experimental + std_residual
    au_m = au_experimental - std_residual
    r_p = function(au_p[i], i+1, bestfit)
    r_m = function(au_m[i], i+1, bestfit)
    error = [r_p['time'], r_m['time']]
    value = np.array([t_half, np.abs(t_half - min(error)), np.abs(t_half - max(error))])
    value = value * 1/60.0
    experimental_time_to_half_max['Layer {} lower error'.format(i+1)] = value[1]
    experimental_time_to_half_max['Layer {} upper error'.format(i+1)] = value[2]
    upper_error_array.append(value[2])
    lower_error_array.append(value[1])
p = pd.DataFrame(
        {"time to half max (hours)": t_half_array / 60.0,
        "upper error": upper_error_array,
        "lower error": lower_error_array
         })
p.index.name = "Layer"
p.to_csv(os.path.join(fittingfolder, 'experimental_time_to_half_max.csv'))
