import pylab as plt
import glob
import seaborn as sns
from Model import *
from filelocations import *


# Description: This script processes the dose response experimental data as well as the fitting
# results obtained from the Hyak cluster. Plots the experimental and expected fits
# using Seaborn.

################################### Output Files ###################################

# Processed csv's containing the first 48 parameter
# fits obtained from the Hyak cluster
best1_location = os.path.join(fittingfolder, 'bestfits_steadystate_indexed_1.csv')
best2_location = os.path.join(fittingfolder, 'bestfits_steadystate_indexed_2.csv')
best3_location = os.path.join(fittingfolder, 'bestfits_steadystate_indexed_3.csv')

# Raw csv's containing the all parameter
# fits obtained from the Hyak cluster
best3_raw_location = os.path.join(fittingfolder, 'bestfits_steadystate_indexed_3_raw.csv')
best1_raw_location = os.path.join(fittingfolder, 'bestfits_steadystate_indexed_1_raw.csv')

# Statistics on the parameter
# fits obtained from the Hyak cluster
bestfits_summary_1_location = os.path.join(fittingfolder, "BestFitsSummary_1.csv")
bestfits_summary_2_location = os.path.join(fittingfolder, "BestFitsSummary_2.csv")
bestfits_summary_3_location = os.path.join(fittingfolder, "BestFitsSummary_3.csv")
dose_response_fit_parameter_summary_location = os.path.join(parameterfolder, 'dose_response_fit_summary.csv')

# Data that summarized the experimental data
# for dose response experiments 1 through 3
dose_response_alldata_location = os.path.join(finalfolder, 'steady_state_alldata.csv')
dose_response_experimental_means_location = os.path.join(finalfolder, 'steady_state_alldata.csv')
dose_response_experimental_std_location = os.path.join(finalfolder, 'steady_state_std.csv')

# Data obtained from using the parameters
# and the cascade model. Used to plot figure.
dose_response_fit_location = os.path.join(fittingfolder, 'dose_response_fit.csv')

# Location of the dose response plot as a pdf
dose_response_plot_location = os.path.join(fittingfolder, 'dose_response_figure.pdf')


############ Process Hyak Fitting Result ###############

exp1 = os.path.join(hyakfolder, 'HyakFittingScripts_CascadeAnalysis_20151009', 'output_fittingscript_steadystate_100115')
exp2 = os.path.join(hyakfolder, 'HyakFittingScripts_CascadeAnalysis_20151209', 'output_fittingscript_steadystate_20151209_exp1')
exp3 = os.path.join(hyakfolder, 'HyakFittingScripts_CascadeAnalysis_20151209', 'output_fittingscript_steadystate_20151209_exp2')

exp1_solutions = glob.glob(os.path.join(exp1, "solution*.csv"))
exp2_solutions = glob.glob(os.path.join(exp2, "solution*.csv"))
exp3_solutions = glob.glob(os.path.join(exp3, "solution*.csv"))

def getBestFits(list_of_files):
    best = []
    for f in list_of_files:
        df = pd.read_csv(f, header=[0,1], index_col=0)
        lastrow = df.iloc[-1]
        best.append(lastrow)
    best = pd.DataFrame(best)
    best = best.reset_index()
    return best

best1 = getBestFits(exp1_solutions)
best2 = getBestFits(exp2_solutions)
best3 = getBestFits(exp3_solutions)

best1.to_csv(best1_location)
best2.to_csv(best2_location)
best3.to_csv(best3_location)


### Note: There is a very clear bimodal distribution within exp_3_solutions
### such that a large subset has a much higher score. I will be removing these
### solutions. This results in 48 scores for best3, 48 scores for best2, and
### 100 scores for best1. I will take the first 48 scores of best1 [arbitrary selection]
### resulting in 144 total fits to be used in downstream analysis.

best3.to_csv(best3_raw_location)
best3_raw = best3
best3 = best3[best3['score']['0'] < 1]
best3.to_csv(best3_location)

best1.to_csv(best1_raw_location)
best1_raw = best1
best1 = best1.iloc[:48]
best1.to_csv(best1_location)

best144 = pd.concat([best1, best2, best3])

def importData(file_location, plot=True):
    exp = pd.read_csv(file_location)
    AU = exp['FL1.Amean']
    treatment = exp['treatment']
    exp_dose = pd.pivot_table(exp, index='treatment', columns=['strain'], values='FL1.Amean')
    if plot:
        ax = exp_dose.plot(logx=True, logy=False)
        ax.set_xlabel("Beta-estradiol (nM)")
        ax.set_ylabel("Fluorescence (AU)")
        ax.set_title("Experimental Does Response")
        ax
    return exp_dose

best1.describe().to_csv(bestfits_summary_1_location)
best2.describe().to_csv(bestfits_summary_2_location)
best3.describe().to_csv(bestfits_summary_3_location)
m = [best1.describe().loc['mean'], best2.describe().loc['mean'], best3.describe().loc['mean']]
m = pd.DataFrame(m)
param_summary = pd.concat([pd.DataFrame(m.mean()), pd.DataFrame(m.std())], axis=1)
param_summary.columns = ['mean', 'std']
param_summary = param_summary.transpose()
u = param_summary['u']
b = param_summary['b']
param_summary.drop('u', axis=1, level=0, inplace=True)
param_summary.drop('b', axis=1, level=0, inplace=True)
param_summary['u','0'] = u['0']
param_summary['b','0'] = b['1']
param_summary = param_summary.transpose()
param_summary.to_csv(dose_response_fit_parameter_summary_location)


#################### Process Dose Response Data #######################
plotting = False
exp1 = importData(os.path.join(expdatafolder, '20150905_184902_CascadeDoseResponse_Final.csv'), plot=plotting)
exp2 = importData(os.path.join(expdatafolder, '20151208_195615_doseresponse_plate1.csv'), plot=plotting)
exp3 = importData(os.path.join(expdatafolder, '20151208_211350_doseresponse_plate2.csv'), plot=plotting)

exp1['exp'] = 1
exp2['exp'] = 2
exp3['exp'] = 3
allexp = pd.concat([exp1, exp2, exp3])
allexp['treatment'] = allexp.index
means = allexp.groupby('treatment').mean()
std = allexp.groupby('treatment').std()
means.drop('exp', axis=1, inplace=True)
std.drop('exp', axis=1, inplace=True)

allexp.to_csv(dose_response_alldata_location)
means.to_csv(dose_response_experimental_means_location)
std.to_csv(dose_response_experimental_std_location)

best_params = param_summary.transpose()

best_params_mn = best_params.loc['mean']
best_params_std = best_params.loc['std']
a = best_params_mn.a
k = best_params_mn.k
b = best_params_mn.b
u = best_params_mn.u
n = best_params_mn.n
p = (list(a), list(k), float(b), float(u), float(n))
refinedmodel = ssModel(0.006*2**np.linspace(0,20,100), *p) * 1600.0
# TODO: print residuals for steady state
# for g in allexp.groupby('exp'):
#     expected = ssModel(np.array(g[1].index), *p) * 1600.0
#     print expected
m = ssModel(means.index, *p) * 1600.0
plt.figure()
m.join(means).plot(logx=True)
means.plot(logx=True)
dose_response_fit = pd.DataFrame({'treatment': refinedmodel.index,
                                  '6059': refinedmodel[2],
                                  '6325': refinedmodel[3],
                                  '6326': refinedmodel[4],
                                  '6327': refinedmodel[5]})
dose_response_fit.set_index('treatment')
dose_response_fit.columns.name = 'layers'
dose_response_fit_copy = dose_response_fit.copy()
dose_response_fit_copy.drop('treatment', axis=1, inplace=True)
dose_response_fit_copy.index.name = 'treatment (Be uM)'
dose_response_fit_copy.columns = ['Strain {}'.format(x) for x in [6059, 6325, 6326, 6327]]
dose_response_fit_copy.to_csv(dose_response_fit_location)



#################### Plot Dose Response Data #######################
def plot(f, means, std):
    plt.xscale('log')
    plt.plot(f['treatment'], f['6059'], marker=None, c='black')
    plt.plot(f['treatment'], f['6325'], marker=None, c='r')
    plt.plot(f['treatment'], f['6326'], marker=None, c='blueviolet')
    plt.plot(f['treatment'], f['6327'], marker=None, c='mediumseagreen')
    capsize = 4
    capthick = 2
    fmt = 'o'
    elinewidth=0.0
    ms = 7.5
    plt.errorbar(means.index, means[6059],
                 yerr=std[6059],
                 fmt=fmt, c='black', capthick=capthick,
                 capsize=capsize, elinewidth=elinewidth,
                 ms=ms)
    plt.errorbar(means.index, means[6325],
                 yerr=std[6325],
                 fmt=fmt, c='r', capthick=capthick,
                 capsize=capsize, elinewidth=elinewidth,
                 ms=ms)
    plt.errorbar(means.index, means[6326],
                 yerr=std[6326],
                 fmt=fmt, c='blueviolet', capthick=capthick,
                 capsize=capsize, elinewidth=elinewidth,
                 ms=ms)
    plt.errorbar(means.index, means[6327],
                 yerr=std[6327],
                 fmt=fmt, c='mediumseagreen', capthick=capthick,
                 capsize=capsize, elinewidth=elinewidth,
                 ms=ms)
    plt.xlim(0.006, 1000)
    plt.ylim(0,20000)

sns.set_palette("Reds")
sns.set_context('poster')
sns.set_style('whitegrid', {'axes.linewidth': 10.0,})
sns.set_style('ticks')
plot(dose_response_fit, means, std)
sns.despine()
plt.savefig(dose_response_plot_location, dpi=400)
sns.axes_style()