import pylab as plt
import seaborn as sns

from Model import *
from filelocations import *


plt.clf()
plt.close()
##################### Projection Monte Carlo ######################
print "Generating data for time-to-half-max plot"
### This generates data for the time-to-half-max plot
overwrite = False
numlayers = 30
outfile=os.path.join(montecarlofolder, "bootstrapped_kinetics_extra.csv")
inputfile = os.path.join(fittingfolder, 'bestfits_kinetics_indexed.csv')

# bootstrap_kinetics_results = resample_parameters_and_evaluate(numlayers,
#                                                                   input_filename=inputfile,
#                                                                   num_rounds=100,
#                                                                   overwrite=True,
#                                                                   outfile=outfile,
#                                                                   dosemodel=False,
#                                                                   cleanup=True)
# for i in range(10):
#     print i
#     bootstrap_kinetics_results = resample_parameters_and_evaluate(numlayers,
#                                                                   input_filename=inputfile,
#                                                                   num_rounds=100,
#                                                                   overwrite=False,
#                                                                   outfile=outfile,
#                                                                   dosemodel=False,
#                                                                   cleanup=True)

print "Data generation for time-to-half-max plot complete"

##################### Final Data Cleanup and Summary ######################
read_p = pd.read_csv(outfile, header=[0,1], index_col=0)
ttss = read_p['ttss']
ttss = ttss.apply(lambda x: x * 25.)
ttss.columns.name = 'layer'
ttss.index.name = 'sample'

ttss = ttss.stack('layer')
ttss = ttss.reset_index()
ttss.columns = ['sample', 'layer', 'time-to-half-max (min)']
grouped = ttss.groupby('layer')
summary = grouped.describe()
summary.drop('sample', inplace=True, axis=1)
summary = summary.unstack(level=1)
summary['layer'] = np.array([int(x) for x in summary.index])
summary = summary.set_index('layer')
summary = summary.sort()
summary.to_csv(os.path.join(finalfolder, "time-to-half-maximal_summary.csv"))
ttss.to_csv(os.path.join(finalfolder, "time-to-half-maximal_row.csv"))
ttss.columns = ['sample', 'layer', 'time-to-half-max']
sns.boxplot(x='layer', y='time-to-half-max', data = ttss, fliersize=0.5)

plt.savefig(os.path.join(finalfolder, "time-to-half-maximal_boxplot.png"), format='png', dpi=400)
plt.savefig(os.path.join(finalfolder, "time-to-half-maximal_boxplot.pdf"), format='pdf', dpi=400)