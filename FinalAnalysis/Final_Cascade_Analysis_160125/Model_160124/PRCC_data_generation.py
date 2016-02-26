from filelocations import *
from Model import *
import pandas as pd

numlayers = 7
param_range = {
       'a': (0.3929,1.408),
       'k': (0.3817,1.041),
       'u': (.2073,.5547),
        'b': (.09550,.1430),
        'n': (1.454,2.015)
      }

param_range_summary = pd.DataFrame(param_range).transpose()
param_range_summary.columns = ['low', 'high']
param_range_summary.index.name = 'parameter'

################################### Output Files ###################################

''' File for storing the range used to generate the PRCC data '''
param_range_summary.to_csv(os.path.join(montecarlofolder, 'parameter_range_for_prcc.csv'))


''' File for storing the data for the PRCC analysis '''
filename_raw = os.path.join(montecarlofolder, "PRCC_uniform_raw.csv")

################################### Data Generation ###################################
# resample_parameters_and_evaluate(numlayers,
#                                  input_filename = None,
#                                  resample = False,
#                                  parameter_bounds = param_range,
#                                  num_rounds=100,
#                                  outfile=filename_raw,
#                                  overwrite=True,
#                                  n_bounds=None,
#                                  n_force=None,
#                                  dosemodel=False,
#                                  cleanup=True)
# for i in range(99):
#     resample_parameters_and_evaluate(numlayers,
#                                      input_filename = None,
#                                      resample = False,
#                                      parameter_bounds = param_range,
#                                      num_rounds=100,
#                                      outfile=filename_raw,
#                                      overwrite=False,
#                                      n_bounds=None,
#                                      n_force=None,
#                                      dosemodel=False,
#                                      cleanup=True)
