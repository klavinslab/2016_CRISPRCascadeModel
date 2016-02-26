import os

################################### Define Filenames ###################################

root = os.path.dirname(__file__)
datafolder = os.path.join(root, "Data")
hyakfolder = os.path.join(datafolder, "hyakfitting_results")
finalfolder = os.path.join(root, "FinalData")
parameterfolder = os.path.join(finalfolder, 'parameters')
expdatafolder = os.path.join(datafolder, 'ExperimentalData')
fittingfolder = os.path.join(finalfolder, 'fitting')
montecarlofolder = os.path.join(finalfolder, "MonteCarlo")