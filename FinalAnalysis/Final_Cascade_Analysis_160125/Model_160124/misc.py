def generateBootstrapForModel(infile, outfile, reps=0, rounds=10, overwrite=False):
    for i in range(reps):
        modelbootstrap(infile, overwrite=overwrite, rounds=rounds, outfile=outfile)

def makeSummary(outfile):
    bootmodelresults = pd.read_csv(outfile, header=[0,1], index_col=0)
    summary = bootmodelresults.describe()
    summary.to_csv(os.path.join(os.path.dirname(outfile), "summary_" + os.path.basename(outfile)))

##################### Kinetics Model Bootstrap ######################
k_model_boot=os.path.join(montecarlofolder, "modelbootstrap_kinetics.csv")
generateBootstrapForModel(kinetics_new_filename, k_model_boot, reps=0, rounds=0, overwrite=False)
makeSummary(k_model_boot)