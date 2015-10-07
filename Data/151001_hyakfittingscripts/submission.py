##############################

#submission code by Willy Voje 151001

##############################

import os

ROOT = os.getcwd()



def create_parrel_submission_script(writepath, workingpath, sdoutpath,
                                    nameofexperiment, cores, nodes, backfill=False):
    """ -> .sh file for submission
    The call file command should be a string of txt that oulines
    what exactly is being called and how it is being modife
    This script writes a .sh file that will actually do the work of
    submitting a job to hyak
    """
    filepath = os.path.join(writepath, 'parallel_sql_job.sh')
    #open file
    f = open(filepath, 'wb')

    f.write("#!/bin/bash")
    # Defines the job jame RENAME FOR YOUR JOB
    f.write("\n#PBS -N " + nameofexperiment)

    # EDIT FOR YOUR JOB
    # For 16 core nodes.
    # Nodes should _never_ be > 1.
    if backfill:
        f.write("\n#PBS -l nodes=" + str(nodes) + ":ppn=" + str(cores) +
                ",mem=20gb,feature=" + str(cores) + "core -q bf")
    else:
        f.write("\n#PBS -l nodes=" + str(nodes) + ":ppn=" + str(cores) +
                ",mem=20gb,feature=" + str(cores) + "core")

    # WALLTIME DEFAULTS TO ONE HOUR - ALWAYS SPECIFY FOR LONGER JOBS
    # If the job doesn't finish in 10 minutes, cancel it
    f.write("\n#PBS -l walltime=48:00:00")

    # Put the STDOUT and STDERR from jobs into the below directory
    f.write("\n#PBS -o " + sdoutpath)
    # Put both the stderr and stdout into a single file
    f.write("\n#PBS -j oe")

    # Sepcify the working directory for this job bundle
    f.write("\n#PBS -d " + workingpath)

    # If you can't run as many tasks as there are cores due to memory constraints
    # you can simply set HYAK_SLOTS to a number instead.
    # HYAK_SLOTS=4
    f.write("\nHYAK_SLOTS=`wc -l < $PBS_NODEFILE`")

    f.write("\nmodule load parallel_sql")
    f.write("\nmodule load anaconda_2.3")
    f.write("\nparallel_sql --sql -a parallel --exit-on-term -j $HYAK_SLOTS")

    f.close()
    os.chmod(filepath, 0777)


def new_directory(directory):
    """
    (str)->directory
    simple code to make a directory in the root if it does not exist
    """
    if not os.path.exists(directory):
        os.makedirs(directory)


number_of_replicates = 100
# Make the list
with open('sql_stuff_to_add', 'wb') as submissionfile:

    script_names = ['fittingscript_combined_100115.py',
                    'fittingscript_kinetics_100115.py',
                    'fittingscript_steadystate_100115.py']

    for script_name in script_names:

        base_name = script_name.split('.py')[0]

        output_dirname = "output_" + base_name
        sdout_dirmane = "error_out_" + base_name

        # Make the output path

        output_path = os.path.join(ROOT, output_dirname)
        new_directory(output_path)

        # Make the STDOUT path

        sdout_path = os.path.join(ROOT, sdout_dirmane)
        new_directory(sdout_path)

        for i in range(number_of_replicates):
            output_location = os.path.join(output_path, "output-" + str(i))
            submissionfile.write("python " + script_name + " " + output_path + "\n")

# Make the submission script
create_parrel_submission_script(writepath=ROOT, workingpath=ROOT, sdoutpath=sdout_path,
                                nameofexperiment='test', cores=12, nodes=1,
                                backfill=True)



