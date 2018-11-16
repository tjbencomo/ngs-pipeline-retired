'''
Author: Tomas Bencomo
Script to run the Data Preprocessing Workflow. 2 run options. Don't specify a directory and simply
run this python script from whichever directory you want the cromwell-execution and log files to be stored.
Option 2 is you can run the script from anywhere and specify the directory where you'd like the logs and output
files to be saved with the command line argument -d (example: python data_preprocessing.py -d directory/to/save/logs)

This python script launches an sbatch job that runs the cromwell process, creating sbatch jobs as the workflow calls for them.
This sbatch job must stay alive the entire time cromwell is running the workflow, so it should be set for a very long time.
It does not need many resources as it is just monitoring the processes.
'''

import os
import subprocess
import argparse

def createParser():
    parser = argparse.ArgumentParser(description='Script to launch Cromwell Data Preprocessing Workflow')
    parser.add_argument('-d', '--directory', type=str, nargs=1, 
                        help='Specify directory to launch workflow from, saving execution files to this directory',
                        dest='directory')

    args = parser.parse_args()
    return args

def parseArgs():
    args = createParser()

    if args.directory is None:
        directory = os.getcwd()
    else:
        directory = ''.join(args.directory)

    if os.path.isdir(directory) is False:
        raise ValueError("Directory does not exist!")
    
    return {'directory' : directory}

def launchJob(args):
    # environment variables
    HOME = os.environ['HOME']
    USER = os.environ['USER']
    SCRATCH = os.environ['SCRATCH']
    PI_SCRATCH = os.environ['PI_SCRATCH']

    cromwell_log_directory = os.path.join(PI_SCRATCH, "cromwell-monitor-logs")
    if os.path.isdir(cromwell_log_directory) is False:
        os.mkdir(cromwell_log_directory)

    os.chdir(args['directory'])

    # sbatch directive variables
    jobName = 'cromwell_data_preprocessing'
    outputFile = '{}'.format(os.path.join(cromwell_log_directory, 'cromwell_data_preprocessing.%j.out'))
    errorFile = '{}'.format(os.path.join(cromwell_log_directory, 'cromwell_data_preprocessing.%j.err'))
    nodes = '1'
    memory = '8000' #in MB
    time = '3-00:00:00'
    email = '{}@stanford.edu'.format(USER)
    mailType = 'END'
    
    # command line execution variables
    softwareDirectory = '/home/groups/carilee/software'
    configPath = os.path.join(softwareDirectory, 'variant-discovery-pipeline', 'your.conf')
    jarPath = os.path.join(softwareDirectory, 'cromwell-35.jar')
    scriptPath = os.path.join(softwareDirectory, 'variant-discovery-pipeline', 'scripts', 'data_preprocessing.wdl')
    inputPath = os.path.join(softwareDirectory, 'variant-discovery-pipeline', 'inputs', 'data_preprocessing_inputs.json')
    wrap = 'java -Dconfig.file={} -jar {} run {} -i {}'.format(configPath, jarPath, scriptPath, inputPath)

    sbatchCommand = 'sbatch --job-name={} --output={} --error={} --nodes={} --mem={} --time={} --mail-user={} --mail-type={} --wrap="{}"'.format(jobName, outputFile, errorFile, nodes, memory, time, email, mailType, wrap)

    stdout = subprocess.run(sbatchCommand, shell=True, stdout=subprocess.PIPE)
    print(stdout.stdout)

def main():
    args = parseArgs()
    launchJob(args)

if __name__ == "__main__":
    main()