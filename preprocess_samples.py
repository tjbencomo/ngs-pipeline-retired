import os
import subprocess
import argparse
import gzip
import re

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
        raise ValueError("{} directory does not exist!")
    
    return {'directory' : directory}

def getSamples(directory):
    '''
    Preconditions: 
        All files to be preprocessed end with .fastq.gz
        File pairs are formated such that pairs are distinguished by the suffix _X.fastq.gz
        such that CTR119_1.fastq.gz and CTR119_2.fastq.gz are a pair
    '''
    files = [f for f in os.listdir(directory) if os.path.isfile(os.path.join(directory, f)) and f.endswith('.fastq.gz')]

    samples = {}
    for f in files:
        if '_1.' in f:
            sampleName = f[0:f.find('_1.')]
            match = f[0:f.rfind('_1.')] + "_2" + f[f.find('.'):]
            if match in files:
                samples[sampleName] = {'fastq1' : os.path.join(directory, f), 'fastq2' : os.path.join(directory, match)} 

    return samples

def getReadGroupInfo(fastqFile, sample):
    '@HISEQ_173_C3LAWACXX:1:1101:1127:2116#0/1'
    f = gzip.open(fastqFile, 'r')
    headerLine = f.readline().decode("utf-8") 
    headerLine = re.split('_|:|#', headerLine)

    readGroupID = headerLine[2] + "." + headerLine[3]
    readGroupPlatformUnit = readGroupID + "." + sample
    readGroupPlatform = "ILLUMINA"

    return readGroupID, readGroupPlatformUnit, readGroupPlatform

def launchCromwellJob(sample, sampleFiles, directory):
    '''
    Preconditions:
        sampleFiles : Dict with 2 key values: fastq1, fastq2 containing file paths to input fastq.gz files
    '''
    USER = os.environ['USER']
    PI_SCRATCH = os.environ['PI_SCRATCH']

    PREPROCESSING_INPUT_FILE = '/home/groups/carilee/software/variant-discovery-pipeline/inputs/preprocessing_template.json'

    # Create folder to store cromwell logs
    cromwell_log_directory = os.path.join(PI_SCRATCH, "cromwell-monitor-logs")
    if os.path.isdir(cromwell_log_directory) is False:
        os.mkdir(cromwell_log_directory)

	# Set working directory to cromwell log directory - cromwell log subfolders will get created here
    os.chdir(cromwell_log_directory)
	
	# Modify input file
    inputPath = os.path.join(directory, 'preprocessing_{}.json'.format(sample))
    # Get Read group info
    ID, PU, PL = getReadGroupInfo(sampleFiles['fastq1'], sample)
    replacements = {'SAMPLE_NAME_HERE': sample, 'FASTQ1_HERE': sampleFiles['fastq1'], 
                    'FASTQ2_HERE':sampleFiles['fastq2'], 'OUTPUT_DIRECTORY_HERE' : directory,
                    'READ_GROUP_HERE' : ID, 'PLATFORM_UNIT_HERE': PU, 'PLATFORM_HERE': PL}

    with open(PREPROCESSING_INPUT_FILE) as infile, open(inputPath, 'w') as outfile:
        for line in infile:
            for src, target in replacements.items():
                line = line.replace(src, target)
            outfile.write(line)

    # print("Created inputs file")

    # sbatch directive variables
    jobName = 'cromwell_preprocessing'
    outputFile = '{}'.format(os.path.join(cromwell_log_directory, 'cromwell_preprocessing.%j.out'))
    errorFile = '{}'.format(os.path.join(cromwell_log_directory, 'cromwell_preprocessing.%j.err'))
    nodes = '1'
    memory = '32000' #in MB
    time = '0-06:00:00'
    email = '{}@stanford.edu'.format(USER)
    mailType = 'END'

    # command line execution variables
    softwareDirectory = '/home/groups/carilee/software'
    configPath = os.path.join(softwareDirectory, 'variant-discovery-pipeline', 'your.conf') # currently not using as SLURM isn't working
    jarPath = os.path.join(softwareDirectory, 'cromwell-35.jar')
    scriptPath = os.path.join(softwareDirectory, 'variant-discovery-pipeline', 'scripts', 'preprocessing.wdl')
    #wrap = 'java -jar {} run {} -i {}'.format(jarPath, scriptPath, inputPath)
    wrap = 'java -Dconfig.file={} -jar {} run {} -i {}'.format(configPath, jarPath, scriptPath, inputPath)

    # launch process
    sbatchCommand = 'sbatch --job-name={} --output={} --error={} --nodes={} --mem={} --time={} --mail-user={} --mail-type={} --wrap="{}"'.format(jobName, outputFile, errorFile, nodes, memory, time, email, mailType, wrap)

    stdout = subprocess.run(sbatchCommand, shell=True, stdout=subprocess.PIPE)
    print(stdout.stdout)

def launchJobs(samples, directory):
    for sample in samples:
        launchCromwellJob(sample, samples[sample], directory)

def launchBashJob(sample, sampleFiles, directory):
    '''
    Preconditions:
        sampleFiles : Dict with 2 key values: fastq1, fastq2 containing file paths to input fastq.gz files
    '''
    # Check if a directory to store all the customized sbatch scripts exists - if not make one
    # Create customized sbatch script for this sample
    # Submit job

def main():
    args = parseArgs()
    # print("Parsed args")
    samples = getSamples(args['directory'])
    # print("Parsed samples")
    launchJobs(samples, args['directory'])

if __name__ == '__main__':
    main()
