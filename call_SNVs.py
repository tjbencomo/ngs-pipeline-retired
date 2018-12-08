
'''
Author: Tomas Bencomo
Run multiple samples through the SNV Calling Pipeline
Inputs:
    -d directory
    -i input file
'''

import os
import argparse
import subprocess

def createParser():
    parser = argparse.ArgumentParser(description='Script to launch Cromwell Data Preprocessing Workflow')
    parser.add_argument('-d', '--directory', type=str, nargs=1, 
                        help='Specify directory to launch workflow from, saving execution files to this directory',
                        dest='directory')
    parser.add_argument('-i', '--inputFile', type=str, nargs=1,
                        help='File containing list of file names to process. All files must be from the same directory, 1 filename per line. Specify full path for inputFile',
                        dest='inputFile')

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

    if args.inputFile is None:
        inputFile = None
    else:
        inputFile = ''.join(args.inputFile)
        if os.path.isfile(inputFile) is False:
            raise ValueError("{} file does not exist!".format(inputFile))
    
    return {'directory' : directory, 'inputFile' : inputFile}

def getSamples(directory, files=None):
    if files is None:
        # Checks that file exists, is a .bam file, and only two .'s in name
        # Two .'s indicates it is the final bam file as analysis ready bam's
        # appear in the following format: sample_name.ref_assembly.bam
        files = [f for f in os.listdir(directory) if os.path.isfile(os.path.join(directory, f)) and f.endswith('.bam') and f.count('.') == 2]
    
    normalTag = 'CTR'
    tumorTag = 'T'

    samples = {}
    for f in files:
        file_id = f.split('.')[0]
        if normalTag in file_id:
            id = file_id[len(normalTag):]
            if id in samples:
                samples[id]['Normal'] = f
            else:
                samples[id] = {'Normal' : f}
        elif tumorTag in file_id:
            id = file_id[len(tumorTag):]
            if id in samples:
                samples[id]['Tumor'] = f
            else:
                samples[id] = {'Tumor' : f}

    return samples

def createCustomizedFile(template_file, customized_file, replacements):
    with open(template_file) as infile, open(customized_file, 'w') as outfile:
        for line in infile:
            for src, target in replacements.items():
                line = line.replace(src, target)
            outfile.write(line)

    return customized_file

def launchJob(sample, sampleFiles, directory, cores="1", memory="48000", normalTag='CTR', tumorTag='T'):
    # Check if a directory to store all the customized sbatch scripts exists - if not make one
    scriptsDirectory = os.path.join(directory, 'jobScripts')
    if os.path.isdir(scriptsDirectory) is False:
        os.mkdir(scriptsDirectory)
    # Create customized sbatch script for this sample
    TEMPLATE_SCRIPT = "/home/groups/carilee/software/variant-discovery-pipeline/SNV_calling_template.sh"
    USER = os.environ['USER']
    customizedScript = os.path.join(scriptsDirectory, 'SNV_calling_{}.sh'.format(sample))
    logFile = os.path.join(directory, 'SNV_calling_log.txt')
    
    '''This needs to be modified!'''
    #pon_filepath = '/scratch/groups/carilee/cromwell-test/mutect-tests/5WESpon.vcf.gz'
    pon_filepath = '/scratch/groups/carilee/forTomas/CollagenFQData/alldata/PanelofNormals/28WESpon.vcf.gz'
    replacements = {'_SAMPLE_NAME_': str(sample), '_CORES_': cores, '_MEMORY_': memory,
                    '_TUMOR_BAM_' : sampleFiles['Tumor'],
                    '_TUMOR_NAME_' : "{}{}".format(tumorTag, sample),
                    '_NORMAL_BAM_' : sampleFiles['Normal'],
                    '_NORMAL_NAME_' : "{}{}".format(normalTag, sample),
                    '_USER_': USER, '_LOG_FILE_': logFile,
					'_WORKING_DIRECTORY_': directory, '_PANEL_OF_NORMALS_': pon_filepath}
	
    createCustomizedFile(TEMPLATE_SCRIPT, customizedScript, replacements)
    stdout = subprocess.run('sbatch {}'.format(customizedScript), shell=True, stdout=subprocess.PIPE)
    print(stdout.stdout)

def launchJobs(samples, directory):
    print("DID YOU MODIFY THE PANEL OF NORMALS PATH") 
    for sample in samples:
        launchJob(sample, samples[sample], directory)

def read_input_file(filepath):
    files = []
    with open(filepath, 'r') as f:
        files = f.readlines()
        files = [f.rstrip('\n') for f in files]
    return files

def main():
    args = parseArgs()

    if args['inputFile'] is not None:
        files = read_input_file(args['inputFile'])
    else:
        files = None

    samples = getSamples(args['directory'], files)
    #print(samples)
    #print(len(samples))
    launchJobs(samples, args['directory'])

if __name__ == '__main__':
    main()
