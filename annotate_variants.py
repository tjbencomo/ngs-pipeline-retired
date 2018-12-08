'''
Author: Tomas Bencomo
Script to generate bash scripts that will annotate a given vcf file
Input:
-d : input directory containing the files to annotate
-i : file with selected files to annotate. Optional. If not given, uses all 
    .vcf files in directory
'''

import argparse
import os
import subprocess
import call_SNVs

def createParser():
    parser = argparse.ArgumentParser()
    parser.add_argument('-d', '--directory', type=str, nargs=1)
    parser.add_argument('-i', '--inputFile', type=str, nargs=1)

    args = parser.parse_args()
    return args

def parseArgs():
    args = createParser()

    if args.directory is None:
        directory = os.getcwd()
    else:
        directory = ''.join(args.directory)

    if os.path.isdir(directory) is False:
        raise ValueError("{} does not exist!".format(directory))

    if args.inputFile is None:
        inputFile = None
    else:
        inputFile = ''.join(args.inputFile)
        if os.path.isfile(inputFile) is False:
            raise ValueError("{} does not exist!".format(inputFile))

    return {'directory' : directory, 'inputFile' : inputFile}

def getSamples(directory, files=None):
    if files is None:
        files = [f for f in os.listdir(directory) if os.path.isfile(directory, f)
                and f.endswith('.vcf.gz')]
    
    samples = {}
    for f in files:
        name = f.split('_')[0]
        samples[name] = os.path.join(directory, f)

    return samples

def read_input_file(filepath):
    files = []
    with open(filepath, 'r') as f:
        files = f.readlines()
        files = [f.rstrip('\n') for f in files]
    return files

def createCustomizedFile(template_file, customized_file, replacements):
    with open(template_file) as infile, open(customized_file, 'w') as outfile:
        for line in infile:
            for src, target in replacements.items():
                line = line.replace(src, target)
            outfile.write(line)
    return customized_file

def launchJob(sample, inputVCF, directory, buildVersion='hg19'):
    scriptsDirectory = os.path.join(directory, 'jobScripts')
    if os.path.isdir(scriptsDirectory) is False:
        os.mkdir(scriptsDirectory)
    softwareDirectory = '/home/groups/carilee/software/variant-discovery-pipeline/'
    TEMPLATE_SCRIPT = os.path.join(softwareDirectory, 'annotate_variants_template.sh')

    USER = os.environ['USER']
    customizedScript = os.path.join(scriptsDirectory, 
                                    'annotate_variants_{}.sh'.format(sample))
    replacements = {'_USER_' : USER, '_WORKING_DIRECTORY_' : directory,
                    '_SAMPLE_NAME_' : sample, '_INPUT_VARIANT_FILE_' : inputVCF,
                    '_OUTPUT_PREFIX_' : "{}_annotated".format(sample),
                    '_BUILD_VERSION_' : buildVersion} 
    createCustomizedFile(TEMPLATE_SCRIPT, customizedScript, replacements)
    stdout = subprocess.run('sbatch {}'.format(customizedScript), shell=True,
                            stdout=subprocess.PIPE)
    print(stdout.stdout)

def launchJobs(samples, directory):
    for s in samples:
        launchJob(s, samples[s], directory)

def main():
    args = parseArgs()
    if args['inputFile'] is not None:
        files = read_input_file(args['inputFile'])
    else:
        files = None
    samples = getSamples(args['directory'], files)
    launchJobs(samples, args['directory']) 

if __name__ == '__main__':
    main()

