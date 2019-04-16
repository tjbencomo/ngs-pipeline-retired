'''
Wrapper script to automate preprocessing workflow submission.
User enters an input file containing sample information which
the script then processes and launches a preprocessing job for
each sample. The script also writes output information so that
the variant-calling script can pick up from there.

The input file is expected to have the following columns:
Sample, Type, FASTQ1, FASTQ2, OutputDirectory, SamplePrefix, Library

Read group info:
ID: FLOWCELL.LANE
PL: ILLUMINA
LB: Taken from CSV

Reference file format:
Reference, Location
'''

import os
import json
import re
import argparse
import csv
import gzip
import subprocess
import pandas as pd

# 1 Read user input file
# Check each sample for requirements?
# Create pipeline dir
# Write jobs and submit

# TODO
# 1) Ask Michael about read groups and figure out how to read them
# 2) Make extract_read_info return dataframe with ID, PL, PU columns - Done
# 3) Join these cols to df in main() - Done
# 4) Create content directory - take as command line argument - Done
# 5) JSON input generation - Done
# 6) SBATCH Launch function either individual files or command line
# 7) Info to pick up with variant-calling

def parse_args():
    parser = argparse.ArgumentParser(description = 'Initiate pipeline')
    parser.add_argument('samples', metavar='sample', type=str, nargs=1,
                            help='tabular file with sample info')
    parser.add_argument('ref_files', metavar='refs', type=str, nargs=1,
                            help='csv file with reference filepaths')
    parser.add_argument('pipeline_dir', metavar='workdir', type=str, nargs=1,
                            help='where the pipeline will save execution logs')
    args = parser.parse_args()
    return {'input': args.samples[0], 
            'refs' : args.ref_files[0], 
            'workdir' : args.pipeline_dir[0] }

def verify_workdir(workdir):
    if os.path.isdir(workdir) is False:
        print("{} does not exist. Creating new directory for logs...".format(workdir))
        try:
            os.mkdir(workdir)
            print("workdir created!")
        except OSError:
            print("Creation of {} workdir failed. Input an existing folder".format(workdir))

def read_input(filepath):
    if os.path.isfile(filepath) is False:
        raise ValueError("{} does not exist!".format(filepath))
    if filepath.endswith('.xlsx'):
        df = pd.read_excel(filepath)
    else:
        df = pd.read_csv(filepath)
    
    df.loc[(df['SamplePrefix'].isnull()) & (df['Type'] == 'Tumor'), 'SamplePrefix'] = 'T'
    df.loc[(df['SamplePrefix'].isnull()) & (df['Type'] == 'Normal'), 'SamplePrefix'] = 'CTR'
    df['Sample'] = df['Sample'].astype(str)
    return df

def check_refs(refs):
    '''
    Ensure all the required refs are in the refs dict
    '''
    keys = ['known_indels_sites_indices', 'dbSNP_vcf', 'dbSNP_vcf_index',
            'ref_dict', 'ref_fasta_index', 'known_indels_sites_VCFs', 
            'ref_name', 'ref_fasta', 'ref_sa', 'ref_amb', 'ref_bwt', 
            'ref_ann', 'ref_pac', 'cromwell_jar', 'email']
    missing = set(keys) - set(refs.keys())
    if len(missing) > 0:
        print(missing)
    return set(refs.keys()) == set(keys)

def read_references(filepath):
    if os.path.isfile(filepath) is False:
        raise ValueError("{} does not exist!".format(filepath))
    refs = {}
    with open(filepath, newline='') as f:
        f.readline() # chew through colnames
        reader = csv.reader(f)
        for row in reader:
            if row[0] not in refs:
                refs[row[0]] = row[1]
            elif row[0] in refs:
                if refs[row[0]] is list:
                    refs[row[0]].append(row[1])
                else:
                    temp = refs[row[0]]
                    refs[row[0]] = []
                    refs[row[0]].append(temp)
                    refs[row[0]].append(row[1])
    if check_refs(refs):
        return refs
    else:
        raise ValueError("Reference file didn't contain all the refs. Check documentation")

def check_requirements(df):
    '''
    Currently just checking that no missing data in input file
    Future should check READ TAGS in FASTQ files
    '''
    return df.isnull().sum() == 0

def read_ID_fastq(filepath):
    with open(filepath, 'r') as f:
        '''@<flowcell ID>:<lane>:<tile>:<x-pos>:<y-pos>'''
        line = f.readline().rstrip('\n')
        line = re.split(':|#', line)
        flowcell = line[0].strip('@')
        lane = line[1]
        return '{}.{}'.format(flowcell, lane)
    raise ValueError("ID Extraction failed for {}".format(filepath))

def read_ID_fastqgz(filepath):
    with gzip.open(filepath, 'rb') as f:
        '''@<flowcell ID>:<lane>:<tile>:<x-pos>:<y-pos>'''
        line = f.readline().decode('utf-8').rstrip('\n')
        line = re.split(':|#', line)
        flowcell = line[0].strip('@')
        lane = line[1]
        return '{}.{}'.format(flowcell, lane)
    raise ValueError("ID Extraction failed for {}".format(filepath))

def extract_read_info(fastq_files):
    read_IDs = []
    for filepath in fastq_files:
        if filepath.endswith('.gz'):
            ID = read_ID_fastqgz(filepath)
        else:
            ID = read_ID_fastq(filepath)
        read_IDs.append(ID)
    return read_IDs

def create_json(sample, fastq1, fastq2, output_dir, read_group, platform, refs):
    bwa = "mem -K 100000000 -p -v 3 -t $SLURM_CPUS_ON_NODE -Y"
    compression = 5
    inputs = {
            "PreProcessingForVariantDiscovery.sample_name" : sample,
            "PreProcessingForVariantDiscovery.input_fastq_1" : fastq1,
            "PreProcessingForVariantDiscovery.input_fastq_2" : fastq2,
            "PreProcessingForVariantDiscovery.output_directory" : output_dir,
            "PreProcessingForVariantDiscovery.read_group" : read_group,
            "PreProcessingForVariantDiscovery.platform" : platform,
            "PreProcessingForVariantDiscovery.known_indels_sites_indices" : refs['known_indels_sites_indices'],
            "PreProcessingForVariantDiscovery.compression_level" : compression,
            "PreProcessingForVariantDiscovery.dbSNP_vcf": refs['dbSNP_vcf'],
            "PreProcessingForVariantDiscovery.dbSNP_vcf_index" : refs['dbSNP_vcf_index'],
            "PreProcessingForVariantDiscovery.ref_dict" : refs['ref_dict'],
            "PreProcessingForVariantDiscovery.ref_fasta_index" : refs['ref_fasta_index'],
            "PreProcessingForVariantDiscovery.known_indels_sites_VCFs" : refs['known_indels_sites_VCFs'],
            "PreProcessingForVariantDiscovery.ref_name" : refs['ref_name'],
            "PreProcessingForVariantDiscovery.bwa_commandline" : bwa,
            "PreProcessingForVariantDiscovery.ref_fasta" : refs['ref_fasta'],
            "PreProcessingForVariantDiscovery.SamToFastqAndBwaMem.ref_sa" : refs['ref_sa'],
            "PreProcessingForVariantDiscovery.SamToFastqAndBwaMem.ref_amb" : refs['ref_amb'],
            "PreProcessingForVariantDiscovery.SamToFastqAndBwaMem.ref_bwt" : refs['ref_bwt'],
            "PreProcessingForVariantDiscovery.SamToFastqAndBwaMem.ref_ann" : refs['ref_ann'],
            "PreProcessingForVariantDiscovery.SamToFastqAndBwaMem.ref_pac" : refs['ref_pac']
            }
    filepath = os.path.join(output_dir, '{}_preprocess.json'.format(sample))
    if os.path.exists(filepath):
        raise ValueError("Already existing file at {}".format(filepath))
    with open(filepath, 'w') as f:
        json.dump(inputs, f)
    return filepath

def create_job(sample, inputs_filepath, cromwell_path, workdir, output_dir, email):
    job_filepath = os.path.join(output_dir, '{}_preprocess.sh'.format(sample))
    with open(job_filepath, 'w') as f:
        f.write("#!/bin/bash\n")
        f.write('#SBATCH --job-name=variant-discovery-pipeline\n')
        f.write('#SBATCH --output={}_preprocess.out\n'.format(os.path.join(output_dir, sample)))
        f.write('#SBATCH --error={}_preprocess.err\n'.format(os.path.join(output_dir, sample)))
        f.write('#SBATCH --nodes=1\n')
        f.write('#SBATCH --mem=4000\n')
        f.write('#SBATCH --cpus-per-task=1\n')
        f.write('#SBATCH --time=1-12:00:00\n')
        f.write('#SBATCH --mail-type=FAIL\n')
        f.write('#SBATCH --mail-user={}\n'.format(email))
        f.write('#SBATCH --workdir={}\n'.format(workdir))
        f.write('ml java\n')
        f.write('set -e\n')
        pipe_dir= os.path.dirname(os.path.realpath(__file__))
        config_filepath = os.path.join(pipe_dir, 'your.conf')
        preproc_wdl= os.path.join(pipe_dir, 'scripts', 'preprocessing.wdl')
        options_filepath = os.path.join(pipe_dir, 'options.json')
        f.write('java -Dconfig.file={} -jar {} run {} -i {} -o {}'.format(config_filepath,
                                                                            cromwell_path,
                                                                            preproc_wdl,
                                                                            inputs_filepath,
                                                                            options_filepath))
        return job_filepath

# def launch_pipeline(row, refs, workdir):
def launch_pipeline(row, refs, workdir):
    # Create json input file and save
    # print(row)
    # print()
    sample_name = row['SamplePrefix'] + row['Sample']
    print(sample_name)
    read_group = row['ID']
    platform = 'ILLUMINA'
    json_filepath  = create_json(sample_name, row['FASTQ1'], row['FASTQ2'], 
                                    row['OutputDirectory'], read_group, platform, refs)
    # # Run sbatch
    job_filepath = create_job(sample_name, json_filepath, refs['cromwell_jar'],
                                workdir, row['OutputDirectory'], refs['email'])
    # print("Submitting job script for {}".format(sample_name))
    # subprocess.run(['sbatch', job_filepath])
    # print("Job submitted for {}".format(sample_name))

def test(row):
    print(row)

def main():
    args = parse_args()
    verify_workdir(args['workdir'])
    df = read_input(args['input'])
    refs = read_references(args['refs'])
    if check_requirements(df) is False:
        raise ValueError("Input file incorrect, check to make sure nothing is missing!")
    df['ID'] = extract_read_info(df['FASTQ1'])
    df.apply(lambda x: launch_pipeline(x, refs, args['workdir']), axis=1)

if __name__ == '__main__':
    main()

