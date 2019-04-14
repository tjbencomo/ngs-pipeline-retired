'''
Wrapper script to automate preprocessing workflow submission.
User enters an input file containing sample information which
the script then processes and launches a preprocessing job for
each sample. The script also writes output information so that
the variant-calling script can pick up from there.

The input file is expected to have the following columns:
Sample, Type, FASTQ1, FASTQ2, OutputDirectory, SamplePrefix, Library
'''

import os
import json
import re
import pandas as pd

# 1 Read user input file
# Check each sample for requirements?
# Create pipeline dir
# Write jobs and submit

# TODO
# 1) Ask Michael about read groups and figure out how to read them
# 2) Make extract_read_info return dataframe with ID, PL, PU columns
# 3) Join these cols to df in main()
# 4) Create content directory - take as command line argument
# 5) JSON input generation
# 6) SBATCH Launch function either individual files or command line
# 7) Info to pick up with variant-calling

def read_input(filepath):
    if os.path.isfile(filepath):
        raise ValueError("{} does not exist!".format(filepath))
    if filepath.endswith('.xlsx'):
        df = pd.read_excel(filepath)
    else:
        df = pd.read_csv(filepath)
    return df

def check_requirements(df):
    '''
    Currently just checking that no missing data in input file
    Future should check READ TAGS in FASTQ files
    '''
    return df.isnull().sum() == 0

def extract_read_info(fastq_files):
    read_tags = []
    for filepath in fastq_files:
        with open(filepath, 'r') as f:
            '''@<instrument>:<run number>:<flowcell ID>:<lane>:<tile>:<x-pos>:<y-pos>'''
            x = f.read_line()
            x = re.split('@|:', x)
            read_tags.append(x)
    return read_tags

def 

def main():
    # get args
    filepath = 'command_args_here'
    df = read_input(filepath)
    if check_requirements(df) is False:
        raise ValueError("Input file incorrect, check to make sure nothing is missing!")
    

    pass

if __name__ == '__main__':
    main()

