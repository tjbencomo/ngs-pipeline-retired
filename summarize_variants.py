'''Author: Tomas Bencomo
Summarizes variants for one or more variant call format files.
'''
import os
import argparse

import vcf as pyvcf
import pandas as pd

# VCF INFO FIELD keys from ANNOVAR
GENE_NAME = 'Gene.refGene'
MUTATION_TYPE = 'ExonicFunc.refGene'
VARIANT_TYPE = 'Func.refGene'

def build_position_summary(vcf_filepath, variant_summary=None):
    """Summarize variants by position.
    Returns variant summary as dict with position for keys
    """
    vcf_reader = pyvcf.Reader(open(vcf_filepath, 'r'))
    if variant_summary is None:
        variant_summary = {}
    for record in vcf_reader:
        variant_function = ''.join(record.INFO[VARIANT_TYPE])
        if variant_function in ('UTR3', 'UTR5', 'exonic'):
            position = '{}_{}'.format(record.CHROM, record.POS)
            variant_summary = update_position_entry(variant_summary, record,
                                                    position, variant_function)
    # No need to add missing entries as there are none
    return variant_summary

def update_position_entry(variant_summary, record, position, variant_type):
    """Add variant to the summary by position and return updated dict"""
    if position not in variant_summary:
        gene = ''.join(record.INFO[GENE_NAME])
        variant_summary[position] = {}
        variant_summary[position]['Gene'] = gene
        variant_summary[position]['Mutation'] = []
        variant_summary[position]['Function'] = []
        variant_summary[position]['Sample'] = []
    if variant_type in ('UTR3', 'UTR5'):
        variant_summary[position]['Function'].append(variant_type)
    elif variant_type == 'exonic':
        mut_type = ''.join(record.INFO[MUTATION_TYPE])
        variant_summary[position]['Function'].append(mut_type)
    mutation = '{}/{}'.format(record.REF, record.ALT[0])
    variant_summary[position]['Mutation'].append(mutation)
    sample_pair = []
    if len(record.samples) != 2:
        raise ValueError('Annotated VCF has more than 2 samples!')
    variant_summary[position]['Sample'].append('{}:{}'.format(
                                                    record.samples[0].sample,
                                                    record.samples[1].sample))
    return variant_summary

def summarize_by_position(vcfs):
    """Summarize list of vcfs by position and return summary dict"""
    variant_summary = {}
    for vcf in vcfs:
        variant_summary = build_position_summary(vcf, variant_summary)
    return variant_summary

def build_gene_summary(vcf_filepath, variant_summary=None):
    """Summarize variants by gene and returns variant summary as dict with genes for keys"""
    vcf_reader = pyvcf.Reader(open(vcf_filepath, 'r'))
    if variant_summary is None:
        variant_summary = {}
    for record in vcf_reader:
        variant_function = ''.join(record.INFO[VARIANT_TYPE])
        if variant_function in ('UTR3', 'UTR5', 'exonic'):
            gene = ''.join(record.INFO[GENE_NAME])
            variant_summary = update_gene_entry(variant_summary, record, gene, 
                                                variant_function)
    variant_summary = add_missing_gene_entries(variant_summary)
    return variant_summary

def update_gene_entry(variant_summary, record, gene, variant_type):
    """Add variant to the summary by gene and return updated dict"""
    if gene not in variant_summary:
        variant_summary[gene] = {}
        variant_summary[gene]['Position'] = []
        variant_summary[gene]['Mutation'] = []
        variant_summary[gene]['Sample'] = []
    if variant_type in ('UTR3', 'UTR5'): 
        if variant_type in variant_summary[gene]:
            variant_summary[gene][variant_type] += 1
        else:
            variant_summary[gene][variant_type] = 1
    elif variant_type == 'exonic':
        mutation_type = ''.join(record.INFO[MUTATION_TYPE])
        if mutation_type in variant_summary[gene]:
            variant_summary[gene][mutation_type] += 1
        else:
            variant_summary[gene][mutation_type] = 1
    variant_summary[gene]['Position'].append('{}_{}'.format(record.CHROM, 
                                                            record.POS))
    variant_summary[gene]['Mutation'].append('{}/{}'.format(record.REF,
                                                            record.ALT[0]))
    sample_pair = []
    if len(record.samples) != 2:
        raise ValueError('Annotated VCF has more than 2 samples!')
    variant_summary[gene]['Sample'].append('{}:{}'.format(
                                                    record.samples[0].sample,
                                                    record.samples[1].sample))
    return variant_summary
                                                            
def add_missing_gene_entries(variant_summary):
    """Add dict keys for any types of variants not found and set values to 0.
    Return updated dict with all needed keys
    """
    required_entries = ['UTR3', 'UTR5',
                        'frameshift_deletion','frameshift_insertion',
                        'nonframeshift_deletion','nonframeshift_insertion',
                        'nonframeshift_substitution','nonsynonymous_SNV','stopgain',
                        'stoploss','synonymous_SNV','unknown']
    for gene in variant_summary:
        for entry in required_entries:
            if entry not in variant_summary[gene]:
                variant_summary[gene][entry] = 0
    return variant_summary

def summarize_by_gene(vcfs):
    """Summarize list of vcf filenames by gene and return summary dict"""
    variant_summary = {}
    for vcf in vcfs:
        variant_summary = build_gene_summary(vcf, variant_summary)
    return variant_summary

def save_summary(variant_summary, filename):
    summary_table = build_summary_table(variant_summary)
    summary_table.to_csv(filename)

def build_summary_table(variant_summary):
    """Convert summary dict to pandas DataFrame - return DataFrame with genes as indices"""
    df = pd.DataFrame.from_dict(variant_summary, orient='index')
    return df

def test():
    data_directory = ('/scratch/groups/carilee/forTomas/CollagenFQData/alldata'
                        '/analysis-ready-bams')
    vcf = '119_annotated.hg19_multianno.vcf'
    vcf_filepath = os.path.join(data_directory, vcf)

    #summary = build_gene_summary(vcf_filepath)
    summary = build_position_summary(vcf_filepath)
    df = build_summary_table(summary)
    print(df.head())

def parseArgs():
    """Parse command line arguments and return args dict"""
    parser = arparse.ArgumentParser(description='Summarize annovar annotated vcfs')
    parser.add_argument('-I', '--input', type=str, nargs=1, dest='input_file')
    parser.add_argument('-d', '--directory', type=str, nargs=1, dest='directory')
    parser.add_argument('O', '--output', type=str, nargs=1, dest='output_file')
    args = parser.parse_args()

    if None in (args.input_file, args.output_file):
        raise ValueError('Missing input or output file! Check your parameters')

    return {'input_file' : ''.join(args.input_file),
            'output_file' : ''.join(args.output_file),
            'directory' : args.directory}

def read_input_file(input_file, directory=None):
    """Parses filenames from input file and returns as list of files"""
    files = []
    with open(input_file, 'r') as f:
        files = f.readlines()
        files = [line.rstrip('\n') for line in lines]
        if directory is not None:
            files = [os.path.join(directory, f) for f in files]
    return files

def normal_behavior():
    args = parseArgs()
    vcfs = read_input_file(args['input_file'], args['directory'])
   
    print('Building gene based summary file')
    gene_summary = summarize_by_gene(vcfs)
    gene_summary_filename = '{}_{}.csv'.format(args['output_file'], 
                                                'GENE')
    save_summary(gene_summary, gene_summary_filename)
    
    print('Building position based summary file')
    position_summary = summarize_by_position(vcfs)
    position_summary_filename = '{}_{}.csv'.format(args['output_file'], 
                                                    'COORDINATE')
    save_summary(position_summary, position_summary_filename)
    print('VCFs summarized!')

def main():
    test()

if __name__ == '__main__':
    main()


