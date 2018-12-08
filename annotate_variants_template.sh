#!/bin/bash
#SBATCH --job-name=annotate_variants
#SBATCH --output=/scratch/groups/carilee/cromwell-monitor-logs/annotate_variants.%j.out
#SBATCH --error=/scratch/groups/carilee/cromwell-monitor-logs/annotate_variants.%j.err
#SBATCH --nodes=1
#SBATCH --mem=4000
#SBATCH --cpus-per-task=1
#SBATCH --time=0-01:00:00
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=_USER_@stanford.edu 
#SBATCH --workdir=_WORKING_DIRECTORY_

ml system singularity

gatk_path=/home/groups/carilee/software/gatk4-latest.simg
ref_assembly=/home/groups/carilee/refs/hg19/ucsc.hg19.fasta
annotations_directory=/home/groups/carilee/refs/annovar-refs/humandb/

sample_name=_SAMPLE_NAME_
passed_only_vcf=${sample_name}_passed.vcf.gz

#this is the vcf straight from FilterMutectCalls from the Variant Calling Pipeline
input_vcf=_INPUT_VARIANT_FILE_
#prefix name for annovar to create its output files with
output_prefix=_OUTPUT_PREFIX_

#annovar assembly info
build_version=_BUILD_VERSION_

singularity exec -B $PI_SCRATCH ${gatk_path}  gatk SelectVariants -R ${ref_assembly} -V ${input_vcf} -O ${passed_only_vcf} --exclude-filtered

table_annovar.pl ${passed_only_vcf} ${annotations_directory} -buildver ${build_version} -out ${output_prefix} -remove -protocol refGene -operation g -nastring . -vcfinput
