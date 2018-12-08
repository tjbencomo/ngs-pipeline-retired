#!/bin/bash
#SBATCH --job-name=SNV_Calling
#SBATCH --output=/scratch/groups/carilee/cromwell-monitor-logs/SNV_Calling.%j.out
#SBATCH --error=/scratch/groups/carilee/cromwell-monitor-logs/SNV_Calling.%j.err
#SBATCH --nodes=1
#SBATCH --mem=_MEMORY_
#SBATCH --cpus-per-task=_CORES_
#SBATCH --time=0-16:00:00
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=_USER_@stanford.edu
#SBATCH --workdir=_WORKING_DIRECTORY_

module load singularity

sample_name=_SAMPLE_NAME_

tumor_bam=_TUMOR_BAM_
tumor_sample_name=_TUMOR_NAME_
normal_bam=_NORMAL_BAM_
normal_sample_name=_NORMAL_NAME_

reference_fasta=/home/groups/carilee/refs/hg19/ucsc.hg19.fasta
reference_fasta_index=/home/groups/carilee/refs/hg19/ucsc.hg19.fasta.fai
reference_dict=/home/groups/carilee/refs/hg19/ucsc.hg19.dict

panel_of_normals=_PANEL_OF_NORMALS_
germline_resource=/home/groups/carilee/refs/gatk-hg19-refs/af-only-gnomad.hg19.vcf.gz
contamination_germline_resource=/home/groups/carilee/refs/gatk-hg19-refs/small_exac_common_3_hg19.vcf.gz

output_vcf=${sample_name}_somatic_unfiltered.vcf.gz
output_bam=${sample_name}_tumor_normal_m2.bam

log_file=_LOG_FILE_

gatk_path=/home/groups/carilee/software/gatk4-latest.simg

# Stop script if any of the commands error out
set -e

echo "Running on ${SLURMD_NODENAME}"
echo "Variant Calling ${sample_name}"

# Mutect2 Calling to create .vcf file

singularity exec -B $PI_SCRATCH ${gatk_path} gatk Mutect2 \
    -R ${reference_fasta} \
    -I ${tumor_bam} \
    -I ${normal_bam} \
    -tumor ${tumor_sample_name} \
    -normal ${normal_sample_name} \
    -pon ${panel_of_normals} \
    --germline-resource ${germline_resource} \
	--af-of-alleles-not-in-resource 0.0000025 \
	--disable-read-filter MateOnSameContigOrNoMappedMateReadFilter \
    -O ${output_vcf} \
    -bamout ${output_bam}


# GetPileupSummaries to summarize read support
# NEED TO FIND -L File

pileup_summary=${sample_name}_pileupsummaries.table

singularity exec -B $PI_SCRATCH ${gatk_path} gatk  GetPileupSummaries \
    -I ${tumor_bam} \
    -V ${contamination_germline_resource} \
	-L ${contamination_germline_resource} \
    -O ${pileup_summary}

# CalculateContamination computes fraction contamination

contamination_table=${sample_name}_calculatecontamination.table

singularity exec -B $PI_SCRATCH ${gatk_path} gatk CalculateContamination \
    -I ${pileup_summary} \
    -O ${contamination_table}

# FilterMutectCalls to exclude calls that do not have a high
# degree of confidence

singularity exec -B $PI_SCRATCH ${gatk_path} gatk FilterMutectCalls \
    -V ${output_vcf} \
    --contamination-table ${contamination_table} \
    -O ${sample_name}_somatic_filtered.vcf.gz

echo "Completed ${sample_name}" >> ${log_file}
