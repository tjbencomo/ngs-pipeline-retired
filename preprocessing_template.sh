#!/bin/bash
#SBATCH --job-name=preprocessing
#SBATCH --output=/scratch/groups/carilee/cromwell-monitor-logs/preprocessing.%j.out
#SBATCH --error=/scratch/groups/carilee/cromwell-monitor-logs/preprocessing.%j.err
#SBATCH --nodes=1
#SBATCH --mem=_MEMORY_
#SBATCH --cpus-per-task=_CORES_
#SBATCH --time=0-08:00:00
#SBATCH --mail-type=END
#SBATCH --mail-user=_USER_@stanford.edu
#SBATCH --workdir=_WORKING_DIRECTORY_

# Make sure script stops if any of the commands or piping fails

set -o pipefail
set -e

module load singularity
module load biology samtools
module load biology bwa

REF_FASTA=/home/groups/carilee/refs/hg19/ucsc.hg19.fasta
REF_FASTA_INDEX=/home/groups/carilee/refs/hg19/ucsc.hg19.fasta.fai
REF_DICT=/home/groups/carilee/refs/hg19/ucsc.hg19.dict

DB_SNP=/home/groups/carilee/refs/hg19/dbsnp_138.hg19.vcf
DB_SNP_INDEX=/home/groups/carilee/refs/hg19/dbsnp_138.hg19.vcf.idx

KNOWN_INDEL_SITES=("/home/groups/carilee/refs/hg19/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf" "/home/groups/carilee/refs/hg19/1000G_phase1.indels.hg19.sites.vcf")

SAMPLE_NAME=_SAMPLE_NAME_

REF_NAME=_REFERENCE_NAME_

LOG_FILE=_LOG_FILE_

FASTQ1=_FASTQ1_FILE_
FASTQ2=_FASTQ2_FILE_

READ_GROUP=_READ_GROUP_
PLATFORM_UNIT=_PLATFORM_UNIT_
PLATFORM=_PLATFORM_

BWA_VERSION=0.7.17-r1188

GATK_PATH=/home/groups/carilee/software/gatk4-latest.simg

BASE_FILE_NAME=${SAMPLE_NAME}.${REF_NAME}

echo "Preprocessing ${SAMPLE_NAME}"

# FastqToSam Convert to Unaligned BAM

UNALIGNED_BAM=${BASE_FILE_NAME}.unaligned.bam

singularity exec ${GATK_PATH} gatk FastqToSam \
            --FASTQ=${FASTQ1} \
            --FASTQ2=${FASTQ2} \
            --OUTPUT=${UNALIGNED_BAM} \
            --PLATFORM_UNIT=${PLATFORM_UNIT} \
            --PLATFORM=${PLATFORM} \
            --SAMPLE_NAME=${SAMPLE_NAME} \
            -RG=${READ_GROUP}

# SamToFastqAndBwaMem Align reads and save to BAM file

BWA_COMMANDLINE="mem -K 100000000 -p -v 3 -t ${SLURM_CPUS_ON_NODE} -Y"

UNMERGED_BAM=${BASE_FILE_NAME}.unmerged.bam

singularity exec ${GATK_PATH} gatk SamToFastq \
            --INPUT=${UNALIGNED_BAM} \
            --FASTQ=/dev/stdout \
            --INTERLEAVE=true \
            --INCLUDE_NON_PF_READS=true \
    | \
    bwa ${BWA_COMMANDLINE} ${REF_FASTA} /dev/stdin - 2> >(tee ${UNMERGED_BAM}.bwa.stderr.log >&2) \
    | \
    samtools view -1 - > ${UNMERGED_BAM}

# MergeBamAlignment Combine aligned and unaligned bams to retain all metadata

ALIGNED_BAM=${BASE_FILE_NAME}.aligned.unsorted.bam

singularity exec ${GATK_PATH} gatk MergeBamAlignment \
            --VALIDATION_STRINGENCY SILENT \
            --EXPECTED_ORIENTATIONS FR \
            --ATTRIBUTES_TO_RETAIN X0 \
            --ALIGNED_BAM ${UNMERGED_BAM} \
            --UNMAPPED_BAM ${UNALIGNED_BAM} \
            --OUTPUT ${ALIGNED_BAM} \
            --REFERENCE_SEQUENCE ${REF_FASTA} \
            --PAIRED_RUN true \
            --SORT_ORDER "unsorted" \
            --IS_BISULFITE_SEQUENCE false \
            --ALIGNED_READS_ONLY false \
            --CLIP_ADAPTERS false \
            --MAX_RECORDS_IN_RAM 2000000 \
            --ADD_MATE_CIGAR true \
            --MAX_INSERTIONS_OR_DELETIONS -1 \
            --PRIMARY_ALIGNMENT_STRATEGY MostDistant \
            --PROGRAM_RECORD_ID "bwamem" \
            --PROGRAM_GROUP_VERSION "${BWA_VERSION}" \
            --PROGRAM_GROUP_COMMAND_LINE "${BWA_COMMANDLINE}" \
            --PROGRAM_GROUP_NAME "bwamem" \
            --UNMAPPED_READ_STRATEGY COPY_TO_TAG \
            --ALIGNER_PROPER_PAIR_FLAGS true \
            --UNMAP_CONTAMINANT_READS true

# MarkDuplicates Identify duplicate reads within bam file

MARKED_BAM=${BASE_FILE_NAME}.aligned.unsorted.duplicates_marked.bam
METRICS_FILE=${BASE_FILE_NAME}.duplicate_metrics

singularity exec ${GATK_PATH} gatk MarkDuplicates \
            --INPUT ${ALIGNED_BAM} \
            --OUTPUT ${MARKED_BAM} \
            --METRICS_FILE ${METRICS_FILE} \
            --VALIDATION_STRINGENCY SILENT \
            --OPTICAL_DUPLICATE_PIXEL_DISTANCE 100 \
            --ASSUME_SORT_ORDER "queryname" \
            --CREATE_MD5_FILE true

# SortAndFixTags Sort Bam file by coordinate and correct tag errors

SORTED_BAM=${BASE_FILE_NAME}.aligned.duplicate_marked.sorted.bam

singularity exec ${GATK_PATH} gatk SortSam \
            --INPUT ${MARKED_BAM} \
            --OUTPUT /dev/stdout \
            --SORT_ORDER "coordinate" \
            --CREATE_INDEX false \
            --CREATE_MD5_FILE false \
| \
singularity exec ${GATK_PATH} gatk SetNmMdAndUqTags \
    --INPUT /dev/stdin \
    --OUTPUT ${SORTED_BAM} \
    --CREATE_INDEX true \
    --CREATE_MD5_FILE true \
    --REFERENCE_SEQUENCE ${REF_FASTA}

# BaseRecalibrator Recompute base quality scores according to technical replicate info
# Cromwell workflow splits this up among several worker tasks - for simplicity
# this is not parallelized

RECALIBRATION_REPORT_FILE=${BASE_FILE_NAME}.recal_data.csv
MILLS_INDELS=/home/groups/carilee/refs/hg19/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf
PHASE_1_INDELS=/home/groups/carilee/refs/hg19/1000G_phase1.indels.hg19.sites.vcf

singularity exec  ${GATK_PATH} gatk BaseRecalibrator \
    -R ${REF_FASTA} \
    -I ${SORTED_BAM} \
    --use-original-qualities \
    -O ${RECALIBRATION_REPORT_FILE} \
    --known-sites ${DB_SNP} \
    --known-sites ${MILLS_INDELS} \
    --known-sites ${PHASE_1_INDELS}

# ApplyBQSR Replace base quality scores with recomputed scores from BaseRecalibrator

ANALYSIS_READY_BAM=${BASE_FILE_NAME}.bam

singularity exec ${GATK_PATH} gatk ApplyBQSR \
    -R ${REF_FASTA} \
    -I ${SORTED_BAM} \
    -O ${ANALYSIS_READY_BAM} \
    -bqsr ${RECALIBRATION_REPORT_FILE} \
    --static-quantized-quals 10 --static-quantized-quals 20 --static-quantized-quals 30 \
    --add-output-sam-program-record \
    --create-output-bam-md5 \
    --use-original-qualities

echo "${SAMPLE_NAME} completed" >> ${LOG_FILE}
