workflow PreProcessingForVariantDiscovery {
    String sample_name
    String ref_name

    File input_fastq_1
    File input_fastq_2
    
    File ref_fasta
    File ref_fasta_index
    File ref_dict

    String bwa_commandline
    Int compression_level

    File dbSNP_vcf
    File dbSNP_vcf_index
    Array[File] known_indels_sites_VCFs
    Array[File] known_indels_sites_indices

    String gatk_path

    String base_file_name = sample_name + "." + ref_name
    String output_directory

    String read_group

    call GetBwaVersion {}

    call FastqToBam {
        input:
            fastq_1 = input_fastq_1,
            fastq_2 = input_fastq_2,
            sample_name = sample_name,
            read_group = read_group,
            output_file_name = base_file_name + ".unaligned.bam",
            gatk_path = gatk_path,
            output_directory = output_directory
    }

    call SamToFastqAndBwaMem {
        input:
            input_bam = FastqToBam.output_bam,
            bwa_commandline = bwa_commandline,
            output_bam_name = base_file_name + ".unmerged",
            output_directory = output_directory,
            ref_fasta = ref_fasta,
            ref_fasta_index = ref_fasta_index,
            ref_dict = ref_dict,
            gatk_path = gatk_path
    }
}

task GetBwaVersion {
    command <<<
        module load biology bwa
        bwa 2>&1 | \
        grep -e '^Version' | \
        sed 's/Version: //'
    >>>
    output {
        String version = read_string(stdout())
    }
}

task FastqToBam {
    File fastq_1
    File fastq_2
    
    String output_file_name
    String output_directory

    String sample_name
    String read_group

    String gatk_path

    command <<<
        module load system singularity
        singularity exec ${gatk_path} gatk FastqToSam \
            --FASTQ=${fastq_1} \
            --FASTQ2=${fastq_2} \
            --OUTPUT=${output_directory}${output_file_name} \
            --SAMPLE_NAME=${sample_name} \
            -RG="${read_group}"
    >>>
    output {
        File output_bam = output_directory + output_file_name
    }
}

task SamToFastqAndBwaMem {
    # Not sure what - hyphens do in bash script - assuming that it is used
    # to mark certain command line parameters as null
    File input_bam
    String bwa_commandline
    String output_bam_name
    String output_directory
    File ref_fasta
    File ref_fasta_index
    File ref_dict

    File ref_amb
    File ref_ann
    File ref_bwt
    File ref_pac
    File ref_sa

    String gatk_path

    command <<<
        set -o pipefail
        set -e

        module load system singularity
        module load biology bwa
        module load biology samtools

        singularity exec ${gatk_path} gatk SamToFastq \
            --INPUT=${input_bam} \
            --FASTQ=/dev/stdout \
            --INTERLEAVE=true \
            --INCLUDE_NON_PF_READS=true \
        | \
        bwa ${bwa_commandline} ${ref_fasta} /dev/stdin - 2> >(tee ${output_directory}${output_bam_name}.bwa.stderr.log >&2) \
        | \
        samtools view -1 - > ${output_directory}${output_bam_name}.bam        
    >>>
    output {
        File output_bam = "${output_directory}${output_bam_name}.bam"
        File bwa_stderr_log = "${output_directory}${output_bam_name}.bwa.stderr.log"
    }
}