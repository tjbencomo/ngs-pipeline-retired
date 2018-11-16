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

    call MergeBamAlignment {
        input:
            unaligned_bam = FastqToBam.output_bam,
            bwa_commandline = bwa_commandline,
            bwa_version = GetBwaVersion.version,
            aligned_bam = SamToFastqAndBwaMem.output_bam,
            output_bam_name = base_file_name + ".aligned.unsorted",
            output_directory = output_directory,
            ref_fasta = ref_fasta,
            ref_fasta_index = ref_fasta_index,
            ref_dict = ref_dict,
            gatk_path = gatk_path
    }

    call MarkDuplicates {
        input:
            input_bam = MergeBamAlignment.output_bam,
            output_bam_name = base_file_name + ".aligned.unsorted.duplicates_marked",
            metrics_filename = base_file_name + ".duplicate_metrics",
            output_directory = output_directory,
            gatk_path = gatk_path
    }

    call SortAndFixTags {
        input:
            input_bam = MarkDuplicates.output_bam,
            output_bam_name = base_file_name + ".aligned.duplicate_marked.sorted",
            output_directory = output_directory,
            ref_dict = ref_dict,
            ref_fasta = ref_fasta,
            ref_fasta_index = ref_fasta_index,
            gatk_path = gatk_path
    }

    call CreateSequenceGroupingTSV {
        input:
            # this is to avoid having to use the gzip library
            # and modify the GATK workflow code to handle byte data
            ref_dict = "/home/groups/carilee/refs/hg19/ucsc.hg19.dict"
    }

    scatter (subgroup in CreateSequenceGroupingTSV.sequence_grouping) {
        call BaseRecalibrator {
            input:
                input_bam = SortAndFixTags.output_bam,
                input_bam_index = SortAndFixTags.output_bam_index,
                recalibration_report_filename = base_file_name + ".recal_data.csv",
                sequence_group_interval = subgroup,
                dbSNP_vcf = dbSNP_vcf,
                dbSNP_vcf_index = dbSNP_vcf_index,
                known_indels_sites_VCFs = known_indels_sites_VCFs,
                known_indels_sites_indices = known_indels_sites_indices,
                ref_dict = ref_dict,
                ref_fasta = ref_fasta,
                ref_fasta_index = ref_fasta_index,
                gatk_path = gatk_path,
        }
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

task MergeBamAlignment {
    File unaligned_bam
    String bwa_commandline
    String bwa_version
    File aligned_bam
    String output_bam_name
    String output_directory
    File ref_fasta
    File ref_fasta_index
    File ref_dict
    # Added in as it seems MergeBamAlignment can't identify the .gz version of dict
    # Look into this in the future
    File ref_dict2 = "/home/groups/carilee/refs/hg19/ucsc.hg19.dict"

    String gatk_path

    command <<<
        module load system singularity
        singularity exec ${gatk_path} gatk MergeBamAlignment \
            --VALIDATION_STRINGENCY SILENT \
            --EXPECTED_ORIENTATIONS FR \
            --ATTRIBUTES_TO_RETAIN X0 \
            --ALIGNED_BAM ${aligned_bam} \
            --UNMAPPED_BAM ${unaligned_bam} \
            --OUTPUT ${output_directory}${output_bam_name}.bam \
            --REFERENCE_SEQUENCE ${ref_fasta} \
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
            --PROGRAM_GROUP_VERSION "${bwa_version}" \
            --PROGRAM_GROUP_COMMAND_LINE "${bwa_commandline}" \
            --PROGRAM_GROUP_NAME "bwamem" \
            --UNMAPPED_READ_STRATEGY COPY_TO_TAG \
            --ALIGNER_PROPER_PAIR_FLAGS true \
            --UNMAP_CONTAMINANT_READS true
    >>>
    output {
        File output_bam = "${output_directory}${output_bam_name}.bam"
    }
}

task MarkDuplicates {
    File input_bam
    String output_bam_name
    String metrics_filename
    String output_directory

    String gatk_path

    command <<<
        module load system singularity
        singularity exec ${gatk_path} gatk MarkDuplicates \
            --INPUT ${input_bam} \
            --OUTPUT ${output_directory}${output_bam_name}.bam \
            --METRICS_FILE ${output_directory}${metrics_filename} \
            --VALIDATION_STRINGENCY SILENT \
            --OPTICAL_DUPLICATE_PIXEL_DISTANCE 100 \
            --ASSUME_SORT_ORDER "queryname" \
            --CREATE_MD5_FILE true
    >>>
    output {
        File output_bam = "${output_directory}${output_bam_name}.bam"
        File duplicate_metrics = "${output_directory}${metrics_filename}"
    }
}

task SortAndFixTags {
    File input_bam
    String output_bam_name
    String output_directory
    File ref_dict
    File ref_fasta
    File ref_fasta_index
    String gatk_path
    
    command <<<
        set -o pipefail
        module load singularity
        singularity exec ${gatk_path} gatk SortSam \
            --INPUT ${input_bam} \
            --OUTPUT /dev/stdout \
            --SORT_ORDER "coordinate" \
            --CREATE_INDEX false \
            --CREATE_MD5_FILE false \
        | \
        singularity exec ${gatk_path} gatk SetNmMdAndUqTags \
            --INPUT /dev/stdin \
            --OUTPUT ${output_directory}${output_bam_name}.bam \
            --CREATE_INDEX true \
            --CREATE_MD5_FILE true \
            --REFERENCE_SEQUENCE ${ref_fasta}
    >>>
    output {
        File output_bam = "${output_directory}${output_bam_name}.bam"
        File output_bam_index = "${output_directory}${output_bam_name}.bai"
        File output_bam_md5 = "${output_directory}${output_bam_name}.bam.md5"
    }
}

task CreateSequenceGroupingTSV {
    File ref_dict

    command <<<
        python <<CODE
        with open("${ref_dict}", "r") as ref_dict_file:
            sequence_tuple_list = []
            longest_sequence = 0
            for line in ref_dict_file:
                if line.startswith("@SQ"):
                    line_split = line.split("\t")
                    # (Sequence_Name, Sequence_Length)
                    sequence_tuple_list.append((line_split[1].split("SN:")[1], int(line_split[2].split("LN:")[1])))
            longest_sequence = sorted(sequence_tuple_list, key=lambda x: x[1], reverse=True)[0][1]
            # We are adding this to the intervals because hg38 has contigs named with embedded colons (:) and a bug in 
            # some versions of GATK strips off the last element after a colon, so we add this as a sacrificial element.
            hg38_protection_tag = ":1+"
            # initialize the tsv string with the first sequence
            tsv_string = sequence_tuple_list[0][0] + hg38_protection_tag
            temp_size = sequence_tuple_list[0][1]
            for sequence_tuple in sequence_tuple_list[1:]:
                if temp_size + sequence_tuple[1] <= longest_sequence:
                    temp_size += sequence_tuple[1]
                    tsv_string += "\t" + sequence_tuple[0] + hg38_protection_tag
                else:
                    tsv_string += "\n" + sequence_tuple[0] + hg38_protection_tag
                    temp_size = sequence_tuple[1]
            # add the unmapped sequences as a separate line to ensure that they are recalibrated as well
            with open("sequence_grouping.txt","w") as tsv_file:
                tsv_file.write(tsv_string)
                tsv_file.close()

            tsv_string += '\n' + "unmapped"

            with open("sequence_grouping_with_unmapped.txt","w") as tsv_file_with_unmapped:
                tsv_file_with_unmapped.write(tsv_string)
                tsv_file_with_unmapped.close()
        CODE
    >>>
    output {
        Array[Array[String]] sequence_grouping = read_tsv("sequence_grouping.txt")
        Array[Array[String]] sequence_grouping_with_unmapped = read_tsv("sequence_grouping_with_unmapped.txt")
    }
}

task BaseRecalibrator {
    File input_bam
    File input_bam_index
    String recalibration_report_filename
    Array[String] sequence_group_interval
    File dbSNP_vcf
    File dbSNP_vcf_index
    Array[File] known_indels_sites_VCFs
    Array[File] known_indels_sites_indices
    File ref_dict
    File ref_fasta
    File ref_fasta_index

    String gatk_path

    command <<<
        module load singularity
        singularity exec ${gatk_path} gatk BaseRecalibrator \
            -R ${ref_fasta} \
            -I ${input_bam} \
            --use-original-qualities \
            -O ${recalibration_report_filename} \
            --known-sites ${dbSNP_vcf} \
            --known-sites ${sep=" --known-sites " known_indels_sites_VCFs} \
            -L ${sep=" -L " sequence_group_interval}
    >>>
    output {
        File recalibration_report = "${recalibration_report_filename}"
    }
}




# TODO
# 2) Do BaseRecalibrator Task
# 3) Do ApplyBQSR Task
# 4) Look into runtime parameters for tasks
# 5) Modify data_processing.py for individual samples
# 6) Should dict not be the .gz version? - Consider just unzipping all the ref files
# 7) Figure out how SAM/BAMs work and the what the dict is for
# 8) Document Code and Put README on Github