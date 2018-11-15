# TASK DEFINITIONS
# NOTES:
# output_bam_basename is not just the filename. It is the entire
# path AND the basename without the extension (.whatever)
task findPairs {
  String directory
  File pairScript
  command <<<
    python ${pairScript} -I ${directory} -f fastq.gz
  >>>
  output {
    File pairings = stdout()
    
  }
}

task bwaAlign {
  # bwa and samtools must be already installed. On sherlock use module load bwa
  File ref_dict
  File ref_fasta
  File ref_fasta_index
  File ref_sa
  File ref_amb
  File ref_bwt
  File ref_ann
  File ref_pac
  
  Array[String] inputFiles
  String bwa_commandline
  String output_bam_basename = inputFiles[2]

  command <<<
    bwa ${bwa_commandline} ${ref_fasta} ${inputFiles[0]} ${inputFiles[1]} | samtools view -1 -T ${ref_fasta}
  >>>
  output {
    File outputBam = stdout()
    String outputBamName = output_bam_basename + '.bam'
    String output_file_basename = output_bam_basename
  }
}

task markDuplicates {
  File inputBam
  String output_bam_basename
  String metrics_filename
  String gatk_path

  command <<<
    ml system singularity
    singularity exec ${gatk_path} gatk MarkDuplicates \
      -I ${inputBam} -M ${metrics_filename} -O ${output_bam_basename}.bam \
	    --VALIDATION_STRINGENCY SILENT \
      --OPTICAL_DUPLICATE_PIXEL_DISTANCE 100 \ #use this as our sequencing data is from unpatterned machines - default value
	    --ASSUME_SORT_ORDER "queryname" \
      --CREATE_MD5_FILE true
  >>>
  output {
    File output_bam = "${output_bam_basename}.bam"
    File duplicate_metrics = "${metrics_filename}"
  }
}

task sortAndFixTags {
  # Using SetNmMdAndUqTags instead of SetNmAndUqTags
  # which was specified in the generic broad workflow
  # SetNmAndUqTags is depreciated and both serve the same purpose
  File input_bam
  String output_bam_basename
  File ref_dict
  File ref_fasta
  File ref_fasta_index
  String gatk_path

  command <<<
    set -o pipefail
    ml system singularity
    singularity exec ${gatk_path} gatk SortSam \
      --INPUT ${input_bam} \
      --OUTPUT /dev/stdout \
      --SORT_ORDER "coordinate" \
      --CREATE_INDEX false \
      --CREATE_MD5_FILE false \
      | \
    singularity exec ${gatk_path} gatk SetNmMdAndUqTags \
      --INPUT /dev/stdin \
      --OUTPUT ${output_bam_basename}.bam \
      --CREATE_INDEX true \
      --CREATE_MD5_FILE true \
      --REFERENCE_SEQUENCE ${ref_fasta}
  >>>
  output {
    File output_bam = "${output_bam_basename}.bam"
    File output_bam_index = "${output_bam_basename}.bai"
    File output_bam_md5 = "${output_bam_basename}.bam.md5"
  }
}

# WORKFLOW DEFINITION
workflow PreprocessingForVariantDiscovery {
  File ref_dict
  File ref_fasta
  File ref_fasta_index
  String gatk_path

  call findPairs
  Array[Array[String]] output_table = read_tsv(findPairs.pairings)
  
  scatter (pair in output_table) {
    call bwaAlign {
      input: 
        inputFiles = pair,
        ref_dict = ref_dict,
        ref_fasta = ref_fasta,
        ref_fasta_index = ref_fasta_index
    }
    String base_file_name = bwaAlign.output_file_basename
    call markDuplicates {
      input:
        inputBam = bwaAlign.outputBam,
        output_bam_basename = base_file_name + ".aligned.unsorted.duplicates_marked",
        metrics_filename = base_file_name + ".duplicate_metrics",
        gatk_path = gatk_path
    }
    call sortAndFixTags {
      input:
        input_bam = markDuplicates.output_bam,
        output_bam_basename = base_file_name + ".aligned.duplicate_marked.sorted",
        ref_dict = ref_dict,
        ref_fasta = ref_fasta,
        ref_fasta_index = ref_fasta_index,
        gatk_path = gatk_path
    }
  }  
}
