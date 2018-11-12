# TASK DEFINITIONS
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
      --OPTICAL_DUPLICATE_PIXEL_DISTANCE 2500 \
	    --ASSUME_SORT_ORDER "queryname" \
      --CREATE_MD5_FILE true
  >>>
  output {
    File output_bam = "${output_bam_basename}.bam"
    File duplicate_metrics = "${metrics_filename}"
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
        inputFiles=pair,
        ref_dict=ref_dict,
        ref_fasta=ref_fasta,
        ref_fasta_index=ref_fasta_index
    }
    call markDuplicates {
      input:
        inputBam=bwaAlign.outputBam,
        output_bam_basename=bwaAlign.output_file_basename + ".aligned.unsorted.duplicates_marked",
        metrics_filename=bwaAlign.output_file_basename + ".duplicate_metrics",
        gatk_path=gatk_path
    }
  }  
}
