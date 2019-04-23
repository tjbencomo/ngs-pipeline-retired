workflow VariantCalling {
  String sample_name
  String normal_sample_name
  String tumor_sample_name
  
  File tumor_bam
  File normal_bam

  File ref_fasta
  File ref_fasta_index
  File ref_dict

  File panel_of_normals
  File germline_resource
  File contamination_germline_resource

  String output_directory

   call Mutect2 {
    input:
      ref_fasta = ref_fasta,
      tumor_bam = tumor_bam,
      normal_bam = normal_bam,
      panel_of_normals = panel_of_normals,
      germline_resource = germline_resource,
      tumor_sample_name = tumor_sample_name,
      normal_sample_name = normal_sample_name,
      sample_name = sample_name
  }

  call GetPileupSummaries {
    input:
      tumor_bam = tumor_bam,
      contamination_germline_resource = contamination_germline_resource,
      sample_name = sample_name
  }

  call CalculateContamination {
    input:
      pileup_summary = GetPileupSummaries.pileup_summary,
      sample_name = sample_name

  }

  call FilterMutectCalls {
    input:
      unfiltered_vcf = Mutect2.vcf,
      contamination_table = CalculateContamination.contamination_table,
      sample_name = sample_name
  }

  call CopyOutputs {
    input:
      output_directory = output_directory,
      filtered_vcf = FilterMutectCalls.filtered_vcf,
      output_bam = Mutect2.bam,
      contamination_table = CalculateContamination.contamination_table,
      pileup_summary = GetPileupSummaries.pileup_summary
  }
}

task Mutect2 {
    File ref_fasta
    File tumor_bam
    File normal_bam
    File panel_of_normals
    File germline_resource
    String tumor_sample_name
    String normal_sample_name
    String sample_name

    command {
      set -e
      ml biology gatk
      gatk Mutect2 \
        -R ${ref_fasta} \
        -I ${tumor_bam} \
        -I ${normal_bam} \
        -tumor ${tumor_sample_name} \
        -normal ${normal_sample_name} \
        -pon ${panel_of_normals} \
        --germline-resource ${germline_resource} \
        --af-of-alleles-not-in-resource 0.0000025 \
        --disable-read-filter MateOnSameContigOrNoMappedMateReadFilter \
        -O ${sample_name}_somatic_unfiltered.vcf.gz \
        -bamout ${sample_name}_tumor_normal_m2.bam
    }
    runtime {
      runtime_minutes: "180"
      memory: 32000
    }
    output {
      File vcf = "${sample_name}_somatic_unfiltered.vcf.gz"
      File bam = "${sample_name}_tumor_normal_m2.bam"
    }
}

task GetPileupSummaries {
  File tumor_bam
  File contamination_germline_resource
  String sample_name

  command {
    set -e
    ml biology gatk
    gatk GetPileupSummaries \
      -I ${tumor_bam} \
      -V ${contamination_germline_resource} \
      -L ${contamination_germline_resource} \
      -O ${sample_name}_pileupsummaries.table
  }
  runtime {
    runtime_minutes: "120"
    memory: 24000
  }
  output {
    File pileup_summary = "${sample_name}_pileupsummaries.table"
  }
}

task CalculateContamination {
  File pileup_summary
  String sample_name

  command {
    set -e
    ml biology gatk
    gatk CalculateContamination \
      -I ${pileup_summary} \
      -O ${sample_name}_calculatecontamination.table
  }
  runtime {
    runtime_minutes: "20"
    memory: 8000
  }
  output {
    File contamination_table = "${sample_name}_calculatecontamination.table"
  }
}

task FilterMutectCalls {
  File unfiltered_vcf 
  File contamination_table 
  String sample_name 

  command {
    set -e
    ml biology gatk
    gatk FilterMutectCalls \
      -V ${unfiltered_vcf} \
      --contamination_table ${contamination_table} \
      -O ${sample_name}_somatic_filtered.vcf.gz
  }
  runtime {
    runtime_minutes: "20"
    memory: 8000
  }
  output {
    File filtered_vcf = "${sample_name}_somatic_filtered.vcf.gz"
  }
}

task CopyOutputs {
  String output_directory 
  File filtered_vcf 
  File output_bam 
  File contamination_table 
  File pileup_summary 

  command {
    mkdir ${output_directory}
    cp ${filtered_vcf} ${output_directory}
    cp ${output_bam} ${output_directory}
    cp ${contamination_table} ${output_directory}
    cp ${pileup_summary} ${output_directory}
  }
  runtime {
    runtime_minutes: "240"
    memory: 16000
  }
}
