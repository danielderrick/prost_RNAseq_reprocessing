#!/usr/bin/env cwl-runner
#
# Authors: Allison Creason, John Letaw, Janice Patterson

class: Workflow
cwlVersion: v1.0

doc: "Kallisto Isoform quantification workflow for Paired-end"

requirements:
  - class: MultipleInputFeatureRequirement

inputs: 

  index: File

  TUMOR_FASTQ_1: File

  TUMOR_FASTQ_2: File

outputs:

  OUTPUT:
    type: 
      type: array
      items: File
    linkMerge: merge_flattened
    outputSource: [trim-galore/reports, trim-galore/fastq1-trimmed, trim-galore/fastq2-trimmed, kallisto/abundance_tsv, kallisto/abundance_h5,  samtools/output]

steps:
  trim-galore:
    run: ../tools/trim_galore/trim-galore.cwl
    in:
      paired: { default: true }
      length: { default: 35 }
      fastq1: TUMOR_FASTQ_1
      fastq2: TUMOR_FASTQ_2

    out: [fastq1-trimmed, fastq2-trimmed, reports]

  kallisto:
    run: ../tools/kallisto/kallisto-quant.cwl
    in:
      index: index
      fastq1: trim-galore/fastq1-trimmed
      fastq2: trim-galore/fastq2-trimmed
      outputDir: { default: "out" }
      bias: { default: true }
      bootstrap: { default: 30 }
      plaintext: { default: false }
      fusion: { default: false }
      single: { default: false }
      frstranded: { default: false }
      rfstranded: { default: false }
      threads: { default: 1 }
      pseudobam: { default: true }

    out: [abundance_tsv, abundance_h5, pseudosam]

  samtools:
    run: ../tools/samtools/samtools.cwl
    in:
      outbam: { default: true }
      header: { default: true }
      threads: { default: 8 } 
      output_name: { default: "pseudo.bam" }
      input: kallisto/pseudosam

    out: [output]

