#!/usr/bin/env cwl-runner
#
# Author: Allison Creason

cwlVersion: v1.0
class: CommandLineTool
baseCommand: ["kallisto","quant"]


doc: "kallisto quantification algorithm"
hints:
  DockerRequirement:
    dockerPull: alliecreason/kallisto:0.43.0

requirements:
  - class: InlineJavascriptRequirement

inputs:
  index:
    type: File
    inputBinding:
      position: 1
      prefix: "--index"

  fastq1:
    type: File
    inputBinding:
      position: 2 

  fastq2:
    type: File
    inputBinding:
      position: 3

  outputDir:
    type: string
    inputBinding:
      position: 1
      prefix: "--output-dir"

  bias:
    type: boolean
    inputBinding:
      position: 1
      prefix: "--bias"

  bootstrap:
    type: int
    inputBinding:
      position: 1
      prefix: "--bootstrap-samples="
      separate: false

  plaintext:
    type: boolean
    inputBinding:
      position: 1
      prefix: "--plaintext"

  fusion:
    type: boolean
    inputBinding:
      position: 1
      prefix: "--fusion"

  single:
    type: boolean
    inputBinding:
      position: 1
      prefix: "--single"

  frstranded:
    type: boolean
    inputBinding:
      position: 1
      prefix: "--fr-stranded"

  rfstranded:
    type: boolean
    inputBinding:
      position: 1
      prefix: "--rf-stranded"

  threads:
    type: int
    inputBinding:
      position: 1
      prefix: "--threads="
      separate: false

  pseudobam:
    type: boolean
    inputBinding:
      position: 1
      prefix: "--pseudobam"

outputs:
  abundance_tsv:
    type: File
    outputBinding:
      glob: $(inputs.outputDir + "/*.tsv")
  abundance_h5:
    type: File
    outputBinding:
      glob: $(inputs.outputDir + "/*.h5")

  pseudosam:
    type: stdout


stdout: pseudo.sam
