#!/usr/bin/env cwl-runner
#
# Authors: Allison Creason, John Letaw, Janice Patterson

cwlVersion: v1.0
class: CommandLineTool
baseCommand: ["trim_galore"]

doc: "Trim Galore"

hints:
    DockerRequirement:
        dockerPull: alliecreason/trim_galore:1.4.3

requirements:
    - class: InlineJavascriptRequirement

inputs:
  paired:
    type: boolean
    inputBinding:
      position: 1
      prefix: "--paired"

  length:
    type: int
    inputBinding:
      position: 1
      prefix: "--length"

  fastq1:
    type: File
    inputBinding:
      position: 2

  fastq2:
    type: File
    inputBinding:
      position: 3

outputs:
   fastq1-trimmed:
     type: File
     outputBinding:
       glob: "*val_1.fq.gz"

   fastq2-trimmed:
     type: File
     outputBinding:
       glob: "*val_2.fq.gz"

   reports:
     type:
       type: array
       items: File
     outputBinding:
       glob: "*trimming_report.txt"
