#!/usr/bin/env cwl-runner
#
# Authors: Allison Creason, John Letaw, Janice Patterson

cwlVersion: v1.0
class: CommandLineTool
baseCommand: ["samtools", "view"]

doc: "samtools view for sam to bam"

hints:
    DockerRequirement:
        dockerPull: alliecreason/samtools:1.6.0

requirements:
    - class: InlineJavascriptRequirement

inputs:
  outbam:
    type: boolean
    inputBinding:
      position: 2
      prefix: "-b"

  header:
    type: boolean
    inputBinding:
      position: 2
      prefix: "-h"

  threads:
    type: int
    inputBinding:
      position: 2
      prefix: "-@"

  input:
    type: File
    inputBinding:
      position: 4

  output_name:
    type: string
    inputBinding:
      position: 3
      prefix: "-o"

outputs:
  output:
    type: File
    outputBinding:
      glob: "*.bam"
