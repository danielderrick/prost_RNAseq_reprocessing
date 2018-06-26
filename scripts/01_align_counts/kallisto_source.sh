#!/bin/bash
#
# Authors: Allison Creason, John Letaw, Janice Patterson


# This script is meant to function within the BioCoders environment.
BIOCODERS="/home/exacloud/lustre1/BioCoders"
source $BIOCODERS/USE_APPS --python3
export PATH="$PATH:$BIOCODERS/Applications/anaconda2/bin"
source activate cwltool

# Create these before running the script.
PROJ="/home/exacloud/lustre1/HeiserLab/derrickd/RNASeq/korkola"
RAW="${PROJ}/raw"
LOGS="${PROJ}/logs"
OUTPUT="${PROJ}/output"
SCRIPTS="${PROJ}/scripts"

# Inputs
INDEX="/home/exacloud/lustre1/BioCoders/DataResources/Transcriptomes/Human/gencode_release24/kallisto_index/gencodev24_kallistoindex"
WORKFLOW="${SCRIPTS}/workflows/kallisto-quant-paired-workflow.cwl"

