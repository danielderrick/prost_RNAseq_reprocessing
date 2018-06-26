#!/bin/bash
#
# Authors: Allison Creason, John Letaw, Janice Patterson

SAMPLE=$1
echo $SAMPLE
echo $INDEX
echo $OUTPUT/$SAMPLE
INPUT="input.json"

# Set up env vars
source /home/exacloud/lustre1/HeiserLab/derrickd/RNASeq/korkola/scripts/kallisto_source.sh

# Generate input json file
if [ ! -e $OUTPUT/$SAMPLE ]; then
    mkdir -p $OUTPUT/$SAMPLE
fi

$SCRIPTS/generate_job.py $SAMPLE $INDEX --raw $RAW > $OUTPUT/$SAMPLE/$INPUT

# Run cwl workflow
cwltool --tmp-outdir-prefix /mnt/scratch/ --outdir $OUTPUT/$SAMPLE --rm-container --rm-tmpdir $WORKFLOW $OUTPUT/$SAMPLE/$INPUT

find /mnt/scratch/ -user derrickd -exec rm -r {} \;

