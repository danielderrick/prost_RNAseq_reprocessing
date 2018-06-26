#!/usr/bin/env python

import os
import argparse
import json

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--raw", default="/home/exacloud/lustre1/HeiserLab/derrickd/RNASeq/korkola/raw")
    parser.add_argument("sample")
    parser.add_argument("index")
    args = parser.parse_args()

    fastqs = os.listdir(args.raw)
    for fastq in fastqs:
        if (args.sample in fastq and "R1.fastq.gz" in fastq):
            read1 = fastq
        elif (args.sample in fastq and "R2" in fastq):
            read2 = fastq

    job = {
        "TUMOR_FASTQ_1" : {
            "class": "File",
            "path": "%s/%s" % (args.raw, read1)
        },
        "TUMOR_FASTQ_2" : {
            "class": "File",
            "path": "%s/%s" % (args.raw, read2)
        },
        "index": {
            "class": "File",
            "path": "%s" % (args.index)
        }
    }

    print json.dumps(job)
