# prostate_line_panel_RNAseq

This directory contains the scripts and intermediate files used to reprocess 
RNAseq data from a panel of 20 prostate cell lines

## Contents

* counts/: gene expression matrices; fpkm and log2(fpkm)

* reports/: summary report pdf and fastQC reports

* plots/: pdf copies of plots in reports/

* scripts/01_align_counts/: scripts used to align the FASTQ files to the Gencode
V24 transcriptome using the Kallisto pseudo-alignment software. These scripts 
are to be run on Exacloud. Submit.txt is the master file to generate and submit 
slurm batch jobs.

* scripts/02_summarize_counts/: R scripts used to import the Kallisto abundance 
estimates, normalize samples for library size differences, summarize transcripts
to the gene-level, and add gene symbol annotation.

* data/: the Kallisto .h5 abundance estimates
