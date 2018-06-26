#!/bin/bash

#SBATCH -N 1
#SBATCH -n 1
#SBATCH --time=24:00:00
#SBATCH --mem 15GB
#SBATCH --cpus-per-task 8
#SBATCH --job-name=run_kallisto
#SBATCH --mail-type=ALL
#SBATCH --gres disk:335
#SBATCH -o /home/exacloud/lustre1/HeiserLab/derrickd/RNASeq/korkola/logs/%j_dockerkallisto.%N.out
#SBATCH -e /home/exacloud/lustre1/HeiserLab/derrickd/RNASeq/korkola/logs/%j_dockerkallisto.%N.err

bash /home/exacloud/lustre1/HeiserLab/derrickd/RNASeq/korkola/scripts/run_kallisto.sh $1
