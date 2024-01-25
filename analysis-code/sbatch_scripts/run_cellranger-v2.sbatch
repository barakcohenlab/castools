#!/usr/bin/env bash

#SBATCH -v -c 8 --mem=64000 -o $2".out" -e $2".err"

fastqfolder=$1
fastq=$2
sample=$3
cells=$4

ml cellranger/6.0.1
cellranger count --id=$sample \
                 --fastqs=$fastqfolder \
                 --transcriptome=/scratch/bclab/aramu/dat/cellranger-index/refdata-gex-GRCh38-2020-A/ \
                 --sample=$fastq \
                 --localcores=8 \
                 --expect-cells=$cells
