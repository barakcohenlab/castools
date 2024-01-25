#!/usr/bin/env bash

#SBATCH -v -n 1 --mem=16000 -o techcor.out -e techcor.err


for f1 in `/bin/ls -1 test-different-cutoffs/old/*_denoise*replicates*tsv`
do
   echo $f1
   Rscript compute_technical_correlation.R $f1
done
