#!/usr/bin/env bash

#SBATCH -v -n 1 --mem=16000 -o SL4-1.out -e SL4-1.err

min_cells_per_trip=10
min_umi_per_cell=5
max_umi_per_cell=1000000
min_trip_per_cell=5

python3 /scratch/bclab/aramu/CAS/castools/castools/oinb_estimation.py SL4-1_cells_umis_tripbc_intranscriptome_denoise.tsv SL4-1_cells_umis_tripbc_intranscriptome_denoise_oinboutput  $min_cells_per_trip $min_umi_per_cell $max_umi_per_cell $min_trip_per_cell
python3 ../castools/castools/plot_bulk_sc.py CAS_LP1_pool2_expression.tsv SL4-1_cells_umis_tripbc_intranscriptome_denoise_oinboutput.tsv SL4-1_cells_umis_tripbc_intranscriptome_denoise_bulk_sc
python3 /scratch/bclab/aramu/CAS/castools/castools/oinb_estimation.py SL4-2_cells_umis_tripbc_intranscriptome_denoise.tsv SL4-2_cells_umis_tripbc_intranscriptome_denoise_oinboutput  $min_cells_per_trip $min_umi_per_cell $max_umi_per_cell $min_trip_per_cell
python3 ../castools/castools/plot_bulk_sc.py CAS_LP1_pool2_expression.tsv SL4-2_cells_umis_tripbc_intranscriptome_denoise_oinboutput.tsv SL4-2_cells_umis_tripbc_intranscriptome_denoise_bulk_sc
