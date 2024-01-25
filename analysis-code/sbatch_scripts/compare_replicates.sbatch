#!/usr/bin/env bash

#SBATCH -v -n 1 --mem=8000 -o compare-replicates.out -e compare-replicates.err


for min_trip_per_cell in 5
do
for min_cells_per_trip in 5
do
for min_readcount in 10
do
for min_umi_per_cell in 50
do
    id=$min_trip_per_cell"_"$min_cells_per_trip"_"$min_readcount"_"$min_umi_per_cell
    echo ${id}
    #python3 /scratch/bclab/aramu/CAS/castools/castools/compare_replicates.py SL4-1_oinboutput_"$id".tsv SL4-2_oinboutput_"$id".tsv SL4_nodenoise_comparereplicates_"$id"
    python3 /scratch/bclab/aramu/CAS/castools/castools/compare_replicates.py SL4-1_intranscriptome_denoise_"$id"_oinboutput.tsv SL4-2_intranscriptome_denoise_"$id"_oinboutput.tsv SL4_denoise_comparereplicates_"$id"
done
done
done
done
