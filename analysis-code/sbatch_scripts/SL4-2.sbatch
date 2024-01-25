#!/usr/bin/env bash

#SBATCH -v -n 1 --mem=16000 -o SL4-2.out -e SL4-2.err

#min_cells_per_trip=40
min_trip_per_cell=5
max_umi_per_cell=1000000
min_cells_per_trip=5
min_readcount=10
min_umi_per_cell=50
#for min_cells_per_trip in 30 40 50 100
#do
id=$min_trip_per_cell"_"$min_cells_per_trip"_"$min_readcount"_"$min_umi_per_cell
echo ${id}

DIR=/scratch/bclab/aramu/CAS/full-run-101921/Combined/Cohen/
python3 /scratch/bclab/aramu/CAS/castools/castools/parse_reads.py ${DIR}/Cohen_SZ-SL4-2_SIC_70_ATGTCCCTAA_S2_R1_001.fastq.gz ${DIR}/Cohen_SZ-SL4-2_SIC_70_ATGTCCCTAA_S2_R2_001.fastq.gz | sort -u > SL4-2_cells_umis_tripbc.tsv
#python3 /scratch/bclab/aramu/CAS/castools/castools/error_correction.py SL4-2_cells_umis_tripbc.tsv SL4-2
#python3 /scratch/bclab/aramu/CAS/castools/castools/tripBC_ec.py SL4-2_ecoutput.tsv  CAS_LP1_pool2_DNA_barcodes.txt SL4-2_tripbcecoutput
#python3 /scratch/bclab/aramu/CAS/castools/castools/translate_cellbcs.py SL4-2_tripbcecoutput.tsv > SL4-2_tripbcecoutput_translated_inwhitelist.tsv
#python3 /scratch/bclab/aramu/CAS/castools/castools/oinb_estimation.py SL4-2_tripbcecoutput_translated_inwhitelist.tsv SL4-2_nodenoise_"$id"_oinboutput  $min_cells_per_trip $min_umi_per_cell $max_umi_per_cell $min_trip_per_cell
#python3 ../castools/castools/plot_bulk_sc.py CAS_LP1_pool2_expression.tsv SL4-2_nodenoise_"$id"_oinboutput SL4-2_nodenoise_"$id"_plots


#python3 ../castools/castools/compute_cellbc_overlap.py ./SL5-1-merged-v2/outs/filtered_feature_bc_matrix/barcodes.tsv.gz SL4-2_cells_umis_tripbc.tsv > SL4-2_cells_umis_tripbc_intranscriptome.tsv
#python3 /scratch/bclab/aramu/CAS/castools/castools/denoise_tBCs.py --quad=SL4-2_cells_umis_tripbc_intranscriptome.tsv --bulkbcs=CAS_LP1_pool2_expression.tsv --exp=SL4-2_cells_umis_tripbc_intranscriptome --rep=_denoise_"$id" --minreadcount=$min_readcount --hamming_lookup=hamming_lookup.pkl
#python3 /scratch/bclab/aramu/CAS/castools/castools/oinb_estimation.py SL4-2_cells_umis_tripbc_intranscriptome_denoise.tsv SL4-2_intranscriptome_denoise_"$id"_oinboutput  $min_cells_per_trip $min_umi_per_cell $max_umi_per_cell $min_trip_per_cell
#python3 ../castools/castools/plot_bulk_sc.py CAS_LP1_pool2_expression.tsv SL4-2_intranscriptome_denoise_"$id"_oinboutput.tsv SL4-2_intranscriptome_denoise_"$id"_plots
#done
