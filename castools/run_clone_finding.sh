#!/bin/bash
#
#SBATCH --job-name=clonal_trip
#SBATCH --mail-type=ALL
#SBATCH --mail-user=szhao@wustl.edu
#
#SBATCH --ntasks=4
#SBATCH --mem-per-cpu=8000
#SBATCH --output=report_network.out

ml miniconda3
eval "$(conda shell.bash hook)"
conda activate mascot 

python3 ~/castools/castools/find_trip_clones.py --trio --min_weight --name