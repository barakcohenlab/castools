import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import statsmodels.api as sm
import argparse
from datetime import datetime
import csv


def extract_scTrip_fast(path,min_umi=25, max_umi = 800000):
    '''
    Function to extract 
    '''
    tripData = [] #Key is the cell id, value is a set of UMIs
    with open(path) as cellumis_tripbc_fh:
        for line in cellumis_tripbc_fh:
            trio = line.rstrip("\n").split('\t')
            tripData.append(trio)
    trio_to_process = {}
    # Before doing anything we need to filter out the cells with min_umi.
    for line in tripData:
        if len(line) == 4:
            cell,umi,trip_bc, _ = line
        else: 
            cell,umi,trip_bc = line
        if cell not in trio_to_process:
            trio_to_process[cell] = []
        trio_to_process[cell].append([umi,trip_bc])
    print(f'length of trio_to_process before is {len(trio_to_process)}')
    keys_to_remove = []
    for key in trio_to_process.keys():
        if len(trio_to_process[key]) < min_umi or len(trio_to_process[key]) > max_umi:
            keys_to_remove.append(key)
    for k in keys_to_remove:
        del trio_to_process[k]
    print(f'length of trio_to_process after is {len(trio_to_process)}')
    # Now create a new dictionary:
    trip_counts = {}
    for key in trio_to_process.keys():
        cell = key
        for umi_trip_pair in trio_to_process[key]:
            umi, trip_bc = umi_trip_pair
            if trip_bc not in trip_counts:
                trip_counts[trip_bc] = {}
            if cell not in trip_counts[trip_bc]:
                trip_counts[trip_bc][cell] = []
            trip_counts[trip_bc][cell].append(umi)
    trip_cells_umi = {}
    for trip_bc in trip_counts:
        counts = []
        for cell in trip_counts[trip_bc]:
            count = len(trip_counts[trip_bc][cell])
            if count >=0:
                # Here I am artificially creating a zero-inflated negative binomial 
                # By left shift the distribution by 1
                counts.append(count-1)
        # Only record tripBC with 500 measurement
        if len(counts) > 500:
            trip_cells_umi[trip_bc] = counts
    return trip_cells_umi

def get_oinb_estimate(scTRIP):
    keys = scTRIP.keys()
    key_list = []
    mu_list = []
    alpha_list = []
    mean_list = []
    auc_list = []
    var_list = []
    for key in keys:
        print(f'we are dealing with {key}')
        counts = list(scTRIP[key])
        if max(counts) != 0:
            key_list.append(key)
            counts = pd.Series(counts)
            res = sm.ZeroInflatedNegativeBinomialP(counts, np.ones_like(counts)).fit( maxiter=200) 
            alpha_list.append(res.params['alpha'])
            mu_list.append(res.params['const'])
            # Get the original list with ones 
            zero_counts = [a for a in counts if a == 0]
            non_zero_counts = [a for a in counts if a > 0]
            zero_counts = zero_counts + 1
            original_counts = zero_counts + non_zero_counts
            mean_list.append(np.mean(original_counts))
            var_list.append(np.var(original_counts))
            auc_list.append(get_auc(original_counts))
    pop_df = pd.DataFrame([key_list, mean_list, var_list, auc_list, mu_list, alpha_list])
    pop_df = pop_df.transpose()
    pop_df.columns = ['tBC', 'mean', 'var','auc', 'mu', 'alpha']
    return pop_df

def get_auc(list):
    '''
    Helper function that helps with getting empirical AUC for a given histogram
    '''
    values, bins, _ = plt.hist(list)
    area = sum(np.diff(bins)*values)
    return area

# Wrapper function for OINF estimation from the distributions
def oinb_estimation(path, filename):
    '''
    Wrapper function for returning a tsv file with columns:
    'tripBC', 'mean', 'var', 'auc', 'mu', 'alpha'
    '''
    trip_counts = extract_scTrip_fast(path)
    scTRIP_stats = get_oinb_estimate(trip_counts)
    with open(filename + '.tsv', 'w', newline= '') as f_output:
        tsv_output = csv.writer(f_output, delimiter='\t')
        tsv_output.writerows(scTRIP_stats)
# Main function 
def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        'trio',
        help="Path to the TRIOs"
    )
    parser.add_argument(
        'file_name',
        help=" FileName that you want to save the file with, no suffix",
        default = f'{datetime.now().strftime("%d/%m/%Y %H:%M:%S")} scTRIP_stats'
    )
    args = parser.parse_args()
    oinb_estimation(args.trio, args.file_name)


if __name__ == "__main__":
    main()