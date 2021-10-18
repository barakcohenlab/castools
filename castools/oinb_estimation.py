import pandas as pd
import logging
import matplotlib.pyplot as plt
import numpy as np
import statsmodels.api as sm
import argparse
import sys
from datetime import datetime


def extract_scTrip_fast(path, filename, min_umi=25, max_umi = 800000,
                        min_cells = 500, min_trip_percell = 5):
    '''
    Function to extract
    '''
    logger.info(f"Minimum UMI per cell {min_umi}")
    logger.info(f"Maximum UMI per cell {max_umi}")
    logger.info(f"Minimum trip per cell {min_trip_percell}")
    tripData = [] #Key is the cell id, value is a set of UMIs
    with open(path) as cellumis_tripbc_fh:
        for line in cellumis_tripbc_fh:
            trio = line.rstrip("\n").split('\t')
            tripData.append(trio)
    trio_to_process = {} # key is cell barcode, value is [umi, trip_bc]
    # Before doing anything we need to filter out the cells with min_umi.
    cell_to_trip = {}
    cell_to_umi = {}
    for line in tripData:
        if len(line) == 4:
            cell,umi,trip_bc, _ = line
        else:
            cell,umi,trip_bc = line
        if cell not in trio_to_process:
            trio_to_process[cell] = []
            cell_to_trip[cell] = set()
            cell_to_umi[cell] = set()
        trio_to_process[cell].append([umi,trip_bc])
        cell_to_trip[cell].add(trip_bc)
        cell_to_umi[cell].add(umi)
    cell_to_normalizationfactor = {}
    numis_in_cell = []
    for cell in cell_to_trip:
        logger.debug(f"cell, normalization factor  {cell}, {float(len(cell_to_umi[cell]))/len(cell_to_trip[cell])}")
        cell_to_normalizationfactor[cell] = float(len(cell_to_umi[cell]))/len(cell_to_trip[cell])
        numis_in_cell.append(len(cell_to_umi[cell]))
    logger.info(f'length of trio_to_process before is {len(trio_to_process)}')
    keys_to_remove = []
    logger.info(f"max umi_per_cell {max([len(x) for x in trio_to_process.values()])}")
    logger.info(f"mean umi_per_cell {np.mean([len(x) for x in trio_to_process.values()])}")
    logger.info(f"median umi_per_cell {np.median([len(x) for x in trio_to_process.values()])}")
    logger.info(f"std umi_per_cell {np.std([len(x) for x in trio_to_process.values()])}")
    trips_in_cell = []
    for key in trio_to_process.keys():
        logging.debug(key, len(trio_to_process[key]))
        trips_in_cell.append(len(cell_to_trip[key]))
        if len(trio_to_process[key]) < min_umi or len(trio_to_process[key]) > max_umi or len(cell_to_trip[key]) < min_trip_percell:
            keys_to_remove.append(key)
    plt.hist(np.log10(trips_in_cell))
    plt.xlabel("log10 trips per cell")
    plt.savefig(filename + "_trips_per_cell.png")
    plt.figure(1)
    plt.hist(np.log10(numis_in_cell))
    plt.savefig(filename + "_umis_per_cell.png")
    numis_in_cell.sort(reverse = True)
    print(np.arange(len(numis_in_cell)), numis_in_cell)
    plt.figure(2)
    plt.scatter(np.arange(len(numis_in_cell)), numis_in_cell)
    plt.savefig(filename + "_rank_umis_per_cell.png")
    plt.figure(3)
    plt.scatter(np.arange(len(numis_in_cell)), np.log10(numis_in_cell))
    plt.savefig(filename + "_rank_log_umis_per_cell.png")
    for k in keys_to_remove:
        del trio_to_process[k]
    logger.info(f'length of trio_to_process after is {len(trio_to_process)}')
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
    logger.info(f"Minimum number of cells per trip {min_cells}")
    for trip_bc in trip_counts:
        counts = []
        for cell in trip_counts[trip_bc]:
            count = len(trip_counts[trip_bc][cell])
            if count >= 0:
                # Here I am artificially creating a zero-inflated negative binomial 
                # By left shift the distribution by 1
                counts.append(count/cell_to_normalizationfactor[cell])
        # Only record tripBC with 500 measurement
        if len(counts) >= min_cells:
            trip_cells_umi[trip_bc] = counts
    return trip_cells_umi

def get_oinb_estimate(scTRIP):
    keys = scTRIP.keys()
    key_list = []
    mu_list = []
    alpha_list = []
    mean_list = []
    median_list = []
    auc_list = []
    var_list = []
    ncells_list = []
    counts_list = []
    sum_list = []
    for key in keys:
        logger.info(f'we are dealing with {key}')
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
            zero_counts = [x+ 1 for x in zero_counts] 
            #original_counts = zero_counts + non_zero_counts
            original_counts = counts #Keep normalized coutns
            ncells_list.append(len(counts))
            counts_list.append(list(counts))
            mean_list.append(np.mean(original_counts))
            median_list.append(np.median(original_counts))
            sum_list.append(np.sum(original_counts))
            var_list.append(np.var(original_counts))
            auc_list.append(get_auc(original_counts))
    pop_df = pd.DataFrame([key_list, mean_list, median_list, var_list, auc_list, mu_list, alpha_list, ncells_list, counts_list, sum_list])
    pop_df = pop_df.transpose()
    pop_df.columns = ['tBC', 'mean', 'median', 'var','auc', 'mu', 'alpha', 'ncells', 'counts', 'sum']
    return pop_df

def get_auc(list):
    '''
    Helper function that helps with getting empirical AUC for a given histogram
    '''
    values, bins, _ = plt.hist(list)
    area = sum(np.diff(bins)*values)
    return area

# Wrapper function for OINF estimation from the distributions
def oinb_estimation(path, filename, args):
    '''
    Wrapper function for returning a tsv file with columns:
    'tripBC', 'mean', 'var', 'auc', 'mu', 'alpha'
    '''
    trip_counts = extract_scTrip_fast(path, filename, args.min_umi_percell, args.max_umi_percell, args.min_cells, args.min_trip_percell)
    scTRIP_stats = get_oinb_estimate(trip_counts)
    scTRIP_stats.to_csv(filename + '.tsv', index = None, sep = '\t')

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
    parser.add_argument(
        'min_cells',
        help="Filter TRIP barcodes that are not seen in atleast min_cells",
        type = int,
        default = 500
    )
    parser.add_argument(
        'min_umi_percell',
        help="Filter cells that have less than max_umi_percell UMIs",
        type = int,
        default = 25
    )
    parser.add_argument(
        'max_umi_percell',
        help="Filter cells that have more than max_umi_percell UMIs",
        type = int,
        default = 800000
    )
    parser.add_argument(
        'min_trip_percell',
        help="Filter cells that have less than min_trip_percell TRIP barcodes",
        type = int,
        default = 5
    )
    args = parser.parse_args()
    oinb_estimation(args.trio, args.file_name, args)


if __name__ == "__main__":
    logger = logging.getLogger('oinb_log')
    logger.setLevel(logging.DEBUG)

    # create console handler and set level to debug
    ch = logging.StreamHandler()
    ch.setLevel(logging.DEBUG)

    # create formatter
    formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')

    # add formatter to ch
    ch.setFormatter(formatter)

    # add ch to logger
    logger.addHandler(ch)
    main()
