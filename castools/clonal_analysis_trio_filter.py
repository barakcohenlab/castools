import os,sys
import logging
from numpy.core.defchararray import count
import pandas as pd
import numpy as np
import gzip
from scipy import stats
import pickle
import argparse
from datetime import datetime

def calculate_hamming(rBC1, good_list):
    for rBC in good_list:
        hamming = hammingDist(rBC1, rBC)
        if hamming < 3:
            return rBC
        else:
            return 0


def hammingDist(str1, str2):
    '''
    Calculating hamming distance of two strings
    https://www.geeksforgeeks.org/hamming-distance-two-strings/
    '''
    i = 0
    count = 0
    while(i < len(str1)):
        if(str1[i] != str2[i]):
            count += 1
        i += 1
    return count


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
    temp_raw_trio = pd.read_csv(path, sep = '\t')
    # Filter out the cells that are not in the final_trio_df
    filtered_key_list = []
    for _, value in trip_counts.items():
        for key, _ in value.items():
            filtered_key_list.append(key)
    filtered_key_list = list(set(filtered_key_list))
    temp_raw_trio = temp_raw_trio[temp_raw_trio['cellBC'].isin(filtered_key_list)]
    # Save the filtered trio file
    return temp_raw_trio

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
    trip_counts = extract_scTrip_fast(args.trio, args.file_name, args.min_umi_percell, args.max_umi_percell, args.min_trip_percell)
    trip_counts.to_csv(f'{args.file_name}_filtered_trio.csv', sep = '\t', index = False)

if __name__ == "__main__":
    logger = logging.getLogger('filter_trio_log')
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
