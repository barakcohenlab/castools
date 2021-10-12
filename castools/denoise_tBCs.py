import os,sys
import pandas as pd
import numpy as np
import pickle
import argparse


def calculate_hamming(rBC1, good_list):
    for rBC in good_list:
        hamming = hammingDist(rBC1, rBC)
        if hamming < 2:
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

def filter_based_on_umi(quint_df, filter_BC_list, min_count = 2):
    '''
    Input：quint_df, a five-column table of cellBC, prom_id, pBC, rBC, umi, counts, and clusters the cell come from.
    Output: the unique quint that are supported by more than 1 read. 
    '''
    cell_list = list(set(quint_df['cellBC'].values))
    pop_list = []
    for cells in cell_list:
        cell_df = quint_df[quint_df['cellBC'] == cells]
        cell_list = list(set(quint_df['cellBC'].values))
        # Get the high count slice
        high_count_slice = cell_df[cell_df['count'] > min_count]
        # Filter the high count slice based on the filter_BC_list
        high_count_slice = high_count_slice[high_count_slice['tBC'].isin(filter_BC_list)]
        # Return a list for the high count slice
        high_count_tBC = set(high_count_slice['tBC'].values)
        for _, row in cell_df.iterrows():
            if row['tBC'] in high_count_tBC:
                pop_list.append([row['cellBC'], row['tBC'] , row['umi'], row['count'] ])
            else:
                possible_rBC = calculate_hamming(row['tBC'], high_count_tBC)
                if possible_rBC != 0:
                    pop_list.append([row['cellBC'], row['umi'], row['tBC'] , row['count']])
    return pd.DataFrame(pop_list, columns = ['cellBC', 'umi' , 'tBC', 'count'])

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--quad', help='tsv file contains the quad information', required=True)
    parser.add_argument('--bulkbcs', help='path to a tsv file that contains the tripBC that can be measured in bulk', required=True)
    parser.add_argument('--rep', help='output rep name', required=True)
    parser.add_argument('--exp', help = 'name of the experiment', required = True)
    # Grab input arguments
    args= parser.parse_args()
    # Read in the Quad that contains the [cellBC, umi, tripBC, count], note even this is a tsv file
    # The delimiter is a space
    prom_quint = pd.read_csv(args.quad, sep= ' ', names = ['cellBC', 'umi', 'tBC', 'count'])
    # Read in the bulk barcodes
    bulk_bcs_exp = pd.read_csv(args.bulkbcs, sep = '\t')
    # Process the bulk BCs to make it into a list
    bulk_bcs = bulk_bcs_exp['tBC'].to_list()
    # Denoise the cell based data
    pop_df = filter_based_on_umi(prom_quint, bulk_bcs,2)
    # Save it to a tsv file
    pop_df.to_csv(args.exp + args.rep + '.tsv', sep = '\t', index = False)
if __name__ == '__main__':
    main()
