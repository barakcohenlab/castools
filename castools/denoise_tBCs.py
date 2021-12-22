import os,sys
import pandas as pd
import numpy as np
import pickle
import argparse



write_hamming_lookup = False
def calculate_hamming(hamming_lookup, rBC1, good_list):
    for rBC in good_list:
        pair1 = ":".join(sorted([rBC, rBC1]))
        if pair1 in hamming_lookup:
                hamming = hamming_lookup[pair1]
        else:
            hamming = hammingDist(rBC1, rBC)
            hamming_lookup[pair1] = hamming
            write_hamming_lookup = True #some new pairs added
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

def filter_based_on_umi(quint_df, filter_BC_list, hamming_lookup, min_count = 2):
    '''
    Input：quint_df, a five-column table of cellBC, prom_id, pBC, rBC, umi, counts, and clusters the cell come from.
    Output: the unique quint that are supported by more than 1 read. 
    '''
    print("minimum read count is > ", min_count, file = sys.stderr)
    print("Reading hamming lookup: ", hamming_lookup, file = sys.stderr)
    if hamming_lookup is not None:
        with open(hamming_lookup, 'rb') as hl_fh:
            hamming_lookup = pickle.load(hl_fh)
    else:
        hamming_lookup = {}
        write_hamming_lookup = True
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
                pop_list.append([row['cellBC'], row['umi'] , row['tBC'], row['count'] ])
            else:
                possible_rBC = calculate_hamming(hamming_lookup, row['tBC'], high_count_tBC)
                if possible_rBC != 0:
                    pop_list.append([row['cellBC'], row['umi'], possible_rBC , row['count']])
                    #add this denoised tripBC to the high count set since it looks good
                    high_count_tBC.add(row['cellBC'])
    if write_hamming_lookup:
        print("Saving hamming lookup to hamming_lookup.pkl", file = sys.stderr)
        with open('hamming_lookup.pkl', 'wb') as fh:
            pickle.dump(hamming_lookup, fh)
    return pd.DataFrame(pop_list, columns = ['cellBC', 'umi' , 'tBC', 'count'])

def main():
    print("Minimum hamming cutoff for tripBC is 2", file = sys.stderr)
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--quad', help='tsv file contains the quad information', required=True)
    parser.add_argument('--bulkbcs', help='path to a tsv file that contains the tripBC that can be measured in bulk', required=True)
    parser.add_argument('--rep', help='output rep name', required=True)
    parser.add_argument('--exp', help = 'name of the experiment', required = True)
    parser.add_argument('--minreadcount', help = 'min read-count for quad', required = True, default = 2, type = int)
    parser.add_argument('--hamming_lookup', help = 'pickled dictionary of hamming distances', required = False)
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
    pop_df = filter_based_on_umi(prom_quint, bulk_bcs, args.hamming_lookup, args.minreadcount)
    # Save it to a tsv file
    pop_df.to_csv(args.exp + args.rep + '.tsv', sep = '\t', index = False)

if __name__ == '__main__':
    main()
