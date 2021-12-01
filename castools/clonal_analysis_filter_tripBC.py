# basic packages
import os,sys
import logging
from importlib import reload 
import argparse
from datetime import datetime
# functional packages
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from kneed import KneeLocator
# Custom packages

def return_filtered_tripBC(trio):
    '''
    This function filters the tripBC dataframe by the trio provided
    '''
    cell_list = list(set(trio['cellBC'].values))
    pop_list_tot = pd.DataFrame(columns = ['cellBC', 'tBC' , 'norm_depth'])
    # Almost solely run this piece of code on one cell only. 
    for cells in cell_list:
        cell_df = trio[trio['cellBC'] == cells]
        # Calculate the normalized depth of each cell.
        # Calculate the normalized depth for each tBC
        count_list = cell_df['norm_depth'].to_list()
        if len(count_list) >0:
            minimum_depth = find_knee(count_list)
            # remove the cells with low depth
            if minimum_depth != 0:
                cell_df  = cell_df[cell_df['norm_depth'] > minimum_depth]
                pop_list_tot = pd.concat([pop_list_tot, cell_df])
    return pop_list_tot    

def find_knee(count_list):
    '''
    helper function to find the knee of the normalized depth
    '''
    x = np.arange(len(count_list))
    kl = KneeLocator(x, count_list, curve = 'convex', direction='decreasing', interp_method='interp1d')
    knee_point = count_list[kl.knee]
    return knee_point

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
    trip_counts = return_filtered_tripBC(args.trio)
    trip_counts.to_csv(f'{args.file_name}_filtered_tBC.csv', sep = '\t', index = False)

if __name__ == "__main__":
    logger = logging.getLogger('filter_tBC_log')
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
