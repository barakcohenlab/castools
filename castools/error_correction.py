"""Error correcting based on hamming distance 

This script processes the TRIP:[CellBC, UMI, TRIPBC]_i and error correct and collapse
based on hamming distance. 

This file takes a list of TRIP:[CellBC, UMI, TRIPBC] that is generated from x.py.
This file will output the collapsed barcodes for further processing.

"""
import numpy as np 
import os,sys
import argparse
import pandas as pd
import csv      
from collections import Counter


def counting_trios(trio):
    '''
    Function returns the counted data from the cells
    '''
    # Return the cell barcode : read_count
    cell_bc_counts_dict = {item[0]:item[-1]for item in trio}
    # Sort the dict into an order dict
    sorted_cell_bc = {k: v for k, v in sorted(cell_bc_counts_dict.items(), key=lambda item: item[1])}
    return sorted_cell_bc

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

def hamming_key_list(key, key_list):
    '''
    Helper function to compare the 
    '''
    for main_key in key_list:
        if hammingDist(key, main_key) <=6:
            return main_key
        else: 
            return key


def collasping_cell_BCs(sorted_cell_bc):
    '''
    Input: sorted dict of cell barcodes
    Output: dictionary of cell barcodes with key(corrected BC):value(matching BCs)
    '''
    # Initiate a list of dictionary to hold the correct call barcodes
    correspondence = {}
    # Get the keys
    search_keys = list(sorted_cell_bc.keys())
    for i in range(len(search_keys)):
        key = search_keys[i]
        if len(correspondence.keys()) == 0:
            correspondence[key] = [key]
        else:
            new_key = hamming_key_list(key, correspondence.keys())
            if new_key in correspondence.keys():
                val = correspondence[new_key]
                new_val = val.append(key)
                correspondence.update(new_key = new_val)
            else:
                correspondence[key] = [key]
        i += 1
    return correspondence

def return_correct_cell_bc(cell_bc, collapsed_cell_bc):
    for key in collapsed_cell_bc.keys():
        if cell_bc in collapsed_cell_bc[key]:
            return key
        else:
            return cell_bc

def correct_cell_bc(trio, collasped_cell_bc):
    '''
    Match the cell Barcode with the collapsed_cell_bc
    '''
    pop_list = []
    for line in trio:
        try:
            cell_bc, umi, trip_bc, count = line
        except:
            print(line)
            continue
        new_cell_bc = return_correct_cell_bc(cell_bc, collasped_cell_bc)
        pop_list.append([new_cell_bc,umi,trip_bc, count])
    return pop_list

def error_correcting(trio_tsv, prefix):
    trio = []
    with open(trio_tsv) as f:
        for line in f:
            trio.append(line.strip().split())
    # Return CellBCs
    sorted_cell_bc = counting_trios(trio)
    collapsed_cell_bc = collasping_cell_BCs(sorted_cell_bc)
    pop_trio = correct_cell_bc(trio, collapsed_cell_bc)
    pop_trio.sort(key = lambda x: x[3]) 
    with open(prefix + "_" + 'output.tsv', 'w', newline= '') as f_output:
        tsv_output = csv.writer(f_output, delimiter='\t')
        tsv_output.writerows(pop_trio)



# Main function 
def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        'trio',
        help="The list of TRIOs"
    )
    parser.add_argument(
        'prefix',
        help="Prefix for the output file",
        default = ""
    )
    args = parser.parse_args()
    error_correcting(args.trio, args.prefix)


if __name__ == "__main__":
    main()
