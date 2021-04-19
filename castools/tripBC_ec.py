import os, sys
import numpy as np
# Plotting stuff
import argparse
import pandas as pd
import pickle

def read_bulk_bc(path):
    bulk_bc_list = []
    with open(path) as f:
        next(f)
        for line in f:
            bulk_bc = line.rstrip('\n')
            if len(bulk_bc) == 16:
                bulk_bc_list.append(bulk_bc)
            else:
                print('not 16:', bulk_bc)
    return bulk_bc_list

def read_in_trio_list(path):
    tripData = []
    with open(path) as cellumis_tripbc_fh:
        for line in cellumis_tripbc_fh:
            trio = line.rstrip("\n").split('\t')
            if len(trio) == 4:
                cell,umi,trip_bc, _ = trio
            else: 
                cell,umi,trip_bc = trio                
            tripData.append([cell,umi,trip_bc])
    return tripData

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


def bulk_error_correction(sc_tripBC_list, bulk_tripBC_list, maximum_error = 8):
    '''
    Correct the scTRIP tripBC with bulk tripBC, with maximum hamming distance. 
    Also return the non-corrected ones.
    '''
    correspondence = dict.fromkeys(sc_tripBC_list)
    i = 0
    j = 0
    for sc_tripBC in sc_tripBC_list:
        #check if it's closer to any of the bulk_tripBCs
        matched_list = {}
        for bulk_tripBC in bulk_tripBC_list:
            # initiate a list by value
            if hammingDist(sc_tripBC, bulk_tripBC) <= maximum_error:
                #print('successful!')
                matched_list[bulk_tripBC] = hammingDist(sc_tripBC,bulk_tripBC)
        # check if we have any match
        if len(matched_list) == 0:
            # print('we see a wrong BC')
            i += 1
        else:
            # sort the matched barcode 
            sorted_matched_list = {k:v for k,v in sorted(matched_list.items(), key = lambda item:item[1])}
            # check if the closest match is unique
            smallest_match = min(sorted_matched_list.values())
            num_smallest_match = list(sorted_matched_list.values()).count(smallest_match)
            if num_smallest_match == 1:
                # unique small match
                correspondence[sc_tripBC] = list(sorted_matched_list.keys())[list(sorted_matched_list.values()).index(smallest_match)]
            else:
                j += 1
        # Count non-match barcodes

    print(f'The total number of not assigned sc tripBC is {i}, and it is {i/len(sc_tripBC_list)}')
    print(f'The total number of ambiguous sc tripBC is {j}, and it is {j/len(sc_tripBC_list)}')
    return correspondence

def write_correct_tripBCs(trios, collapsed_tripBCs):
    '''
    Match new trios with the collapsed tripBCs
    '''
    pop_list = []
    for line in trios:
        try:
            cell_bc, umi, trip_bc= line
        except:
            print(line)
            continue
        if collapsed_tripBCs[trip_bc] != None:
            new_tripBC = collapsed_tripBCs[trip_bc]
        else:
            new_tripBC = trip_bc
        pop_list.append([cell_bc,umi,new_tripBC])        
    return pop_list

def error_correction_with_bulk(trios_path, bulk_path):
    '''
    Function to correct tripBC error with bulk sequenced tripBCs
    Input: path to the trioFiles (with read count, i.e. 4 columns), path to 
    bulk mapped tripBCs (txt file, with header)
    Output: a csv file with corrected tripBC, uncorrected trios remain the same.
    '''
    bulk_bc_list = read_bulk_bc(bulk_path)
    trios = read_in_trio_list(trios_path)
    sc_tripBC_list = [a[2] for a in trios]
    print('we are going to correct sc tripBC')
    corrected_sc_tripBC_dict = bulk_error_correction(sc_tripBC_list, bulk_bc_list)
    corrected_trios = write_correct_tripBCs(trios, corrected_sc_tripBC_dict)
    corrected_trios = pd.DataFrame(corrected_trios, columns= ['cellBC', 'umi', 'tripBC'])
    return corrected_trios


def main():
    parser = argparse.ArgumentParser(description = __doc__ )
    parser.add_argument('trio_path', help = "the path to scTRIP trios")
    parser.add_argument('bulk_scTRIP', help = 'the path to bulk TRIP bcs')
    parser.add_argument('prefix', help = 'prefix to the output file', default = 'scTRIP')
    args = parser.parse_args()
    corrected_bulk_trios = error_correction_with_bulk(
        args.trio_path, args.bulk_scTRIP
    )
    corrected_bulk_trios.to_csv(str(args.prefix) + '.tsv', sep = '\t', index = False)
if __name__ == '__main__':
    main()
