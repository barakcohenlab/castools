#!/bin/bash/python3

import sys
import argparse
import matplotlib.pyplot as plt
import seaborn as sns
from collections import defaultdict
import numpy as np
import itertools

parser = argparse.ArgumentParser(
    description = 'Plot whether the same cell bc/UMIs have the same TRIP bcs')
parser.add_argument(
    'bc_file', help = 'Cell bc-UMI-TRIP bc file with counts')
parser.add_argument(
    'output', help = 'Output figure filename')
args = parser.parse_args()

# initialise dictionary to store cell_bc/UMIs as keys and TRIP bcs as values
cell_bcs_umi = defaultdict(dict)
total_pairs = 0

# read the barcode file
with open(args.bc_file, 'r') as f:
    for line in f:
        cell_bc, umi, trip_bc, count = line.strip('\n').split('\t')
        cell_bcs_umi[(cell_bc, umi)][trip_bc] = int(count)
        total_pairs += 1

def hamdist(str1, str2): # From http://code.activestate.com/recipes/499304-hamming-distance/
    """Count the # of differences between equal length strings str1 and str2"""

    diffs = 0

    for ch1, ch2 in zip(str1, str2):
        if ch1 != ch2:
            diffs += 1
            
    return diffs

def calculate_pairwise_hamming_distances(trip_bcs):
    """ Calculate pairwise hamming distances between list of TRIP bcs"""

    hamming_dists = []

    for bc1, bc2 in itertools.combinations(trip_bcs, 2):
        hamming_dists.append(hamdist(bc1, bc2))

    return hamming_dists

def check_top_barcode(trip_bc_counts):
    """Calculate the proportion of reads the top barcodes takes"""

    max_count = 0

    for bc, count in trip_bc_counts.items():
        if count > max_count:
            max_count = count
            top_bc = bc

    other_counts = 0

    for bc, count in trip_bc_counts.items():
        if bc != top_bc:
            other_counts += count

    return trip_bc_counts[top_bc]/other_counts 

def process_multiple_trip_bcs(trip_bc_counts):
    """Process cell_bc/UMI pairs that have multiple TRIP barcodes assigned"""

    if check_top_barcode(trip_bc_counts) > 0.9:
        return 'mostly_matched'
    else:
        trip_bcs = [k for k, v in trip_bc_counts.items()]
        hamming_dists = calculate_pairwise_hamming_distances(trip_bcs)
        if np.mean(hamming_dists) < 2:
            return 'mostly_matched'
        else:
            return 'not_matched'

# initialise counting of cell_bc/UMIs that match one or multiple TRIP bcs
one_count = 0
mostly_matched = 0
not_matched = 0

# count number of cell_bc/umi pairs that have one vs multiple TRIP bcs
for bc_pair, trip_bcs in cell_bcs_umi.items():
    if len(trip_bcs) == 1:
        one_count += 1
    else:
        matched = process_multiple_trip_bcs(trip_bcs)
        if matched == 'mostly_matched':
            mostly_matched += 1
        elif matched == 'not_matched':
            not_matched += 1

# plot the results
sns.barplot(x = ['one', 'mostly_matched', 'multiple'],
            y = [one_count, mostly_matched, not_matched])
plt.xlabel('Number of TRIP barcodes')
plt.ylabel('Count')
plt.savefig(args.output + '.png')

    