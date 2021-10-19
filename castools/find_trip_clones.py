'''
Input: fully filtered trio -- csv or tsv file with trio and counts; minimum weight -- empirically check, int; file-name --- to save the network
Output: a pickle file contains the list of list of clonal tripBCs
'''

import os,sys
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import networkx as nx
import argparse

def get_clone_info_fast(df, min_weight =10):
    
    total_bcs = list(set(df['tBC'].to_list()))
    # This is hardcoded LP barcode
    total_bcs.remove('AGGTTGCACGACAATC')
    pop_list = []
   # a loop structure as the number of pop_dict stops to increase
    cluster_num = 1
    for idx, tBC in enumerate(total_bcs):
        print(f'we are dealing with {idx}')
        connected_list = check_cell_overlap_known(df, tBC, min_weight)
        connected_list.append(tBC)
        pop_list.append(connected_list)
    return pop_list

def check_cell_overlap_known(df, target_BC, mininum_weight):
    tBC_list = list(set(df['tBC']))
    slice_df = df[df['tBC']==target_BC]
    target_cell_list = set(slice_df['cellBC'].to_list())
    # Now pop the old barcode
    tBC_list.remove(target_BC)
    connected_bc = []
    for idx, test_BC in enumerate(tBC_list):
        slice_df = df[df['tBC']==test_BC]
        #slice_df = slice_df.loc[slice_df['count'] > 10]
        test_cell_list = set(slice_df['cellBC'].to_list())
        overlap = len(target_cell_list.intersection(test_cell_list))
        if overlap> mininum_weight:
            connected_bc.append(test_BC)
    return connected_bc

def collapse_network(coll):
    coll = list_of_b
    edges = []
    for i in range(len(coll)):
        a = coll[i]
        for j in range(len(coll)):
            if i != j:
                b = coll[j]
                if set(a).intersection(set(b)):
                    edges.append((i,j))

    G = nx.Graph()
    G.add_nodes_from(range(len(coll)))
    G.add_edges_from(edges)
    final_network = []
    for c in nx.connected_components(G):
        combined_lists = [coll[i] for i in c]
        flat_list = [item for sublist in combined_lists for item in sublist]
        final_network.append(list(set(flat_list)))
    return final_network


def parse_arguments():
    parser = argparse.ArgumentParser(description='Save a pickle file of a list of list each sublist is a group of tripBCs')
    parser.add_argument("--trio",
                        help = "Trio file contains the counts")
    parser.add_argument("--min_weight",
                    help = "minimum weight required to call share bcs")
    parser.add_argument("--file_name",
                    help = "name of the pickle file to save with")
    args = parser.parse_args()
    return args

def main():
    args = parse_arguments()
    print('Read in arguments fine')
    if args.trio.endswith('.tsv'):
        trio = pd.read_csv('args.trio', sep = '\t')
    elif args.trio.endswith('.csv'):
        trio = pd.read_csv('args.trio', sep = '\t')
    # Construct the 1 to all network 
    net = get_clone_info_fast(trio, args.min_weight)
    # Copplase network
    final_network = collapse_network(net)
    # Save it with pickle
    with open(args.filename + '.pkl', 'wb') as handle:
        pickle.dump(final_network, handle, protocol=pickle.HIGHEST_PROTOCOL)