import os,sys
import logger
import pickle
# Utility Package
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import networkx as nx
import time 

def get_clone_info_fast(df, min_cell_overlap):
    total_bcs = list(set(df['tBC'].to_list()))
    # Remove LP1 barcode 
    if 'AGGTTGCACGACAATC' in total_bcs:
        total_bcs.remove('AGGTTGCACGACAATC')
    pop_list = []
   # a loop structure as the number of pop_dict stops to increase
    cluster_num = 1
    start = time.time()
    for idx, tBC in enumerate(total_bcs):
        end = time.time()
        percent = round(idx/len(total_bcs)*100, 2)
        if idx % 50 == 0:
            print(f'we are dealing with {percent}%, time has ellapsed {round(end-start,2 )}')
        connected_list = check_cell_overlap_known_v2(df, tBC, min_cell_overlap)
        connected_list.append(tBC)
        pop_list.append(connected_list)
    return pop_list

def check_cell_overlap_known_v2(df, target_BC, min_cell_overlap):
    tBC_list = list(set(df['tBC']))
    slice_df = df[df['tBC']==target_BC]
    target_cell_list = set(slice_df['cellBC'].to_list())
    # Now pop the old barcode
    tBC_list.remove(target_BC)
    connectome = {}
    for _, row in df.iterrows():
        tBC = row['tBC']
        cellBC = row['cellBC']
        # Check the logic here 
        if tBC not in connectome:
            connectome[tBC] = []
            if cellBC in target_cell_list:
                connectome[tBC].append(cellBC)
        elif tBC in connectome:
            if cellBC in target_cell_list:
                connectome[tBC].append(cellBC)
    # Counting
    tBC_len_dic = {k:len(v) for k, v in connectome.items()}
    # return keys with more than 20 connections
    connected_bc = [k for k, v in tBC_len_dic.items() if v > min_cell_overlap]
    return connected_bc

def merge_network(coll):
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

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        'trio',
        help="Path to the TRIOs"
    )
    parser.add_argument(
        'min_cell',
        help=" minimum number of cells that overlaps between two tripBCs",
        default = 7
    )
    parser.add_argument(
            'file_name',
            help=" minimum number of cells that overlaps between two tripBCs",
        )

    args = parser.parse_args()
    tBC_lists= get_clone_info_fast(args.trio)
    final_network = merge_network(tBC_lists)
    with open(arg.file_name + '.pkl', 'wb') as f:
        pickle.dump(final_network, f)

if __name__ == "__main__":
    logger = logging.getLogger('network_analysis_log')
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
