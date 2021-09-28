import networkx as nx
import os, sys
from networkx.algorithms.components import connected
import numpy as np
import pandas as pd
import random
import pickle
from collections import Counter
from collections import OrderedDict
import matplotlib.pyplot as plt
import util

def draw_graph(graph, prefix):
    plt.figure(figsize = (12, 8))
    pos = nx.random_layout(graph)
    edges,weights = zip(*nx.get_edge_attributes(graph,'weight').items())
    nx.draw(graph, pos, node_color='k', node_size = 5,
            edgelist=edges, edge_color=weights, width=1.0, edge_cmap=plt.cm.Blues)
    plt.savefig(prefix + "_graph.png")

    plt.figure(figsize = (12, 8))
    pos = nx.random_layout(graph)
    edges,weights = zip(*nx.get_edge_attributes(graph,'weight').items())
    nx.draw_spring(graph, node_color='k', node_size = 5,
            edgelist=edges, edge_color=weights, width=1.0, edge_cmap=plt.cm.Blues)
    plt.savefig(prefix + "_spring_graph.png")

    plt.figure(figsize = (12, 8))
    plt.hist(weights, bins = 200)
    plt.xlabel('Weights')
    plt.yscale('log')
    plt.savefig(prefix + "_weights_hist.png")

    plt.figure(figsize = (12, 8))
    (x,y) = util.ecdf(weights)
    plt.scatter(x = x, y = y)
    plt.ylabel('percentage')
    plt.xlabel('edge weights')
    plt.savefig(prefix + '_weights_cdf.png')
    plt.show()

    graph_filtered = graph
    edge_weights = nx.get_edge_attributes(graph_filtered, 'weight')
    #Only keep edges with atleast weight 2
    graph_filtered.remove_edges_from((e for e, w in edge_weights.items() if w < 2))
    plt.figure(figsize = (12, 8))
    pos = nx.random_layout(graph_filtered)
    edges,weights = zip(*nx.get_edge_attributes(graph_filtered,'weight').items())
    nx.draw(graph_filtered, pos, node_color='k', node_size = 5,
            edgelist=edges, edge_color=weights, width=1.0, edge_cmap=plt.cm.Blues)
    plt.savefig(prefix + "_filtered_w2_graph.png")

    plt.figure(figsize = (12, 8))
    pos = nx.random_layout(graph_filtered)
    edges,weights = zip(*nx.get_edge_attributes(graph_filtered,'weight').items())
    nx.draw_circular(graph_filtered, node_color='k', node_size = 5,
            edgelist=edges, edge_color=weights, width=1.0, edge_cmap=plt.cm.Blues)
    plt.savefig(prefix + "_filtered_w2_graph_circular.png")
    plt.figure(figsize = (12, 8))
    pos = nx.random_layout(graph_filtered)
    edges,weights = zip(*nx.get_edge_attributes(graph_filtered,'weight').items())
    nx.draw_spectral(graph_filtered, node_color='k', node_size = 5,
            edgelist=edges, edge_color=weights, width=1.0, edge_cmap=plt.cm.Blues)
    plt.savefig(prefix + "_filtered_w2_graph_spectral.png")
    plt.figure(figsize = (12, 8))
    pos = nx.random_layout(graph_filtered)
    edges,weights = zip(*nx.get_edge_attributes(graph_filtered,'weight').items())
    nx.draw_spring(graph_filtered, node_color='k', node_size = 5,
            edgelist=edges, edge_color=weights, width=1.0, edge_cmap=plt.cm.Blues)
    plt.savefig(prefix + "_filtered_w2_graph_spring.png")
############################ 
############################ Look at the depth #################################
############################ 
def return_read_count(path):
    depth = []
    with open(path) as cellumis_tripbc_fh:
        for line in cellumis_tripbc_fh:
            trio = line.rstrip("\n").split('\t')
            try:
                if len(trio) == 4:
                    _,_,_,num_of_reads = trio
            except ValueError: 
                print('Not getting the correct input')            
            depth.append(int(num_of_reads))
    return depth


def return_tripBC_num(path):
    tripBC_list = []
    with open(path) as cellumis_tripbc_fh:
        for line in cellumis_tripbc_fh:
            trio = line.rstrip("\n").split(' ')
            try:
                if len(trio) == 4:
                    _,_,tripBC,_ = trio
            except ValueError: 
                print('Not getting the correct input')            
            tripBC_list.append(tripBC)
    tripBC_num = set(tripBC_list)
    return len(tripBC_num)

########################## tripBC read count vs occurrence. ####################

def return_tripBCs_depth():
    '''
    Does it matter? To get a sort of correlation between the sequencing depth and
    '''
    pass

############################ 
############################ Error Correction for TRIP #########################
############################ 
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


def sample_hamming_distance(path, sample_size = 100000):
    '''
    Calculating the average hamming distance among tripBCs in the file
    '''
    tripBC_list = []
    with open(path) as cellumis_tripbc_fh:
        for line in cellumis_tripbc_fh:
            trio = line.rstrip("\n").split(' ')
            try:
                if len(trio) == 4:
                    _,_,tripBC,_ = trio
            except ValueError: 
                print('Not getting the correct input') 
            if len(tripBC) == 16:           
                tripBC_list.append(tripBC)
    tripBC_set = list(set(tripBC_list))
    # Now we sample a pair from the list
    # Initiate a empty list
    hamming_sample_list = []
    j = 0
    for j in range(sample_size):
        pair = random.sample(tripBC_set,2)
        hamming = hammingDist(pair[0], pair[1])
        hamming_sample_list.append(hamming)
        j += 1
    return hamming_sample_list

# Count the number of tripBCs

def get_abundance_tripBC(path):
    tripBC_dict = {}
    with open(path) as cellumis_tripbc_fh:
        for line in cellumis_tripbc_fh:
            trio = line.rstrip("\n").split(' ')
            try:
                if len(trio) == 4:
                    _,_,tripBC,_ = trio
            except ValueError: 
                print('Not getting the correct input') 
            if len(tripBC) == 16:
                if tripBC in tripBC_dict:
                    tripBC_dict[tripBC] +=1
                else:
                    tripBC_dict[tripBC] = 1            
    return tripBC_dict
############################ How to find the nearest neighbor for the tripBCs? 
# To Be Written
############################ A not-replaced version of the code. ###############
# So the simple version of the code is to just do a irreplaceable sampling. 
def hamming_key_list(key, key_list):
    '''
    Helper function to compare the 
    '''
    for main_key in key_list:
        if hammingDist(key, main_key) <=14:
            return main_key
        else: 
            return key

def collapse_tripBC(ordered_tripBCs):
    '''
    Function to collapse the tirpBCs based on the hamming distance 
    (Possibly with the actually hamming distance distribution?)
    Right now I will just go for the empirical one
    Input: an ordered list of tripBCs based occurance. 
    Output: a dictionary with the original tripBC as the key, and the target 
        tripBC as the value
    '''
    # Now we work on getting the without replacement Hamming distance 
    # Initiate a dictionary to hold the mapping
    correspondence = {}
    # Get the keys
    search_keys = ordered_tripBCs
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
    # Use the values as keys
    collapsed_tripBCs = {}
    for key in correspondence.keys():
        new_keys = correspondence[key]
        if new_keys != None:
            for new_key in new_keys:
                collapsed_tripBCs[new_key] = key
        else:
            collapsed_tripBCs[key] = key
    return collapsed_tripBCs

def order_tripBCs(tripBCs):
    '''
    Function to return a list of ordered tripBCs
    '''
    counted_list = Counter(tripBCs)
    # Sort the items by value
    sorted_tripBCs = {k: v for k, v in sorted(counted_list.items(), key=lambda item: item[1])}
    # Convert dict to OrderedDict, while python 3.7 dictionary is now ordered, 
    # But pandas seems to ignore the order, so...
    sorted_tripBCs = OrderedDict(sorted_tripBCs)
    sorted_tripBCs_list = list(sorted_tripBCs.keys())
    return sorted_tripBCs_list

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
        new_tripBC = collapsed_tripBCs[trip_bc]
        pop_list.append([cell_bc,umi,new_tripBC])        
    return pop_list

def error_correction(path):
    trios = []
    tripBCs = []
    # Read in file
    with open(path) as f:
        for line in f:
            trio = line.rstrip("\n").split(' ')
            tripBC = trio[2]
            trios.append(trio[0:3])
            tripBCs.append(tripBC)
    # Order the trip BCs and return the ordered list (unique tirpBCs)
    ordered_tripBCs = order_tripBCs(tripBCs)
    # Error correction without replacement 
    collapsed_tripBCs = collapse_tripBC(ordered_tripBCs)
    # Write the correct tripBCs
    corrected_trios = write_correct_tripBCs(trios, collapsed_tripBCs)
    corrected_trios = pd.DataFrame(corrected_trios, columns= ['cellBC', 'umi', 'tripBC'])
    return corrected_trios

def remove_ambiguous_tripBC(corrected_trios, bulk_bc_list):
    '''
    Helper function to remove ambiguous trip barcodes 
    ''' 
    pop_df = corrected_trios[corrected_trios.tripBC.isin(bulk_bc_list)]
    return pop_df

############################ 
############################ Confidence Graph ##################################
############################ 
def read_in_trio(path):
    tripData = []
    with open(path) as cellumis_tripbc_fh:
        for line in cellumis_tripbc_fh:
            trio = line.rstrip("\n").split()
            if len(trio) == 4:
                cell,umi,trip_bc, _ = trio
            else: 
                cell,umi,trip_bc = trio                
            tripData.append([cell,umi,trip_bc])
    trio_df = pd.DataFrame(tripData, columns= ['cellBC', 'umi', 'tripBC'])
    return trio_df

def hash_BCs(trio):
    '''
    Go through all the trios and hash all the tripBCs.
    
    '''
    tripbc_list = list(set(trio['tripBC'].values))
    tripbc_dict = {}
    for num, bc in enumerate(tripbc_list):
        tripbc_dict[bc] = num
    return tripbc_dict

def read_bulk_bc(path):
    print(os.path.splitext(path)[1])
    if os.path.splitext(path)[1] == '.txt':
        bulk_bc_list = []
        with open(path) as f:
            next(f)
            for line in f:
                bulk_bc = line.rstrip('\n')
                if len(bulk_bc) == 16:
                    bulk_bc_list.append(bulk_bc)
                else:
                    print('not 16:', bulk_bc)
    elif os.path.splitext(path)[1] == '.tsv':
        bulk_bc_df = pd.read_csv(path, sep = '\t')
        bulk_bc_list = list(bulk_bc_df['tBC'].values)
    return bulk_bc_list

########################## UMI Graph ###########

def generate_graph_count_umi(trio, tripbc_dict):
    '''
    Function to create a graph for the barcodes.
    Think about OOP in next iteration

    '''
    #Initiate graph
    G = nx.Graph()
    # Return the list of tripBC values
    nodes = list(tripbc_dict.values())
    tripBC_list = list(tripbc_dict.keys())
    # Initiate Nodes
    G.add_nodes_from(nodes)
    # Create a dictionary to hold edges and weights
    # TODO: break this up into smaller functional units
    edge_dict = {}
    ## Construct a list of beginTRIP:{endtripDictList} where endTRIP:weights
    # Do not drop umi
    print(f'the unique duos have {len(trio)} members', file = sys.stderr)
    print(f'the number of trip barcodes is {len(tripBC_list)} members', file = sys.stderr)
    for pos, start_tripBC in enumerate(tripBC_list):
        print(f'We are on {pos + 1} of {len(tripBC_list)} start Node', file = sys.stderr)
        edge_dict[start_tripBC] = {}
        for end_tripBC in tripBC_list[pos + 1:]: # Compute only top half of the matrix since it's symmetric
            #weight = count_uni_weight_filter_helper(trio, start_tripBC, end_tripBC) # Number of UMI in shared cells
            weight = count_umi_weight_helper(trio, start_tripBC, end_tripBC) # Count number of shared cells
            edge_dict[start_tripBC][end_tripBC] = weight
        #### Now we deal with the edge dict
    for start_node in edge_dict:
        # Get the numerical
        hashed_start = tripbc_dict[start_node]
        for end_node in edge_dict[start_node]:
            hashed_end = tripbc_dict[end_node]
            weight = edge_dict[start_node][end_node]
            # Add edge to the graph
            if weight > 0:
                G.add_edge(hashed_start, hashed_end, weight = weight)
                print("Graph-out", hashed_start, hashed_end, weight, file = sys.stderr)
    return G

def count_umi_weight_helper(trio, start_tripBC, end_tripBC):
    # Slice the total df
    start_trio_df = trio.loc[trio['tripBC'] == start_tripBC]
    end_trio_df = trio.loc[trio['tripBC'] == end_tripBC]
    # Get the cellBC
    start_cellBC = set(start_trio_df['cellBC'].values)
    end_cellBC = set(end_trio_df['cellBC'].values)
    # Intersect 
    common_cellBC = start_cellBC.intersection(end_cellBC)
    #weight = len(start_trio_df[start_trio_df['cellBC'].isin(common_cellBC)]) + len(end_trio_df[end_trio_df['cellBC'].isin(common_cellBC)])
    weight = len(common_cellBC)
    return weight

def count_uni_weight_filter_helper(trio, start_tripBC, end_tripBC):
    '''
    To further filter the tripBC, I implemented a filter on the umi has to been seen at least 
    twice.
    '''
    # Slice the total df
    start_trio_df = trio.loc[trio['tripBC'] == start_tripBC]
    end_trio_df = trio.loc[trio['tripBC'] == end_tripBC]
    # Get the cellBC
    start_cellBC = set(start_trio_df['cellBC'].values)
    end_cellBC = set(end_trio_df['cellBC'].values)
    # Intersect 
    common_cellBC = start_cellBC.intersection(end_cellBC)
    weight = 0
    for cell_bc in common_cellBC:
        start_cell_bc = start_trio_df.loc[start_trio_df['cellBC'] == cell_bc]
        end_cell_bc = end_trio_df.loc[end_trio_df['cellBC'] == cell_bc]
        if (len(start_cell_bc) > 2) and (len(end_cell_bc) > 2):
            # weight += len(start_cell_bc) 
            # weight += len(end_cell_bc)
            weight = (len(start_cell_bc) + len(end_cell_bc))/2
    return weight


def pre_filter(trios, min_umi_per_cell = 25, max_tripBC_per_cell = 100):
    '''
    Input: trios: containing cellBC, umi, tripBC.
        min_umi_per_cell: the minimum number of cells that are used for 
        max_tripBC_per_cell: the maximum number of tripBC that a cell can have. 
    Output: new trios that all the filters are done.
    '''
    # 
    print("min_umi_per_cell", min_umi_per_cell, file = sys.stderr)
    print("max_tripBC_per_cell", max_tripBC_per_cell, file = sys.stderr)
    cell_bc_list = list(set(trios['cellBC'].values))
    filtered_cellBC_list = []
    for cell_bc in cell_bc_list:
        data_slice = trios[trios['cellBC'] == cell_bc]
        if len(data_slice) > min_umi_per_cell:
            if len(set(data_slice['tripBC'])) < max_tripBC_per_cell:
                filtered_cellBC_list.append(cell_bc)
    # filter based on min umi per cell
    umi_filtered_trio = trios[trios.cellBC.isin(filtered_cellBC_list)]
    return umi_filtered_trio




def main():
    trio_dir = sys.argv[1]
    # Read in the LP trios
    trio = read_in_trio(trio_dir)
    # Read in the mapped barcode list
    tripBC_list = sys.argv[2]
    op_prefix = sys.argv[3]
    tripBC = read_bulk_bc(tripBC_list)
    trio_clean = remove_ambiguous_tripBC(trio, tripBC)
    trio_filtered = pre_filter(trio_clean)
    tripbc_dict = hash_BCs(trio_filtered)
    with open(op_prefix + '_has.pickle', 'wb') as handle:
        pickle.dump(tripbc_dict, handle, protocol= pickle.HIGHEST_PROTOCOL)
    graph = generate_graph_count_umi(trio_filtered, tripbc_dict)
    with open(op_prefix + '_graph.pickle', 'wb') as handle:
        pickle.dump(graph, handle, protocol = pickle.HIGHEST_PROTOCOL)
    draw_graph(graph, op_prefix)

if __name__ == "__main__":
    main()
