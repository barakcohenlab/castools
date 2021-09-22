'''
Scratch for MASCOT.
First just write all the utility functions. 
'''
import networkx as nx
import os, sys
from networkx.algorithms.components import connected
import numpy as np
import pandas as pd
import random
import os
from collections import Counter
from collections import OrderedDict
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


def return_tripBC(path):
    tripBC_list = []
    with open(path) as cellumis_tripbc_fh:
        for line in cellumis_tripbc_fh:
            trio = line.rstrip("\n").split('\t')
            try:
                if len(trio) == 4:
                    _,_,tripBC,_ = trio
            except ValueError: 
                print('Not getting the correct input')   
            if len(tripBC) == 16:         
                tripBC_list.append(tripBC)
            else:
                print('wrong bc:' ,tripBC)
    tripBC_set = tripBC_list
    return tripBC_set

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
            trio = line.rstrip("\n").split('\t')
            try:
                if len(trio) == 4:
                    _,_,tripBC,_ = trio
                elif len(trio) == 3:
                    _,_,tripBC = trio
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
        if collapsed_tripBCs[trip_bc] != None:
            new_tripBC = collapsed_tripBCs[trip_bc]
        else:
            new_tripBC = trip_bc
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
###########################
########################### Bulk Error Correction ##############################
##########################

########################## Check hamming distribution
# Considerations for this analysis see notebook: CAS April 12th 2021
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

def check_bulk_hamming(bulk_bc_list, trio, trial = 10000):
    hamming_list = []
    for bulk_bc in bulk_bc_list:
        sc_tripBC_list = random.sample(trio, trial)
        for sc_tripBC in sc_tripBC_list:
            hamming_list.append(hammingDist(bulk_bc, sc_tripBC))
    return hamming_list

########################### Error Correction using bulk RNA-seq ################
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
def error_correction_with_bulk(trios_path, bulk_path):
    '''
    Right now it is a messsssssssss!
    '''
    bulk_bc_list = read_bulk_bc(bulk_path)
    print(f'the total number of bulk TRIP BC is {len(bulk_bc_list)}')
    trios = read_in_trio_list(trios_path)
    sc_tripBC_list = [a[2] for a in trios]
    print('we are going to correct sc tripBC')
    corrected_sc_tripBC_dict = bulk_error_correction(sc_tripBC_list, bulk_bc_list)
    corrected_trios = write_correct_tripBCs(trios, corrected_sc_tripBC_dict)
    corrected_trios = pd.DataFrame(corrected_trios, columns= ['cellBC', 'umi', 'tripBC'])
    corrected_trios = remove_ambiguous_tripBC(corrected_trios,bulk_bc_list)
    return corrected_trios

###########################
# additional code 20210607
###########################

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
            trio = line.rstrip("\n").split('\t')
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

def generate_graph(trio, tripbc_dict):
    '''
    Function to create a graph for the barcodes.
    Think about OOP in next iteration
    '''
    #Initiate graph
    G = nx.Graph()
    # Return the list of tripBC values
    nodes = list(tripbc_dict.values())
    # Initiate Nodes
    G.add_nodes_from(nodes)
    # Create edges
    cellBC_list = list(set(trio['cellBC'].values))
    # Create a dictionary to hold edges and weights
    # TODO: break this up into smaller functional units
    edge_dict = {}
    # This part is really really slow, why?
    for cellBC in cellBC_list:
        same_cell_data = trio.loc[trio['cellBC'] == cellBC]
        # TODO: this is not optimal, I am basically rewriting a graph
        # Construction again, need to be better,
        tripBC_list = list(set(same_cell_data['tripBC'].values))
        # If there is no tripBC?
        # Trick is to drop the used tripBC
        not_look_list = []
        if len(tripBC_list) > 0:
            tripBC = tripBC_list[0]
            not_look_list.append(tripBC)
            tripBC_to_connect = [n for n in tripBC_list if n != tripBC]
            if tripBC not in edge_dict.keys():
                # Initiate another dict
                edge_dict[tripBC] = {}
                # Add all other tripBC as connected 
                for connected_tripBC in tripBC_to_connect:
                    # Slice the df with same cellBC, and the connected tripBC
                    # We divide it by 2 because it's double counted
                    temp_weight = len(same_cell_data.loc[same_cell_data['tripBC'] == connected_tripBC, 'umi'])
                    edge_dict[tripBC][connected_tripBC] = temp_weight
            # Now we deal with the situation wither tripBC in already in the edge_dict
            else: 
                for connected_tripBC in tripBC_to_connect:
                    if connected_tripBC not in edge_dict[tripBC].keys():
                        # Slice the df to add new key to the weight
                        temp_weight = len(same_cell_data.loc[same_cell_data['tripBC'] == connected_tripBC, 'umi'])
                        edge_dict[tripBC][connected_tripBC] = temp_weight                        
                    else:
                        temp_weight = len(same_cell_data.loc[same_cell_data['tripBC'] == connected_tripBC, 'umi'])
                        # Update the weights
                        edge_dict[tripBC][connected_tripBC] += temp_weight
                   
    #### Now we deal with the edge dict
    for start_node in edge_dict.keys():
        end_node_dict = edge_dict[start_node]
        # Get the numerical 
        hashed_start = tripbc_dict[start_node]
        for end_node in end_node_dict.keys():
            hashed_end = tripbc_dict[end_node]
            weight = end_node_dict[end_node]
            # Add edge to the graph
            G.add_edge(hashed_start, hashed_end, weight = weight)
    return G



    

######################### Graph Operation ######################################

# Look at the distributions and maintain true connections 

# Empirical CDF

def ecdf(data):
    """ https://cmdlinetips.com/2019/05/empirical-cumulative-distribution-function-ecdf-in-python/ """
    x = np.sort(data)
    n = x.size
    y = np.arange(1, n+1) / n
    return(x,y)



########################## fast graph ###########

def generate_graph_fast(trio, tripbc_dict):
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
    # Drop umi
    duo = trio.drop(['umi'], axis = 1)
    # print total length
    print(f'the total length of duo is {len(duo)}')
    # Get unique cell --> tripBC counts
    duo = duo.drop_duplicates()
    print(f'the unique duos have {len(duo)} members')
    for pos, start_tripBC in enumerate(tripBC_list):
        temp_tripBC_list = [x for x in tripBC_list if x != start_tripBC]
        print(f'We are on {pos} of {len(tripBC_list)} start Node')
        edge_dict[start_tripBC] = {}
        for end_tripBC in temp_tripBC_list:
            cellBC_start = set(duo.loc[duo['tripBC'] == start_tripBC]['cellBC'].values)
            cellBC_end = set(duo.loc[duo['tripBC'] == end_tripBC]['cellBC'].values)
            weight = len(cellBC_start.intersection(cellBC_end))
            edge_dict[start_tripBC][end_tripBC] = weight
        #### Now we deal with the edge dict
    for start_node in edge_dict.keys():
        end_node_dict = edge_dict[start_node]
        # Get the numerical 
        hashed_start = tripbc_dict[start_node]
        for end_node in end_node_dict.keys():
            hashed_end = tripbc_dict[end_node]
            weight = end_node_dict[end_node]
            # Add edge to the graph
            G.add_edge(hashed_start, hashed_end, weight = weight)
    return G



# The third iteration of the graph generation, this time I would like to add a 

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
    print(f'the unique duos have {len(trio)} members')
    for pos, start_tripBC in enumerate(tripBC_list):
        temp_tripBC_list = [x for x in tripBC_list if x != start_tripBC]
        print(f'We are on {pos} of {len(tripBC_list)} start Node')
        edge_dict[start_tripBC] = {}
        for end_tripBC in temp_tripBC_list:
            weight = count_umi_weight_helper(trio, start_tripBC, end_tripBC)
            edge_dict[start_tripBC][end_tripBC] = weight
        #### Now we deal with the edge dict
    for start_node in edge_dict.keys():
        end_node_dict = edge_dict[start_node]
        # Get the numerical 
        hashed_start = tripbc_dict[start_node]
        for end_node in end_node_dict.keys():
            hashed_end = tripbc_dict[end_node]
            weight = end_node_dict[end_node]
            # Add edge to the graph
            G.add_edge(hashed_start, hashed_end, weight = weight)
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
    weight = len(start_trio_df[start_trio_df['cellBC'].isin(common_cellBC)]) + len(end_trio_df[end_trio_df['cellBC'].isin(common_cellBC)])
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
            weight += len(start_cell_bc) 
            weight += len(end_cell_bc)
    return weight


def pre_filter(trios, min_umi_per_cell = 25, max_tripBC_per_cell = 100):
    '''
    Input: trios: containing cellBC, umi, tripBC.
        min_umi_per_cell: the minimum number of cells that are used for 
        max_tripBC_per_cell: the maximum number of tripBC that a cell can have. 
    Output: new trios that all the filters are done.
    '''
    # 
    cell_bc_list = list(set(trios['cellBC'].values))
    filtered_cellBC_list = []
    for cell_bc in cell_bc_list:
        data_slice = trios[trios['cellBC'] == cell_bc]
        print(f'we are dealing with {cell_bc}')
        if len(data_slice) > min_umi_per_cell:
            if len(set(data_slice['tripBC'])) < max_tripBC_per_cell:
                filtered_cellBC_list.append(cell_bc)
    # filter based on min umi per cell
    umi_filtered_trio = trios[trios.cellBC.isin(cell_bc_list)]
    return umi_filtered_trio
