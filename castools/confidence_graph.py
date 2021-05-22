'''
A script to generate and save a graph representation of all the error corrected TripBCs.
Input: error-corrected trios (both cellBC and tripBC)
Output: a pickle file contains a graph object
---------------
Todo: identify densely connected nodes to determine the original cell with groups of tripBCs

'''
import networkx as nx
import pandas as pd
import argparse
from datetime import datetime
import pickle

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

def build_graph(trio_path, file_name):
    corrected_trio = read_in_trio(trio_path)
    tripbc_dict = hash_BCs(corrected_trio)
    with open(file_name + 'tripbc_hash' + '.pickle', 'rb') as handle:
        pickle.dump(tripbc_dict, handle, protocol = pickle.HIGHEST_PROTOCOL)
    graph = generate_graph(corrected_trio, tripbc_dict)
    with open(file_name + '.pickle', 'rb') as handle:
        pickle.dump(graph, handle, protocol = pickle.HIGHEST_PROTOCOL)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        'trio',
        help="Path to the TRIOs"
    )
    parser.add_argument(
        'file_name',
        help=" FileName that you want to save the file with, no suffix",
        default = f'{datetime.now().strftime("%d/%m/%Y %H:%M:%S")} confidence_gradph'
    )
    args = parser.parse_args()
    build_graph(args.trio, args.file_name)


if __name__ == "__main__":
    main()