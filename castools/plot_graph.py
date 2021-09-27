import argparse
import networkx as nx
import os, sys
from networkx.algorithms.components import connected
import numpy as np
import pandas as pd
import pickle
import matplotlib.pyplot as plt
import util

def parse_arguments():
    parser = argparse.ArgumentParser(description="Filter and plot the graph")
    parser.add_argument(
        'graph_pickle',
        help="Graph file that has been pickled"
    )
    parser.add_argument(
        'cutoff',
        help="Minimum number of cells supporting an edge"
    )
    parser.add_argument(
        'prefix',
        help="Prefix for output files"
    )
    args = parser.parse_args()
    return args

def main():
    args = parse_arguments()
    graph_pickle = args.graph_pickle
    cutoff = args.cutoff
    prefix = args.prefix
    with open(args.graph_pickle, "rb") as fh:
        graph = pickle.load(fh)
    plt.figure(figsize = (12, 8))

    graph_filtered = graph
    edge_weights = nx.get_edge_attributes(graph_filtered, 'weight')
    #Only keep edges with atleast weight 2
    graph_filtered.remove_edges_from((e for e, w in edge_weights.items() if w < int(cutoff)))
    plt.figure(figsize = (12, 8))
    pos = nx.random_layout(graph_filtered)
    edges,weights = zip(*nx.get_edge_attributes(graph_filtered,'weight').items())
    nx.draw(graph_filtered, pos, node_color='k', node_size = 5,
            edgelist=edges, edge_color=weights, width=1.0, edge_cmap=plt.cm.Blues)
    plt.savefig(prefix + "_filtered_w" + str(cutoff) + "_graph.png")

    plt.figure(figsize = (12, 8))
    pos = nx.random_layout(graph_filtered)
    edges,weights = zip(*nx.get_edge_attributes(graph_filtered,'weight').items())
    nx.draw_spring(graph_filtered, node_color='k', node_size = 5,
            edgelist=edges, edge_color=weights, width=1.0, edge_cmap=plt.cm.Blues)
    plt.savefig(prefix + "_filtered_w" + str(cutoff) + "_graph_spring.png")

if __name__ == "__main__":
    main()
