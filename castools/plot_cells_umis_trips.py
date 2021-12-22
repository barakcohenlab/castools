import gzip
import sys
import argparse
import csv
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

def parse_arguments():
    parser = argparse.ArgumentParser(description="Plot bulk means against sc means")
    parser.add_argument(
        'sc',
        help="File with sc cells UMIs and TRIPs"
    )
    parser.add_argument(
        'output',
        help="Prefix for output files"
    )
    args = parser.parse_args()
    return args


def read_sc(sc_f):
    cell_umis = {}
    cell_trips = {}
    with open(sc_f) as sc_fh:
        header = True
        for line in sc_fh:
            if header:
                header = False
                next
            cell, umi, trip = line.split("\t")
            if cell not in cell_umis:
                cell_umis[cell] = set()
                cell_trips[cell] = set()
            cell_umis[cell].add(umi)
            cell_trips[cell].add(trip)
    return {"cell_umis": cell_umis, "cell_trips": cell_trips}

def plot_sc(sc, prefix):
    cell_umis_lengths = [len(sc["cell_umis"][x]) for x in sc["cell_umis"]]
    cell_trips_lengths = [len(sc["cell_trips"][x]) for x in sc["cell_trips"]]
    plt.figure(1)
    plt.hist(np.log10(cell_umis_lengths), bins = 30)
    plt.xlabel("Number of UMIs per cell (log10)")
    plt.savefig(prefix + "_umis_per_cell.png")
    plt.figure(2)
    plt.hist(np.log10(cell_trips_lengths), bins = 30)
    plt.xlabel("Number of TRIPs per cell (log10)")
    plt.savefig(prefix + "_trips_per_cell.png")
    plt.figure(3)
    plt.bar(range(len(cell_umis_lengths)), np.sort(np.log10(cell_umis_lengths)))
    plt.ylabel("Number of UMIs per cell (log10)")
    plt.savefig(prefix + "_umis_per_cell_bar.png")
    plt.figure(4)
    plt.bar(range(len(cell_umis_lengths)), np.sort(np.log10(cell_trips_lengths)))
    plt.ylabel("Number of TRIPs per cell (log10)")
    plt.savefig(prefix + "_trips_per_cell_bar.png")

def main():
    args = parse_arguments()
    sc = read_sc(args.sc)
    plot_sc(sc, args.output)

if __name__ == "__main__":
    main()
