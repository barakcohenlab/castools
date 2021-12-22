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
        help="File with sc cellbc UMI tripbc read-count"
    )
    parser.add_argument(
        'bulk',
        help="File with bulk means"
    )
    parser.add_argument(
        'output',
        help="Prefix for output files"
    )
    args = parser.parse_args()
    return args
def main():
    args = parse_arguments()
    bulk = read_bulk(args.bulk)
    sc = read_sc(args.sc)
    plot_bulk_sc(bulk, sc, args.output)

def read_bulk(bulk_f):
    bulk_df = pd.read_csv(bulk_f, sep = "\t")
    return bulk_df

def compute_sc_means_like_bulk(sc_f, bulk):
    """
    Sum up the number of reads at each Trip barcode and normalize by number of
    cell barcodes.
    """
    trip_nreads = {} # Number of reads per TRIP
    trip_cells = {} # list of cells per trip
    trip_means = {} # Mean number of reads per TRIP (reads normalized by number of cells)
    trip_dna_counts = pd.read_csv("CAS_LP1_TRIP_DNA_counts.tsv", delimiter = "\t", index_col = "tBC")
    sc_means = pd.DataFrame()
    bulk_tBC = bulk['tBC'].to_list()
    with open(sc_f) as sc_fh:
        for line in sc_fh:
            cell, umi, trip, count = line.split("\t")
            if trip not in trip_nreads:
                trip_nreads[trip] = 0
                trip_cells[trip] = set()
            trip_nreads[trip] += int(count)
            trip_cells[trip].add(cell)
    for tBC in trip_nreads:
        if tBC in bulk_tBC:
            #print(tBC, trip_nreads[tBC], trip_nreads[tBC]/len(trip_cells[tBC]))
            trip_means[tBC] = trip_nreads[tBC]/len(trip_cells[tBC])
            sc_means = sc_means.append({"tBC": tBC, "mean": trip_nreads[tBC]/len(trip_cells[tBC]), "mean2": trip_nreads[tBC]/trip_dna_counts.loc[tBC, "norm_count"]}, ignore_index = True)
    return(sc_means)

def plot_bulk_sc(bulk_df, sc_df, prefix):
    from sklearn.linear_model import LinearRegression
    merged = bulk_df.merge(sc_df, how = 'inner', on = 'tBC')
    print(merged.shape)
    print(merged.head(2))
    plt.figure(1)
    plt.scatter(merged['mean'], merged['exp'])
    plt.xlabel("sc mean read-count")
    plt.ylabel("bulk exp")
    plt.savefig(prefix + "_bulk_sc_meanreadcount.png")
    plt.figure(2)
    plt.scatter(merged['mean2'], merged['exp'])
    plt.xlabel("sc mean read-count normalized by bulk DNA counts")
    plt.ylabel("bulk exp")
    plt.savefig(prefix + "_bulk_sc_meanreadcount-2.png")
    print("mu R^2")
    model = LinearRegression(fit_intercept = False)
    reg = model.fit(np.asarray(merged['mean']).reshape(-1, 1), np.asarray(merged['exp']))
    print(reg.score(np.asarray(merged['mean']).reshape(-1, 1), np.asarray(merged['exp'])))
    print("mean pearson")
    print(np.corrcoef(np.asarray(merged['mean']), np.asarray(merged['exp'])))

def main():
    args = parse_arguments()
    bulk = read_bulk(args.bulk)
    sc_means = compute_sc_means_like_bulk(args.sc, bulk)
    plot_bulk_sc(bulk, sc_means, args.output)

if __name__ == "__main__":
    main()
