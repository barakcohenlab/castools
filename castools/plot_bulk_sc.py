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
        'bulk',
        help="File with bulk means"
    )
    parser.add_argument(
        'sc',
        help="File with sc means"
    )
    parser.add_argument(
        'output',
        help="Prefix for output files"
    )
    args = parser.parse_args()
    return args

def read_bulk(bulk_f):
    bulk_df = pd.read_csv(bulk_f, sep = "\t")
    return bulk_df

def read_sc(sc_f):
    sc = pd.read_csv(sc_f, sep = "\t")
    return sc

def plot_bulk_sc(bulk_df, sc_df, prefix):
    merged = bulk_df.merge(sc_df, how = 'inner', on = 'tBC')
    print(merged.shape)
    plt.scatter(merged['mu'], merged['exp'])
    plt.savefig(prefix + "_bulk_sc_mean.png")
    merged.to_csv(prefix + "_merged.tsv", sep = "\t")
    from sklearn.linear_model import LinearRegression
    model = LinearRegression(fit_intercept = False)
    #print(np.asarray(merged['mu']).reshape(-1, 1))
    #print(np.asarray(merged['exp']))
    reg = model.fit(np.asarray(merged['mu']).reshape(-1, 1), np.asarray(merged['exp']))
    print(np.sqrt(reg.score(np.asarray(merged['mu']).reshape(-1, 1), np.asarray(merged['exp']))))

def main():
    args = parse_arguments()
    bulk = read_bulk(args.bulk)
    sc = read_sc(args.sc)
    print(bulk.head(), sc.head())
    plot_bulk_sc(bulk, sc, args.output)

if __name__ == "__main__":
    main()
