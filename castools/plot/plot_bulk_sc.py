import gzip
import sys
import argparse
import csv
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import scipy.stats as stats
from sklearn.linear_model import LinearRegression

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
    merged["ncells_log10"] = np.log10(merged["ncells"])
    print(merged.shape)
    print(merged.head)
    plt.figure(0)
    sns.scatterplot(data = merged, x = 'mu', y = 'exp', hue = 'ncells_log10')
    plt.xlabel("single-cell mu")
    plt.ylabel("bulk exp")
    model = LinearRegression(fit_intercept = True)
    reg = model.fit(np.asarray(merged['mu']).reshape(-1, 1), np.asarray(merged['exp']))
    #plt.title("correlation-mu-unfiltered" + " " + prefix + " " +  str(np.sqrt(reg.score(np.asarray(merged['mu']).reshape(-1, 1), np.asarray(merged['exp'])))))
    plt.title("correlation-mu-unfiltered" + str(stats.pearsonr(merged['mu'], merged['exp'])))
    sns.scatterplot(data = merged, x = 'mu', y = 'exp', hue = 'ncells_log10')
    plt.savefig(prefix + "_bulk_sc_mu.png")
    merged.to_csv(prefix + "_bulk_sc_merged.tsv", sep = "\t")
    print("mu")
    print("correlation-mu-unfiltered", prefix, np.sqrt(reg.score(np.asarray(merged['mu']).reshape(-1, 1), np.asarray(merged['exp']))))
    print("correlation-mu-unfiltered", prefix, stats.pearsonr(merged['mu'], merged['exp'])[0])

    plt.figure(1)
    plt.xlabel("single-cell mean")
    plt.ylabel("bulk exp")
    model = LinearRegression(fit_intercept = True)
    reg = model.fit(np.asarray(merged['mean']).reshape(-1, 1), np.asarray(merged['exp']))
    #plt.title("correlation-mean-unfiltered" + " " + prefix + " " + str(np.sqrt(reg.score(np.asarray(merged['mean']).reshape(-1, 1), np.asarray(merged['exp'])))))
    plt.title("correlation-mean-unfiltered" + str(stats.pearsonr(merged['mean'], merged['exp'])))
    sns.scatterplot(data = merged, x = 'mean', y = 'exp', hue = 'ncells_log10')
    plt.savefig(prefix + "_bulk_sc_mean.png")
    print("mean")
    print("correlation-mean-unfiltered", prefix, np.sqrt(reg.score(np.asarray(merged['mean']).reshape(-1, 1), np.asarray(merged['exp']))))
    print("correlation-mean-unfiltered", prefix, stats.pearsonr(merged['mean'], merged['exp'])[0])

def main():
    args = parse_arguments()
    bulk = read_bulk(args.bulk)
    sc = read_sc(args.sc)
    plot_bulk_sc(bulk, sc, args.output)

if __name__ == "__main__":
    main()
