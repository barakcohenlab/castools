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
        'rep1',
        help="File for replicate1"
    )
    parser.add_argument(
        'rep2',
        help="File for replicate2"
    )
    parser.add_argument(
        'output',
        help="Prefix for output files"
    )
    args = parser.parse_args()
    return args

def read_sc(rep_f):
    rep = pd.read_csv(rep_f, sep = "\t")
    return rep

def plot_rep1_rep2(rep1, rep2, prefix):
    merged = rep1.merge(rep2, how = 'inner', on = 'tBC')
    print(merged.shape)
    plt.scatter(merged['mu_x'], merged['mu_y'])
    plt.savefig(prefix + "_bulk_sc_mean.png")
    merged.to_csv(prefix + "_merged.tsv", sep = "\t")
    from sklearn.linear_model import LinearRegression
    model = LinearRegression(fit_intercept = False)
    #print(np.asarray(merged['mu']).reshape(-1, 1))
    #print(np.asarray(merged['exp']))
    reg = model.fit(np.asarray(merged['mu_x']).reshape(-1, 1), np.asarray(merged['mu_y']))
    print(np.sqrt(reg.score(np.asarray(merged['mu_x']).reshape(-1, 1), np.asarray(merged['mu_y']))))
    print(reg.score(np.asarray(merged['mu_x']).reshape(-1, 1), np.asarray(merged['mu_y'])))

def main():
    args = parse_arguments()
    rep1 = read_sc(args.rep1)
    rep2 = read_sc(args.rep2)
    plot_rep1_rep2(rep1, rep2, args.output)

if __name__ == "__main__":
    main()
