import gzip
import sys
import argparse
import csv
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.linear_model import LinearRegression

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

    #Print R^2, the score from the regression
    plt.figure(0)
    plt.scatter(merged['mu_x'], merged['mu_y'])
    plt.xlabel("rep1 mu")
    plt.ylabel("rep2 mu")
    plt.savefig(prefix + "_bulk_sc_mu.png")
    merged.to_csv(prefix + "_merged.tsv", sep = "\t")
    model = LinearRegression(fit_intercept = False)
    reg = model.fit(np.asarray(merged['mu_x']).reshape(-1, 1), np.asarray(merged['mu_y']))
    print("mu")
    print(reg.score(np.asarray(merged['mu_x']).reshape(-1, 1), np.asarray(merged['mu_y'])))
    plt.close()

    plt.figure(1)
    plt.scatter(merged['mean_x'], merged['mean_y'])
    plt.xlabel("rep1 mean")
    plt.ylabel("rep2 mean")
    plt.savefig(prefix + "_bulk_sc_mean.png")
    model = LinearRegression(fit_intercept = False)
    reg = model.fit(np.asarray(merged['mean_x']).reshape(-1, 1), np.asarray(merged['mean_y']))
    print("mean")
    print(reg.score(np.asarray(merged['mean_x']).reshape(-1, 1), np.asarray(merged['mean_y'])))
    plt.close()

    plt.figure(2)
    plt.scatter(merged['var_x'], merged['var_y'])
    plt.xlabel("rep1 var")
    plt.ylabel("rep2 var")
    plt.savefig(prefix + "_bulk_sc_var.png")
    model = LinearRegression(fit_intercept = True)
    reg = model.fit(np.asarray(merged['var_x']).reshape(-1, 1), np.asarray(merged['var_y']))
    print("var")
    print(reg.score(np.asarray(merged['var_x']).reshape(-1, 1), np.asarray(merged['var_y'])))
    plt.close()

    plt.figure(3)
    plt.scatter(merged['alpha_x'], merged['alpha_y'])
    plt.xlabel("rep1 alpha")
    plt.ylabel("rep2 alpha")
    plt.savefig(prefix + "_bulk_sc_alpha.png")
    model = LinearRegression(fit_intercept = True)
    reg = model.fit(np.asarray(merged['alpha_x']).reshape(-1, 1), np.asarray(merged['alpha_y']))
    print("alpha")
    print(reg.score(np.asarray(merged['alpha_x']).reshape(-1, 1), np.asarray(merged['alpha_y'])))
    plt.close()

    plt.figure(4)
    plt.scatter(merged['mu_x'], merged['alpha_x'])
    plt.xlabel("rep1 mu")
    plt.ylabel("rep1 alpha")
    plt.savefig(prefix + "_sc_mu_alpha_x.png")
    model = LinearRegression(fit_intercept = False)
    reg = model.fit(np.asarray(merged['mu_x']).reshape(-1, 1), np.asarray(merged['alpha_x']))
    print("mu vs alpha x")
    print(reg.score(np.asarray(merged['mu_x']).reshape(-1, 1), np.asarray(merged['alpha_x'])))
    plt.close()

    plt.figure(5)
    plt.scatter(merged['mu_y'], merged['alpha_y'])
    plt.xlabel("rep1 mu")
    plt.ylabel("rep1 alpha")
    plt.savefig(prefix + "_sc_mu_alpha_y.png")
    model = LinearRegression(fit_intercept = False)
    reg = model.fit(np.asarray(merged['mu_y']).reshape(-1, 1), np.asarray(merged['alpha_y']))
    print("mu vs alpha y")
    print(reg.score(np.asarray(merged['mu_y']).reshape(-1, 1), np.asarray(merged['alpha_y'])))
    plt.close()

def main():
    args = parse_arguments()
    rep1 = read_sc(args.rep1)
    rep2 = read_sc(args.rep2)
    plot_rep1_rep2(rep1, rep2, args.output)

if __name__ == "__main__":
    main()
