import gzip
import sys
import argparse
import csv
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scipy.stats as stats
from sklearn.linear_model import LinearRegression

def parse_arguments():
    parser = argparse.ArgumentParser(description="Compare replicate1 mean and variance to replicate2")
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
    merged.to_csv(prefix + "_merged.tsv", sep = "\t")
    model = LinearRegression(fit_intercept = False)
    reg = model.fit(np.asarray(merged['mu_x']).reshape(-1, 1), np.asarray(merged['mu_y']))
    #plt.title(prefix + " mu rep1 vs rep2 " + str(reg.score(np.asarray(merged['mu_x']).reshape(-1, 1), np.asarray(merged['mu_y']))))
    plt.title(prefix + " mu rep1 vs rep2 " + str(stats.pearsonr(merged['mu_x'], merged['mu_y'])))
    plt.savefig(prefix + "_rep1_rep2_mu.png")
    print("mu")
    #print(prefix + " mu rep1 vs rep2 ", reg.score(np.asarray(merged['mu_x']).reshape(-1, 1), np.asarray(merged['mu_y'])))
    print(prefix + " mu rep1 vs rep2 ", stats.pearsonr(merged['mu_x'], merged['mu_y']))
    plt.close()

    plt.figure(1)
    plt.scatter(merged['mean_x'], merged['mean_y'])
    plt.xlabel("rep1 mean")
    plt.ylabel("rep2 mean")
    model = LinearRegression(fit_intercept = False)
    reg = model.fit(np.asarray(merged['mean_x']).reshape(-1, 1), np.asarray(merged['mean_y']))
    #plt.title(prefix + " mean rep1 rep2 " + str(reg.score(np.asarray(merged['mean_x']).reshape(-1, 1), np.asarray(merged['mean_y']))))
    plt.title(prefix + " mean rep1 vs rep2 " + str(stats.pearsonr(merged['mean_x'], merged['mean_y'])))
    plt.savefig(prefix + "_rep1_rep2_mean.png")
    print("mean")
#    print(prefix + " mean rep1 rep2 ", reg.score(np.asarray(merged['mean_x']).reshape(-1, 1), np.asarray(merged['mean_y'])))
    print(prefix + " mean rep1 rep2 ", stats.pearsonr(merged['mean_x'], merged['mean_y']))
    plt.close()

    plt.figure(2)
    plt.scatter(merged['var_x'], merged['var_y'])
    plt.xlabel("rep1 var")
    plt.ylabel("rep2 var")
    model = LinearRegression(fit_intercept = True)
    reg = model.fit(np.asarray(merged['var_x']).reshape(-1, 1), np.asarray(merged['var_y']))
    #plt.title(prefix + " var rep1 rep2 " + str(reg.score(np.asarray(merged['var_x']).reshape(-1, 1), np.asarray(merged['var_y']))))
    plt.title(prefix + " var rep1 vs rep2 " + str(stats.pearsonr(merged['var_x'], merged['var_y'])))
    plt.savefig(prefix + "_rep1_rep2_var.png")
    print("var")
    #print(prefix + " var rep1 rep2 ", reg.score(np.asarray(merged['var_x']).reshape(-1, 1), np.asarray(merged['var_y'])))
    print(prefix + " var rep1 rep2 ", stats.pearsonr(merged['var_x'], merged['var_y']))
    plt.close()

    plt.figure(3)
    plt.scatter(merged['alpha_x'], merged['alpha_y'])
    plt.xlabel("rep1 alpha")
    plt.ylabel("rep2 alpha")
    model = LinearRegression(fit_intercept = True)
    reg = model.fit(np.asarray(merged['alpha_x']).reshape(-1, 1), np.asarray(merged['alpha_y']))
    #plt.title(prefix  + " alpha rep1 rep2 " + str(reg.score(np.asarray(merged['alpha_x']).reshape(-1, 1), np.asarray(merged['alpha_y']))))
    plt.title(prefix + " alpha rep1 vs rep2 " + str(stats.pearsonr(merged['alpha_x'], merged['alpha_y'])))
    plt.savefig(prefix + "_rep1_rep2_alpha.png")
    print("alpha")
    #print(prefix  + " alpha rep1 rep2 ", reg.score(np.asarray(merged['alpha_x']).reshape(-1, 1), np.asarray(merged['alpha_y'])))
    print(prefix + " alpha rep1 rep2 ", stats.pearsonr(merged['alpha_x'], merged['alpha_y']))
    plt.close()

    plt.figure(8)
    plt.scatter(merged['mu_x'], merged['alpha_x'])
    plt.xlabel("rep1 mu")
    plt.ylabel("rep1 alpha")
    model = LinearRegression(fit_intercept = False)
    reg = model.fit(np.asarray(merged['mu_x']).reshape(-1, 1), np.asarray(merged['alpha_x']))
    #plt.title(prefix + " mu vs alpha rep1 " + str(reg.score(np.asarray(merged['mu_x']).reshape(-1, 1), np.asarray(merged['alpha_x']))))
    plt.title(prefix + " mu vs alpha rep1" + str(stats.pearsonr(merged['mu_x'], merged['alpha_x'])))
    plt.savefig(prefix + "_rep1_mualpha.png")
    print("mu vs alpha x")
    #print(prefix + " mu vs alpha rep1 ", reg.score(np.asarray(merged['mu_x']).reshape(-1, 1), np.asarray(merged['alpha_x'])))
    print(prefix + " mu vs alpha rep1", stats.pearsonr(merged['mu_x'], merged['alpha_x']))
    plt.close()

    plt.figure(9)
    plt.scatter(merged['mu_y'], merged['alpha_y'])
    plt.xlabel("rep2 mu")
    plt.ylabel("rep2 alpha")
    model = LinearRegression(fit_intercept = False)
    reg = model.fit(np.asarray(merged['mu_y']).reshape(-1, 1), np.asarray(merged['alpha_y']))
    #plt.title(prefix + " mu vs alpha rep2 " + str(reg.score(np.asarray(merged['mu_y']).reshape(-1, 1), np.asarray(merged['alpha_y']))))
    plt.title(prefix + " mu vs alpha rep2" + str(stats.pearsonr(merged['mu_y'], merged['alpha_y'])))
    plt.savefig(prefix + "_rep2_mualpha.png")
    print("mu vs alpha y")
    #print(prefix + " mu vs alpha rep2 ", reg.score(np.asarray(merged['mu_y']).reshape(-1, 1), np.asarray(merged['alpha_y'])))
    print(prefix + " mu vs alpha rep2", stats.pearsonr(merged['mu_y'], merged['alpha_y']))
    plt.close()

def main():
    args = parse_arguments()
    rep1 = read_sc(args.rep1)
    rep2 = read_sc(args.rep2)
    plot_rep1_rep2(rep1, rep2, args.output)

if __name__ == "__main__":
    main()
