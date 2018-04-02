#!/usr/bin/env python

# Matplotlib histogram of discrete variable
# overlaid with PDF (i.e. PMF) a Poisson distribution
# centered at the same mean

import numpy as np
import scipy.stats as stats
import matplotlib.pyplot as plt
import argparse
from cycler import cycler # For Matplotlib colour cycle

### Define defaults for matplotlib plots
plt.rc('lines', linewidth=1.5, markersize=4)
plt.rc('axes', prop_cycle=(cycler('color', ['k','b','r','orange','c','m','y','g'])), # Color cycle defaults to black
       labelweight='heavy', labelsize=22, titlesize=22) # Default fontsizes for printing
plt.rc('axes.spines', top=False, right=False) # Switch off top/right axes
plt.rc('legend', fontsize=12) # Default fontsizes for printing
plt.rc('xtick', labelsize=16) # Default fontsizes for printing
plt.rc('ytick', labelsize=16) # Default fontsizes for printing
plt.rc('figure', figsize=(11,8.5), titlesize=22, titleweight='heavy') # Default fontsizes for printing

parser = argparse.ArgumentParser()
parser.add_argument("-f","--files",help="Files of data. Will be read in as integers (i.e. discrete variable)",nargs='+',type=str,required=True)
parser.add_argument("-n","--name",help="Name of data. Will be used in row names and appended to filenames", type=str)

args = parser.parse_args()



def combine_files(files):
    """Read in a list of files and return concatenated
       Numpy array."""
    cat = np.array([], dtype=np.int32)
    for fn in files:
        x = np.loadtxt(fn, dtype=np.int32)
        cat = np.concatenate((cat,x))
    return cat

def histogram(values):
    """Histogram an array of values. Bin edges are calculated
       to be +/- 5 from the max and min of the given array
       Usage: histogram(values)
       Returns: Array of normalised density in each bin,
                Array of bin edges,
                Array of matplotlib patch objects"""
    minval, maxval = np.min(values), np.max(values)
    dens, bins, patches = plt.hist(values,\
                                   bins=np.arange(minval-5, maxval+6, 1),\
                                   density=True)
    return dens, bins[:-1], patches # Last bin is closed-ended 

def poisson_pmf(values, xs):
    """Generate an idealised probability mass function
       for a Poisson distribution centered at a given value.
       If more than one value is provided, the mean will
       be used as the centre of the Poisson distribution.
       The PMF is evaluated at the x-points provided.

       Usage: poisson_pmf(values, xs)
       Returns: Array of probability densities at points in xs"""

    xbar = np.mean(values)
    p_obj = stats.poisson(xbar)
    return p_obj.pmf(xs)

def norm_pdf(values, xs):
    """Generate an idealised probability density function
       for a normal distribution centered at the mean of
       given values, and with width std. dev. of the given
       values. The PDF is evaluated at the x-points provided.

       Usage: norm_pdf(values, xs)
       Returns: Array of probability densities at points in xs"""

    xbar = np.mean(values)
    sdev = np.std(values)
    n_obj = stats.norm(loc=xbar, scale=sdev)
    return n_obj.pdf(xs)

def desc_stats(values, name):
    """Writes out some descriptive statistics of the values
       read in. Prepends output with 'name'"""
    xbar = np.mean(values)
    var = np.var(values)
    sigma = np.std(values)
    error = stats.sem(values)
    
    with open("stats_%s.txt" % name,'a') as f:
        f.write("#Description    Mean   Variance   Std. Dev.  Std. Err.\n")
        f.write("%s %7.3f  %7.3f  %7.3f  %7.3f\n" \
                % (name.ljust(len(name)+8), xbar, var, sigma, error))

def overlay_hist_poisson(xs, hist_ys, poisson_ys, name, values):
    """Plot an overlay of histogrammed data and an idealised
       Poisson distribution of the same mean. Skew and kurtosis
       of the data are also printed to the descriptive statistics file."""

    with open("stats_%s.txt" % name,'a') as f:
        f.write("#Description    Skew (+ve = left skewed)   Kurtosis (+ve = longer-tailed)\n")
        f.write("%s %7.3f                 %7.3f\n" \
                % (name.ljust(len(name)+8), stats.skew(values), stats.kurtosis(values)))
    fig = plt.figure()
    ax = fig.gca()
    ax.set_title("Data overlaid with Poisson distribution")
#    ax.set_xticks(xs)
    ax.set_ylabel("Probability density")
    ax.set_xlabel("Counts")
    ax.plot(xs, hist_ys, label="Histogrammed data, mean = %5.2f" % np.mean(values) )
    ax.plot(xs, poisson_ys, label="Poisson distribution", linestyle='-')
    plt.legend(loc="upper left")
    plt.tight_layout()
    plt.savefig("Histogram+Poisson_overlay_%s.png" % name, dpi=300)

def overlay_hist_normal(xs, hist_ys, norm_ys, name, values):
    """Plot an overlay of histogrammed data and an idealised
       normal distribution of the same mean & std. dev."""

    fig = plt.figure()
    ax = fig.gca()
    ax.set_title("Data overlaid with continuous normal distribution")
#    ax.set_xticks(xs)
    ax.set_ylabel("Probability density")
    ax.set_xlabel("Counts")
    ax.plot(xs, hist_ys, label="Histogrammed data, mean = %5.2f" % np.mean(values) )
    ax.plot(xs, norm_ys, label="Normal distribution", linestyle='-')
    plt.legend(loc="upper left")
    plt.tight_layout()
    plt.savefig("Histogram+normal_overlay_%s.png" % name, dpi=300)


### Main below here ###
if __name__ == "__main__":
    data = combine_files(args.files)
    try:
        desc_stats(data, args.name)
        hist_ys, xs, _ = histogram(data)
        ys1 = poisson_pmf(data, xs)
        overlay_hist_poisson(xs, hist_ys, ys1, args.name, data)
        ys2 = norm_pdf(data, xs)
        overlay_hist_normal(xs, hist_ys, ys2, args.name, data)
    except AttributeError:
        desc_stats(data, "File_1")
        hist_ys, xs, _ = histogram(data)
        ys1 = poisson_pmf(data, xs)
        overlay_hist_poisson(xs, hist_ys, ys1, "File_1", data)
        ys2 = norm_pdf(data, xs)
        overlay_hist_normal(xs, hist_ys, ys2, "File_1", data)


