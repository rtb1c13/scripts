#!/usr/bin/env python

# Fitting & plotting of lineshapes.
# Averages across integrals data output by IR_lineshapes.py
# 
# Author: Richard Bradshaw, R.T.Bradshaw@soton.ac.uk

# Requirements: lmfit, argparse, numpy, matplotlib, lmcurvefit module, smooth module

from __future__ import print_function

from lmcurvefit import *
from smooth import *
import numpy as np
import scipy.stats as stats
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import sys, argparse

# Argparser

def parse():
   parser = argparse.ArgumentParser()
   parser.add_argument("-i","--inputs",help="Identifiers for input integral and frequency files. Integral files will be loaded as 'integrals_[name].txt, freqs as [name]_wavenums.txt",nargs='+',type=str,required=True)
   parser.add_argument("-m","--models",help="Models for each peak to be fitted. Defaults to a single Gaussian model",nargs='+',type=str,default=['gaussian'])
   parser.add_argument("-p","--prefixes",help="Prefix used for each peak/model",nargs='+',type=str)
   parser.add_argument("-c","--centres",help="Initial guess for centrepoints of each peak/model. Defaults to 1700 cm-1.",nargs='+',type=float,default=[1700.0])
   parser.add_argument("-s","--sigmas",help="Initial guess for sigmas of each peak/model. Defaults to 5.0 cm-1",nargs='+',type=float,default=[5.0])
   parser.add_argument("-o","--output",help="Filename of fitted model image to output. Defaults to 'fitted_model.png'",type=str)
   parser.add_argument("-sm","--smoothing",help="Flag to perform cubic spline smothing of input data prior to fitting. Defaults to false.",dest='smflag', action='store_true')
   parser.set_defaults(smflag=False)

   if len(sys.argv)==1:
      parser.print_help()
      sys.exit(1)

   args = parser.parse_args()
   return args


# Read in integrals & fields 
def read_files(files):
    intfiles = [ 'integrals_' + f + '.txt' for f in files ]
    freqfiles = [ f + '_wavenums.txt' for f in files ]
    def read_integrals(currfile):
        integrals = 0
        reals, imags = np.loadtxt(currfile, unpack=True)
        integrals += reals + 1j*imags
        return integrals
    def read_freqs(currfile):
        return np.loadtxt(currfile)
    aveints = np.sum(map(read_integrals,intfiles),axis=0) / len(files)
    avefreqs = np.sum(map(read_freqs,freqfiles),axis=0) / len(files)
    return aveints, avefreqs

# Do FFT
def calc_fft(integrals, offset=0, dt=2e-15, t_len=2**12):
    # offset = offset of frequencies in wavenumbers, dt = timestep, t_len = length of fft sample
    # FFT parameters
    Fs = 1./dt # Sampling rate s-1
    fftn = t_len*(2**4) # Samples in FFT signal inc. padding
    idxs = np.arange(fftn) # Sample indices
    T = fftn/Fs # Time length of fft signal in s
    frq = idxs/T # Set of frequencies
    frq = frq/2.99792458e10 # Frq -> Wavenumbers

    # Do fft with 16*padding
    fftout = np.fft.fft(integrals,n=fftn)
    negfft = fftout[fftn//2:] # First half are at -ve frequencies
    posfft = fftout[:fftn//2]
    negfrq = -1*frq[1:fftn//2+1] # Same as arange(-len(fftout)/2,len(fftout)/2)
    posfrq = frq[:fftn//2]
    # Only return real portion for plotting
    return np.concatenate((negfrq[::-1],posfrq)) + offset, np.concatenate((negfft,posfft)).real


# Smooth data

def smooth_data(data, spline_kwargs=None, plot_kwargs=None):
    """Smooths y-values of data with a 1D smoothed spline fit.
       Usage: smooth_data(data, [spline_kwargs, plot_kwargs])
       Returns: smoothed dataset and matplotlib figure of smoothing"""
    
    k = Smoothed(data)
    if (spline_kwargs):
        k.do_smooth(**spline_kwargs)
    else:
        k.do_smooth()
    if (plot_kwargs):
        outfig = k.plot(**plot_kwargs)
    else:
        outfig = k.plot()
    data[1] = k.smoothed_ys
    return data, outfig

# Do fit

def single_fit(data, run_fit=True, fittype='gaussian', prefix='', **kwargs):
    outfit = Fit(data)
    if run_fit:
        fitdict = {'linear': outfit.linear,
                   'gaussian': outfit.gaussian,
                   'lorentzian': outfit.lorentzian,
                   'voigt': outfit.voigt}
        fitdict[fittype](prefix=prefix,**kwargs)
    return outfit

def combine_fits(data,models,params=None):
    # Models is iterable, params is dict of prefix_center etc
    multi = Multifit(*models)
    if params:
        for model in models:
            multi.init_params(prefix=model._mod.prefix,center=params[model._mod.prefix+'center'], sigma=params[model._mod.prefix+'sigma'],amplitude=1500.0)
    else:
        for model in models:
            multi.init_params(prefix=model._mod.prefix)
    multi.make_mod()
    multi.do_multifit(data)
    return multi

# Save plot
def plot_fit(model, title=None, plot_init=False, annotate=True, fn="fitted_model.png", offset=0):
    data_keys = {'color':'black',
                 'linestyle':'dashed',
                 'marker':'',
                 'linewidth':1}
    fit_keys = {'color':'blue',
                'linestyle':'solid',
                'marker':'',
                'linewidth':2}
    fig = plt.figure()
    fig.ax = model.plots(xlabel='Frequency / cm$^{-1}$', ylabel='Intensity (arbitrary units)',
                         data_kws=data_keys, fit_kws=fit_keys)
    # Don't plot initial models by default
    if not plot_init:
        fig.ax.lines.pop(0) # Remove 'init' plot and redo legend
        fig.ax.relim() # Rescale limits of axes without 'init' plot
        fig.ax.autoscale_view()
        fig.ax.legend(['Final fit','Data'])
    else:
        fig.ax.legend(['Initial fit','Final fit','Data'])
    # Set annotation with centre/FWHM values
    if annotate:        
        gs = gridspec.GridSpec(4,1)
        fig.ax.set_position(gs[0:3].get_position(fig))
        fig.ax.set_subplotspec(gs[0:3])
        fig.tight_layout()
        tableax = fig.add_subplot(gs[3])
        tableax.axis('tight')
        tableax.axis('off')
        col_labels=['Centre','FWHM']
        row_labels=[ "Peak %s" % idx for idx,val in enumerate(model._mods,1) ]
        table_vals=[ [model.total_fit.values[idx._mod.prefix+'center'], model.total_fit.values[idx._mod.prefix+'fwhm']] for idx in model._mods ] # Centre/FWHM vals as list of lists by row
        table_vals= [ [ np.round(i,2) for i in inner ] for inner in table_vals ] # Rounding
        the_table = tableax.table(cellText=table_vals,
                                  colWidths = [0.1]*3,
                                  rowLabels=row_labels,
                                  colLabels=col_labels,
                                  loc='center')
    # Set title
    if title:
        fig.ax.set_title(title)
    else:
        fig.ax.set_title("Multiple model fit, $\chi^2$ = %8.5f" % model.total_fit.chisqr)
    fig.ax.set_xlim((offset-100,offset+100))
    fig.savefig(fn,dpi=300)
    
# Set paramdict
def set_params(cents, sigs, prefs=''):
    prefcents, prefsigs = [ i + 'center' for i in prefs ], [ i + 'sigma' for i in prefs ]
    paramdict = {}
    paramdict.update(zip(prefcents,cents))
    paramdict.update(zip(prefsigs,sigs))
    return paramdict


# Write log (fit_report)
def write_log(model,fn="logfile.txt"):
    with open(fn,'w') as f:
        print(model,file=f)

### Main below here ###
if __name__ == "__main__":
    global inp_args
    inp_args = parse()
    integrals, freqs = read_files(inp_args.inputs)
    xs, ys = calc_fft(integrals, offset=np.mean(freqs))
    tot_data = np.vstack((xs,ys))
    np.savetxt("combined_lineshape.txt", zip(xs,ys))
    if inp_args.smflag:
        tot_data, smoothfig = smooth_data(tot_data)
        smoothfig.savefig("Smoothed_average_lineshape.png")
    # Below here needs to be updated with switches based on argparse
    if inp_args.prefixes is not None:
        fits = [ single_fit(tot_data, prefix=z[0], fittype=z[1]) for z in zip(inp_args.prefixes,inp_args.models) ]
        paramdict = set_params(inp_args.centres, inp_args.sigmas, inp_args.prefixes)
    else:
        fits = [ single_fit(tot_data, fittype=m) for m in inp_args.models ]
        paramdict = set_params(inp_args.centres, inp_args.sigmas)
    mm = combine_fits(tot_data, fits, params=paramdict)
    write_log(mm)
    plot_fit(mm,offset=np.mean(freqs))
# Add function to write out values of composite plot, and/or y-max and FWHM, here:

#xs = np.arange(-100.,100.,1.)
#ys = stats.norm.pdf(xs,-20.,10.,)
#noise = ys * np.random.normal(1.0,0.05,len(ys))
#data = np.vstack((xs,noise))
#outfit = Fit(data)
#outfit.gaussian(prefix='gau1_')
#fig = plt.figure()
#data_keys = {'color':'black',
#            'linestyle':'dashed',
#            'marker':''}
#fit_keys = {'color':'blue',
#            'linestyle':'solid',
#            'marker':'',
#            'linewidth':'2'}
#fig.ax = outfit.plots(xlabel='Frequency',ylabel='Intensity',data_kws=data_keys,fit_kws=fit_keys)
#fig.ax.lines.pop(0) # Remove 'init' plot and redo legend
#fig.ax.relim() # Rescale limits of axes without 'init' plot
#fig.ax.autoscale_view()
#fig.ax.legend(['fit_test','data_test'])
#fig.ax.set_title("Center = %6.2f, FWHM = %6.2f, Chi2 = %8.5f" % (outfit.fit.values['gau1_center'], outfit.fit.values['gau1_fwhm'], outfit.fit.chisqr))
#fig.suptitle("Test of single fit")
#fig.savefig("single_fit.png")
#
# Double fit below here

#ys2 = stats.norm.pdf(xs,20.,15.,)
#noise2 = ys2 * np.random.normal(1.0,0.1,len(ys))
#data2 = np.vstack((xs,noise2))
#outfit2 = Fit(data2)
#outfit2.gaussian(prefix='gau2_')
##
#compnoise = noise + noise2
#compdata = np.vstack((xs,compnoise))
#gaufit1 = single_fit(compdata,prefix='gau1_')
#gaufit2 = single_fit(compdata,prefix='gau2_')
#mm = combine_fits(compdata,(gaufit1,gaufit2),params=paramdict)
#write_log(mm)
#
#multi = Multifit(outfit,outfit2)
#multi.init_params(prefix='gau1_',center=-20., sigma=20.)
#multi.init_params(prefix='gau2_',center=40., sigma=10.)
#multi.make_mod()
#multi.do_multifit(compdata)
#fig = plt.figure()
#fig.ax = multi.plots(xlabel='Frequency',ylabel='Intensity',data_kws=data_keys,fit_kws=fit_keys)
#fig.ax.lines.pop(0) # Remove 'init' plot and redo legend
#fig.ax.relim() # Rescale limits of axes without 'init' plot
#fig.ax.autoscale_view()
#fig.ax.legend(['fit_test','data_test'])
#fig.ax.set_title("""Center1 = %6.2f, FWHM = %6.2f\nCenter2 = %6.2f, FWHM = %6.2f, Chi2 = %8.5f"""
#                 % (multi.total_fit.values['gau1_center'], 
#                    multi.total_fit.values['gau1_fwhm'], 
#                    multi.total_fit.values['gau2_center'],
#                    multi.total_fit.values['gau2_fwhm'],
#                    multi.total_fit.chisqr))
#fig.suptitle("Test of multiple fit")
#fig.savefig("multi_fit.png")
