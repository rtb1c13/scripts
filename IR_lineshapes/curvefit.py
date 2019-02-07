#!/usr/bin/env python

# Author: Richard Bradshaw, R.T.Bradshaw@soton.ac.uk
# Module to fit various curves to provided data
# Current available curves: Gaussian, Lorentzian, Voigt
# Reports back various curve properties and optionally chi-squared fit value

# Requirements: Numpy, scipy.optimize v.0.17 or greater

import numpy as np
from scipy.optimize import curve_fit

class Fit():
    """Class to contain methods for curve fitting
       to various models"""
    def __init__(self,data):
        """Initialises data"""
        self.data = data
    def _linfit(self,x,a,b):
        """Linear function"""
        return (a * x) + b
    def linear(self):
        """Linear fit of data"""
        return curve_fit(self.__linfit,*self.data)
    def _gaufit(self,x,h,m,sig):
        """Gaussian function"""
        prefac = ( h / (sig * np.sqrt(2*np.pi)) )
        expo = np.exp( (-(x-m)**2) / (2*sig**2) )
        return prefac * expo
    def gaussian(self):
        """Gaussian fit of data"""
        popt, pcov = curve_fit(self._gaufit,*self.data)
        fwhm = 2*popt[2] * np.sqrt(2*np.log(2))
        return (popt,pcov,fwhm)
    def _lorfit(self,x,h,m,sig):
        """Lorentzian function"""
        return ( h / np.pi ) * ( sig / ( (x-m)**2 + sig**2) )
    def lorentzian(self):
        """Lorentzian fit of data"""
        popt, pcov = curve_fit(self._lorfit,*self.data)
        fwhm = 2*popt[2]
        return (popt,pcov,fwhm)
    def _voifit(self,x,h,m,sig,alpha):
        """Pseudo-Voigt function - linear mix of 
           Gaussian & Lorentzian models"""
        sigg = sig/np.sqrt(2*np.log(2))
        gauprefac =  ((1-alpha)*h) / (sigg * np.sqrt(2*np.pi))
        gauexpo = np.exp( (-(x-m)**2) / (2*sigg**2) )
        gau = gauprefac*gauexpo
        lor = ( (alpha*h) / np.pi ) * ( sig / ( (x-m)**2 + sig**2) )
        return gau + lor
    def voigt(self):
        """Pseudo-Voigt fit of data. Fraction of Lorentzian curve
           limited between 0 .le. alpha .le. 1"""
        minbounds = np.array([-np.inf,-np.inf,-np.inf,0.])
        maxbounds = np.array([np.inf,np.inf,np.inf,1.])
        popt, pcov = curve_fit(self._voifit,*self.data,bounds=(minbounds,maxbounds))
        fwhm = 2*popt[2]
        return (popt,pcov,fwhm)
