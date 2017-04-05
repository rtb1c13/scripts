#!/usr/bin/env python

# Author: Richard Bradshaw, R.T.Bradshaw@soton.ac.uk
# Module to fit various curves to provided x/y data
# Current available curves: Linear, Gaussian, Lorentzian, Voigt 

# Requirements: lmfit, numpy, matplotlib (as dependencies of lmfit)

from lmfit.models import LinearModel,GaussianModel,LorentzianModel,VoigtModel

class FitError(Exception):
    """Exception in lmfit wrapper"""

class Fit():
    """Class to contain methods for curve fitting from
       lmfit package."""
    def __init__(self, data):
        """Usage: Fit(data). Initialises data as xs and ys attributes.
           Data should be a 2d numpy array of x and y values."""
        if len(data) != 2:
            raise FitError("""Your data is formatted incorrectly -  
                              it should be a 2D array of all x-,
                              then all y-values""")
        self.xs = data[0]
        self.ys = data[1]

    def __str__(self):
        """Prints lmfit fit report for the current object"""
        try:
            return self.fit.fit_report()
        except AttributeError:
            return "No fit yet performed for this object."

    def linear(self, **kwargs):
        """Linear fit of data. Usage: self.linear([**kwargs])
           kwargs are passed to a lmfit LinearModel."""
        self._mod = LinearModel(**kwargs)
        pars = self._mod.guess(self.ys, self.xs)
        self.fit = self._mod.fit(self.ys, pars, x=self.xs)

    def gaussian(self, **kwargs):
        """Gaussian fit of data. Usage: self.gaussian([**kwargs])
           kwargs are passed to a lmfit GaussianModel."""
        self._mod = GaussianModel(**kwargs)
        pars = self._mod.guess(self.ys, self.xs)
        self.fit = self._mod.fit(self.ys, pars, x=self.xs)

    def lorentzian(self, **kwargs):
        """Lorentzian fit of datia. Usage: self.lorentzian([**kwargs])
           kwargs are passed to a lmfit LorentzianModel."""
        self._mod = LorentzianModel(**kwargs)
        pars = self._mod.guess(self.ys, self.xs)
        self.fit = self._mod.fit(self.ys, pars, x=self.xs)

    def voigt(self, **kwargs):
        """Voigt fit of data. Usage: self.voigt([**kwargs])
           kwargs are passed to a lmfit VoigtModel."""
        self._mod = VoigtModel(**kwargs)
        pars = self._mod.guess(self.ys, self.xs)
        self.fit = self._mod.fit(self.ys, pars, x=self.xs)

    def plots(self, **kwargs):
        """Returns matplotlib axes with original data, fitted
           function & initial model.
           Usage: self.plots([**kwargs])
           kwargs are passed to lmfit.ModelResult.plot_fit"""
        try:
            return self.fit.plot_fit(**kwargs)
        except AttributeError:
            raise FitError("No fit yet performed to plot!")

class Multifit():
    """Composite model from two or more underlying models,
       passed as Fit objects defined in lmcurvefit. Models in
       Fit objects should have been defined with unique prefixes
       otherwise output in the composite model may be confusing/incorrect."""
    def __init__(self, *args):
        """Usage: Multifit(model1, [model2, model3, ...])
           Reads in individual models to perform a composite fit.
           Models should be read in as Fit objects with their
           own defined individual models already assigned"""
        self._mods = args
        try:
           self._pars = self._mods[0]._mod.make_params()
        except AttributeError:
           raise FitError("""Your starting models should be read in as Fit objects
                             each with a single model already defined.""")
        for fit in self._mods[1:]:
            self._pars.update(fit._mod.make_params())

    def __str__(self):
        """Prints lmfit fit report for the current object"""
        try:
            return self.total_fit.fit_report()
        except AttributeError:
            return "No composite fit yet performed for this object."

    def init_params(self, prefix='', center=0, sigma=10, amplitude=10):
        """Usage: self.init_params([prefix='', center=0, sigma=10, amplitude=10])
           Sets initial guess parameters for the model defined with 'prefix'."""
        self._pars[prefix+'center'].set(center)
        self._pars[prefix+'sigma'].set(sigma)
        self._pars[prefix+'amplitude'].set(amplitude)

    def make_mod(self):
        """Usage: self.make_mod()
           Makes composite model from all models read in."""
        self._compmod = self._mods[0]._mod
        for fit in self._mods[1:]:
            self._compmod += fit._mod

    def do_multifit(self,data):
        """Usage: self.do_multifit(data)
           Performs fitting of data to composite model.
           Data should be a 2D numpy array of x and y values"""
        if len(data) != 2:
            raise FitError("""Your data is formatted incorrectly -  
                              it should be a 2D array of all x-,
                              then all y-values""")
        self.xs = data[0]
        self.ys = data[1]
        try:
           self.total_fit = self._compmod.fit(self.ys, self._pars, x=self.xs)
           self.composite_fits = self.total_fit.eval_components(x=self.xs)
        except AttributeError:
           raise FitError("""You don't seem to have a composite model - run
                             make_mod() first!""")

    def plots(self, **kwargs):
        """Returns matplotlib axes with original data, fitted
           function & initial model.
           Usage: self.plots([**kwargs])
           kwargs are passed to lmfit.ModelResult.plot_fit"""
        try:
            return self.total_fit.plot_fit(**kwargs)
        except AttributeError:
            raise FitError("No fit yet performed to plot!")
