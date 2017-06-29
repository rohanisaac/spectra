"""
Conversion of spectra class to functions
author: Rohan Isaac
"""


import numpy as np
from numpy import sqrt, pi
import re
import pandas as pd
from scipy import signal
from lmfit import Model
from lmfit.models import PolynomialModel
from lmfit.lineshapes import lorentzian, gaussian, voigt
from uncertainties import ufloat


# ---
# Helper functions
# ---
def print_peak_results(self):
    """
    Compute all peak results including additional outputs such as Height
    and FWHM
    """
    params = self.out.params
    print('\tPosition\tHeight\tFWHM')
    for p in range(self.num_peaks):
        center = ufloat(params['p%s_center' % p].value, params['p%s_center' % p].stderr)
        amplitude = ufloat(params['p%s_amplitude' % p].value, params['p%s_amplitude' % p].stderr)
        sigma = ufloat(params['p%s_sigma' % p].value, params['p%s_sigma' % p].stderr)
        height = self.height(amplitude, sigma)
        fwhm = self.fwhm(sigma)
        print('Peak%s\t%s\t%s\t%s' % (p, center, height, fwhm))

def guess_peak_width(self, max_width=None):
    """ Find an initial guess for the peak with of the data imported,
    use in peak finding and model buildings and other major functions,
    probably should call in the constructor

    Parameters
    ----------
    max_width : int (default = total points/5)
        Max width of peaks to search for in points

    Notes
    -------
    Locates the max value in the data
    Finds the peak width associated with this data

    Returns
    -------
    data_max : float
        max intensity of y-data
    data_max_pos : int
        index of max data
    test_peak_width :
        guess for peak width of data

    """
    if max_width is None:
        max_width = self.num_points / 5

    self.data_max = max(self.y)
    self.data_max_pos = np.argmax(self.y)
    self.test_peak_width = self.find_fwhm(self.data_max_pos)

    print("Peak width of about %s (in x-data units)" % self.test_peak_width)
    return self.data_max, self.data_max_pos, self.test_peak_width

def set_peak_width(self, width):
    """
    Sets peak width
    """
    self.test_peak_width = width

def find_fwhm(self, position):
    """ Find the fwhm of a point using a very simplistic algorithm.
    Could return very large width.

    Arguments
    ---------
    position : int
        index of peak of which we are trying to determine fwhm

    Returns
    -------
    fwhm : float
        fwhm of peak in units of x-data

    """
    left = position
    right = position
    half_max = self.y[position] / 2

    # change max_width to function of data set

    # make sure index does not get out of bounds
    while (self.y[left] > half_max and left > 0):
        left = left - 1
    while (self.y[right] > half_max and right < (self.num_points - 1)):
        right = right + 1

    # left = find index to left when height is below half_max
    # right same as above
    # find distance between these two point
    fwhm = self.x[right] - self.x[left]

    return fwhm

def filter_high_freq(self, cof=5):
    """ Filter high frequency y-data using 1-D fft

    Parameters
    ----------
    cof : int (default=5)
        Coefficient above which to truncate Fourier series.

    """
    pass

def calibrate_x(self, m, b):
    """ Applies a linear correction to the x-values """
    # Need to change the active data set
    # Save the old data etc.

    pass

def reset(self):
    """
    Reset the y data to the original value
    """
    self.y = self.y_bak

def height(self, amplitude, sigma):
    """
    Converts amplitude to height

    Factors:
    lorentzian: pi
    gaussain: sqrt(2*pi)
    voigt: sqrt(2*pi)
    """
    return amplitude / (sigma * self.afactor)

def amplitude(self, height, sigma):
    """
    Converts amplitude to height

    Factors:
    lorentzian: pi
    gaussain: sqrt(2*pi)
    voigt: sqrt(2*pi)
    """
    return height * (sigma * self.afactor)

def fwhm(self, sigma):
    """
    Converts sigma to fwhm

    Factors:
    lorentzian: 2.0
    gaussain: 2.354820
    voigt: 3.60131
    """
    return self.wfactor * sigma
