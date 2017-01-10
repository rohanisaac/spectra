from __future__ import division
import numpy as np
from scipy import signal
from numpy import log, exp

def find_peaks(x, y, width=5, threshold=5, limit=20):
    """ Find peaks in active data set using continuous wavelet transformation

    Parameters
    ----------
    width: float (default=5)
        estimate of peak size in x data units
    threshold: float (default=5)
        min percent of max to count as a peak (eg 5 = only peaks above 5
        percent reported)
    limit: int
        max limit of peaks to report (sorted by intensity)

    Returns
    -------
    peak_pos : list
        indices associated with peak positions

    Notes
    -----
    Do I really want to include x data here? Should I report x-positions in
    terms of indices or in terms of x positions
    """
    # scale factor to remove units from x data
    xscale = len(x) / (max(x) - min(x))

    lower = width * xscale * 0.75
    upper = width * xscale * 1.25

    peak_pos = signal.find_peaks_cwt(y, np.arange(lower, upper))

    data_max = np.max(y)
    # remove peaks that are not above the threshold.
    peak_pos = [i for i in peak_pos if
                (y[i] / data_max) > (threshold / 100)]


    # only use the most intense peaks, zip two lists together,
    # make the y-values as the first item, and sort by it (descending)
    peak_pos = [yval for (index, yval) in sorted(zip(y[peak_pos], peak_pos),
                                                 reverse=True)]

    peak_pos = sorted(peak_pos[0:limit])

    return peak_pos


def gaussian(x, x0, amp, fwhm):
    """
    Gaussian peak with
    amp: Amplitude (height)
    x0: Peak center
    fwhm: Full width at half maximum
    """
    return amp*(exp(-4*log(2) * ((x-x0)/fwhm)**2))


def lorentzian(x, x0, amp, fwhm):
    """
    Lorentzian peak with:
    amp: Amplitude (height)
    x0: Peak center
    fwhm: Full width at half maximum
    """
    return amp/(1+(2*(x-x0)/fwhm)**2)


def voigt(x, x0, amp, fwhm, alpha=0.5):
    """
    Normalized pseudo voigt peak with:
    amp: Amplitude (height)
    x0: Peak center
    fwhm: Full width at half maximum
    alpha: fraction lorentzian (1-fraction gaussian)
    """
    return ((1 - alpha)*gaussian(x, x0, amp, fwhm) +
            alpha*lorentzian(x, x0, amp, fwhm))

def guess_peak_width(x, y, max_width=None):
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
        max_width = len(y) / 5

    data_max_pos = np.argmax(y)
    test_peak_width = find_fwhm(x, y, data_max_pos)

    return test_peak_width

def find_fwhm(x, y, position):
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
    half_max = y[position] / 2

    # change max_width to function of data set

    # make sure index does not get out of bounds
    while (y[left] > half_max and left > 0):
        left = left - 1
    while (y[right] > half_max and right < (len(y) - 1)):
        right = right + 1

    # left = find index to left when height is below half_max
    # right same as above
    # find distance between these two point
    fwhm = x[right] - x[left]

    return fwhm


