"""
Attempt a functionalized verison of spectral analysis
"""
from __future__ import division
from helper_functions import savitzky_golay, percent_diff
from scipy import interpolate, signal
import numpy as np

def find_background(x_data, y_data, percent=3, poly_deg=3, smoothing=5):
    """ Finds the background of the the spectrum passed, after smoothing data.
    
    Parameters
    ----------
    y_data: array
        The intensity of spectrum to be analyzed. Equal x-spacing assumed
    percent: int (default = 3)
        The percent that data points in the background can be different from the
        previous one
    poly_deg: integer
        Positve integer greater than 1. High polynomials will results in bad
    
    Returns
    -------
    Array of background spectrum, which can then be subtracted from spectrum. 
    Same size as y_data
    
    """
    print "Finding background ... " 

    # smooth the data
    y_smooth = savitzky_golay(y_data,smoothing,poly_deg)
    
    # start with the first x,y postion
    bg_x = [x_data[0]]
    bg_y = [y_smooth[0]]
    
    # probably could write more efficiently
    
    # if the next point is less percent different from the previous, add it to 
    # the list
    for i in range(1,len(y_data)):
        if percent_diff(y_smooth[i], y_smooth[i-1]) < percent:
            bg_x.append(x_data[i])
            bg_y.append(y_smooth[i])      
    
    # intepolate this 
    bg_spline = interpolate.interp1d(bg_x, bg_y, kind='quadratic')
    
    # return the backgroud y values
    return bg_spline(x_data)

def find_peaks(y_data, lower=1, upper=5, threshold=5, limit=20):
    """ Find peaks in actve data set using continous wavelet 
    transformation from `scipy.signal`. Smooth before search and filter peaks 
    after.
    
    Parameters
    ---------
    y_data: np.array
        Array of y-data values associated with spectra
    lower: int (optional: default 1)
        lower bound of peak size in terms of data points
    upper: int (optional: default 5)
        upper bound of data points spanned by peak
    threshold: int (default: 5)
        min % of peak value has to be to count as a peak
    limit: int (default: 20)
        max limit of peaks to detect (sorted by intensity)  
    
    Returns
    -------
    peak_pos: list
        indicies associated with peak positions
       
    """
    print "Looking for peaks ... "  

    peak_pos = signal.find_peaks_cwt(y_data,np.arange(lower,upper),min_snr=2)
    print "Found ", len(peak_pos), " peaks."
    
    # remove peaks that are not above the threshold.
    peak_pos = [i for i in peak_pos if (y_data/max(y_data)) > (threshold/100)]  
    print "After filtering out peaks below ", threshold, "percent, we have ", \
        len(peak_pos), " peaks."
    
    # zip the intensities and the peak positions together, and sort by
    # intensities in descending order. 
    peak_pos = [y for (x,y) in sorted(zip(y_data[peak_pos],peak_pos),reverse=True)]
    # truncate the list to the required size
    peak_pos = peak_pos[0:limit]
    print "Using ", limit, " peaks at ", peak_pos
    
    return peak_pos