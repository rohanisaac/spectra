"""
Attempt a functionalized verison of spectral analysis
"""
from __future__ import division
from helper_functions import savitzky_golay
from scipy import interpolate

def find_background(x_data, y_data, poly_deg=3, smoothing=5):
    """ Finds the background of the the spectrum passed, after smoothing data.
    
    Paramters
    ---------
    y_data: array
        The intensity of spectrum to be analyzed. Equal x-spacing assumed
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

    bg_x = [x_data[0]]
    bg_y = [y_smooth[0]]
    
    for i in range(1,len(y_data)):
        if percent_diff(y_data[i], y_data[i-1]) < 3:
            bg_x.append(x_data[i])
            bg_y.append(y_smooth[i])      
    
    bg_spline = interpolate.interp1d(bg_x, bg_y, kind='quadratic')
    
    return bg_spline(x_data)

def find_peaks(self, lower=None, upper=None, threshold=5, limit=20):
    """ Find peaks in actve data set using continous wavelet 
    transformation from `scipy.signal`
    
    Variables modified
    ------------------
    peak_pos: list of indicies associated with peak positions
    num_peaks: # of peaks found
    
    Options
    -------
    lower: lower limit of peak size to search for
    upper: upper limit for same (in index values)
    threshold: min % of peak value has to be to count as a peak
    limit: max limit of peaks to detect (sorted by intensity)        
    """
    print "Looking for peaks ... "  
    
    if lower == None:
        lower = self.test_peak_width*1
    if upper == None:
        upper = self.test_peak_width*5

    peak_pos = signal.find_peaks_cwt(self.active,np.arange(lower,upper),min_snr=2)
    
    print "Found ", len(peak_pos), " peaks."
    
    
    # remove peaks that are not above the threshold.
    peak_pos = [i for i in peak_pos if (self.active[i]/self.data_max) > (threshold/100)]  
    
    # remove peaks that are pretty close together?
    # only use the top limit peaks
    # first sort the peak list by size of y
    # print np.argsort(self.s.y[self.peak_pos])
    
    print "After filtering out peaks below ", threshold, "percent, we have ", len(peak_pos), " peaks."
    
    #only use the top limit peaks, zip two list together, 
    # make the y-values as the first item, and sort by it (in descending)
    peak_pos = [y for (x,y) in sorted(zip(self.active[peak_pos],peak_pos),reverse=True)]
    self.peak_pos = peak_pos[0:limit]
    self.num_peaks = len(self.peak_pos)
    print "Using ", self.num_peaks, " peaks at ", self.peak_pos
    
def percent_diff(a, b):
    return 200*(a-b)/(a+b)