"""
Spectra class
------- -----
Analyze spectral data using combination of numpy, scipy, peak-o-mat and some simple algorithms
@author: Rohan Isaac
"""
from __future__ import division
from scipy import signal # for find_peaks_cwt
import numpy as np
import sys  # to add to PATH global variable
import os   # to find working dir
import math # for misc math functions

import savitzky_golay as sg # better smoothing function

# peak-o-mat stuff
PATH = os.getcwd()+'/'
sys.path.append(PATH + 'peak-o-mat-1.1.9/')

from peak_o_mat.model import Model
from peak_o_mat.spec import Spec
from peak_o_mat.fit import Fit

class Spectra:
    """ Stores spectra data in Spec (peak-o-mat) format
    
    Contains
    --------
    1. Raw data (`s` - Spec object)
    2. Cleaned data (with background removed) (`x`, `y` - numpy arrays)
    3. Background (`bg` - Spec object)
    4. Number of peaks found (`num_peaks`)
    5. Model and Model string (`m` and `model_str`)
    
    """
    def __init__(self, filename, remove_spikes=False):
        self.s = Spec(PATH + filename)
        # set model attributes
        self.max_points = len(self.s.y)
        self.guess_peak_width()
        if remove_spikes:
            self.remove_spikes()
        # better to initialize all the data points

    def find_background(self,sub_range=None,poly_deg=3,smoothing=5):
        """ Returns (bgx,bgy) background spectra of (x,y) data passed. 
        
        Procedure
        ---------
        1. Smooths data with `signal.medfilt`
        2. Splits data into subranges
        3. Finds y-value min of subranges, and assigns to mid x-value
        4. Fits n degree polynomial to min values
        5. Returns original x-data and y-values of x-data evaluated with the polynomial
        """
        if sub_range == None:
            sub_range = (self.max_points/20)
            
        # smooth y-data
        #y_data = signal.medfilt(self.s.y,smoothing)
        
        # use same poly order as bg fitting. 
        y_data = sg.savitzky_golay(self.s.y,smoothing,3)
        #y_data = signal.savitzky_golay(self.s.y,smoothing,3)
    
        # find # of sub-ranges/intervals
        intervals = int(math.ceil(len(self.s.y)/sub_range))
        
        # find min of each subinterval and place it at mid point of 
        bg_x = np.zeros(intervals)
        bg_y = np.zeros(intervals)
        for i in range(intervals):
            pos = i*sub_range
            half_range = int(math.floor(sub_range/2))
            bg_x[i] = self.s.x[pos+half_range-1]
            bg_y[i] = min(y_data[pos:(pos+sub_range)])
            
        self.slmin = signal.argrelmin(sg.savitzky_golay(y_data,9,5))
        
        # also find the 10% of min data points in the data set and poly fit them
        
        
        # find the space between peaks and fit these
            
        #set bg_x,bg_y
        bg_poly = np.polyfit(bg_x,bg_y,poly_deg)
        bg = np.poly1d(bg_poly)
        self.bg = Spec(self.s.x,bg(self.s.x),'background')
        
    def subtract_background(self):
        """ Subtract background `bg` from raw spectra `s` """
        self.s = self.s - self.bg

    def find_peaks(self,lower=None,upper=None,threshold=5,limit=20):
        """ Find peaks in data set `s` using continous wavelet transformation (`scipy.signal`) on the y data passed. Updates `peak_pos` to list of indicies associated with peaks. Also counts peaks in `peak_pos` 
        
        Options
        -------
        lower: lower limit of peak size to search for
        upper: upper limit for same (in index values)
        threshold: min % of peak value has to be to count as a peak
        limit: max limit of peaks to detect (sorted by intensity)        
        """
        if lower == None:
            lower = self.test_peak_width*1
        if upper == None:
            upper = self.test_peak_width*5
        print "Lower: " , lower
        self.peak_pos = signal.find_peaks_cwt(self.s.y,np.arange(lower,upper))
        
        # remove peaks that aren't very big
        self.peak_pos = [i for i in self.peak_pos if (self.s.y[i]/self.data_max) > (threshold/100)]  
        
        # remove peaks that are pretty close together?
        
        
        
        # only use the top limit peaks
        # first sort the peak list by size of y
        print np.argsort(self.s.y[self.peak_pos])
        self.peak_pos = self.peak_pos[1:limit]
        self.num_peaks = len(self.peak_pos)

    def build_model(self, max_width=25, peak_type='LO'):
        """ Builds a peak-o-mat model `m` of peaks of default type lorentzian (LO). Uses some basic algorithms to determine initial parameters for amplitude and fwhm (limit on fwhm to avoid fitting background as peaks. """   
        background = 'CB'
        self.model_str = background
        for i in range(self.num_peaks):
            # build model string of peaks
            self.model_str = self.model_str + ' ' + peak_type + str(i+1)
        
        self.m = Model(self.model_str)
        
        params = [{'const':0.0}]
        for i in self.peak_pos:
            # Find basic characteristics (position,height,fwhm) 
            amplitude = self.s.y[i]
            position = self.s.x[i]
            # Find fwhm for each peak
            half_max = amplitude/2
            left = right = i
            # make sure index does not get out of bounds
            while (self.s.y[left]>half_max and left > 0 ):
                left = left-1
            while (self.s.y[right] > half_max and right < (self.max_points-1)):
                right = right + 1
            # left = find index to left when height is below half_max
            # right same as above
            # find distance between these two point
            fwhm = self.s.x[right] - self.s.x[left]
            
            # make sure doesn't blow up
            if fwhm > max_width:
                fwhm = max_width
            
            params.append({'amp': amplitude, 'fwhm': fwhm, 'pos': position })
        
        self.m.set_parameters(params)
        
        # Evaluate the model at the 5 times the precision of the data
        # ???? Is this helpful at all ????? check both ways
        self.x = np.linspace(min(self.s.x),max(self.s.x),5*self.max_points)
        self.y = self.m.evaluate(self.x)
        
    def fit_data(self):
        """ Attempt to fit data using peak-o-mat Fit function generated model. Updated model with fit parameters. """
        self.f = Fit(self.s, self.m)
        self.res = self.f.run()
        self.m.update_from_fit(self.res)
        self.y = self.m.evaluate(self.x)
        
    def output_results(self):
        """ Output fit paramters as csv values with errors"""
        for i in ["pos","fwhm", "area", "amp"]:
            print i
            print self.m.parameters_as_csv(selection=i,witherrors=True)
    
    # ---
    # Helper functions, might make private
    # ---     
    
    def guess_peak_width(self,max_width = 50):
        """ Find an initial guess for the peak with of the data imported, use in peak finding and model buildings
        
        Details
        -------
        Locates the max value in the data
        Finds the peak width associated with this data (using simple technique in model building)
        """
        self.data_max = np.amax(self.s.y)
        self.data_max_pos = np.amax( np.where(self.data_max==self.s.y) )
        self.test_peak_width = self.find_fwhm(self.data_max_pos)
        
        
    def remove_spikes(self,strength = 0.5):
        """ Attempts to remove spikes in data set using a simple test of the pixels around it. Fractional value of strenght needed. Does not modify data, but returns a copy. """
        mean = lambda x,y: (x+y)/2
        y = self.s.y
        data_max = np.amax(y)
        for i in range(1,len(y)-1):
            if (np.abs( y[i] -  mean(y[i-1],y[i+1]) )/ data_max ) > strength:
                y[i] = mean(y[i-1],y[i+1])
        return y
        
    def find_fwhm(self,position, max_width=50):
        """ Find the fwhm of a point using a very simplisitic algorigthm. Works on base data set. Has a hard limit for peak width """
        left = right = position
        half_max = self.s.y[position]/2
        
        # change max_widht to function of data set
        
        # make sure index does not get out of bounds
        while (self.s.y[left] > half_max and left > 0 ):
            left = left-1
        while (self.s.y[right] > half_max and right < (self.max_points-1)):
            right = right + 1
            
        # left = find index to left when height is below half_max
        # right same as above
        # find distance between these two point
        fwhm = self.s.x[right] - self.s.x[left]
            
        # make sure doesn't blow up
        if fwhm > max_width:
            fwhm = max_width
            
        return fwhm
