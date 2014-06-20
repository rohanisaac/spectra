"""
Spectra class
-------------
Analyze spectral data using combination of numpy, scipy, peak-o-mat and 
some simple algorithms
@author: Rohan Isaac
"""

# keep only one spec object, lots of y-data
# print in class functions, not in driver

from __future__ import division
from scipy import signal, interpolate, fftpack
import numpy as np
import sys, os, math

PATH = os.getcwd()+'/'

# better smoothing funciton; should be in np.signal, but import it locally 
# till it appears
from helper_functions import *

# peak-o-mat stuff
sys.path.append(PATH + '../peak-o-mat-1.1.9/')
from peak_o_mat.model import Model
from peak_o_mat.spec import Spec
from peak_o_mat.fit import Fit

class Spectra:
    """ Stores spectra data in Spec (peak-o-mat) format
    
    Contains
    --------   
    origin: A spec object to store the original x,y data, 
        load the file, as well as to use for the fitting routines
    
    The following y-data (the x-data remains constant for the other operations)
    bg: attempt poly background fit
    
    1. Original Data
    2. Active Data
    3. Background
    
    
    4. Number of peaks found (`num_peaks`)
    5. Model and Model string (`m` and `model_str`)
    
    """
    def __init__(self, filename):
        # import data into spec object
        print "Loading file ... " 
        self.origin = Spec(PATH + filename)
        
        
        self.x = self.origin.x  # for ease of use
        self.data_points = len(self.origin.y) 
        self.active = self.origin.y # working data set

        # make a first guess of peak width    
        self.guess_peak_width()
        

    def find_background(self, sub_range=None, poly_deg=3, smoothing=5):
        """ Attempts to find the background of the spectra, 
        and updates the `bg` array
        
        Procedure
        ---------
        1. Smooths y-data with savtzky_golay
        2. Splits data into subranges
        3. Finds y-value min of subranges, and assigns to mid x-value
        4. Fits n degree polynomial to min values
        5. Updates background data with the polynomial evaluated at x-data
        """
        print "Finding background ... " 
        
        if sub_range == None:
            sub_range = (self.data_points/5)
            
        # smooth y-data, maybe add poly order as a parameter
        # smooth_y = savitzky_golay(self.active,smoothing,3)
        # smooth_y = self.filter_high_freq(self.active)
    
        # find # of sub-ranges/intervals
        intervals = int(math.ceil(self.data_points/sub_range))
        
        # find min of each subinterval and place it at mid point of 
        #bg_x = np.zeros(intervals)
        #bg_y = np.zeros(intervals)
        bg_x = []
        bg_y = []
        for i in range(intervals):
            pos = i*sub_range
            half_range = int(math.floor(sub_range/2))
            
            # if y-value is below 20% add to list
            
            bg_x_test = self.x[pos+half_range-1]
            bg_y_test = min(smooth_y[pos:(pos+sub_range)])
            
            if bg_y_test/self.data_max < .3:
                bg_x.append(bg_x_test)
                bg_y.append(bg_y_test)
            
        #self.slmin = signal.argrelmin(sg.savitzky_golay(y_data,9,5))
        # also find the 10% of min data points in the data set and poly fit
            # them
        # find the space between peaks and fit these
            
        #bg_poly = np.polyfit(bg_x,bg_y,poly_deg)    # polynomial coefficeints
            # of fit
        #bg = np.poly1d(bg_poly) # evaluate polynomial
        #self.bg = bg(self.origin.x)
        bg_x_full = [self.x[0]] + bg_x + [self.x[-1]]
        bg_y_full = [smooth_y[0]] + bg_y + [smooth_y[-1]]
        
        # smooth this data
        #bg_y_full = signal.medfilt(bg_y_full, 5)        
        
        bg_spline = interpolate.interp1d(bg_x_full, bg_y_full, kind='quadratic')
        
        self.bg = bg_spline(self.x)
        
    def subtract_background(self):
        """ Subtract background from active spectra """
        print "Subtracting background ... "
        self.active = self.active - self.bg

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

    def build_model(self, peak_type='LO', max_width=None):
        """ Builds a peak-o-mat model of peaks in listed by index in `peak_pos`
        
        Parameters
        ----------
        peak_type : string
            Peaks can be of the following types: 
            (to setup custom peaks and more, see peak-o-mat docs)
            
            | 'LO' : symmetric lorentzian
            | 'GA' : symmetric gaussain
            | 'VO' : voigt profile
            | 'PVO' : psuedo-voigt profile
            | 'FAN' : fano lineshape        
        
        
        Peaks can veof default type lorentzian (LO). Uses some basic algorithms 
        to determine initial parameters for amplitude and fwhm (limit on fwhm 
        to avoid fitting background as peaks. """   
        
        print "Building model ... "
        
        if max_width==None:
            # peaks should be at most one tenth of the data set
            max_width = self.data_points/10
        
        background = 'CB' # since using a seperate background algorithm, always constant background 

        # start building model string        
        model_str = background
        for i in range(self.num_peaks):
            # build model string of peaks
            model_str = model_str + ' ' + peak_type + str(i+1) # number peaks from 1 not 0
        
        print "Model string is ", model_str
        test_model = Model(model_str)
        
        # find initial values for paramters of model (position,height,fwhm) 
        params = [{'const':0.0}] # again background is always zero.
        for i in self.peak_pos:
            amplitude = self.active[i]
            position = self.origin.x[i]
            
            ## seperate into a seperate function
            # Find fwhm for each peak
            half_max = amplitude/2
            left = i
            right = i
            
            # look for left and right postions of when data is below half max; 
            # also make sure index does not get out of bounds
            while (self.active[left]>half_max and left > 0 ):
                left = left - 1
            while (self.active[right] > half_max and right < (self.data_points-1)):
                right = right + 1

            # find distance between these two point
            fwhm = self.origin.x[right] - self.origin.x[left]
            
            # make sure doesn't blow up
            if fwhm > max_width:
                fwhm = max_width
            
            params.append({'amp': amplitude, 'fwhm': fwhm, 'pos': position })
        
        # set this model 
        test_model.set_parameters(params)
        
        # Evaluate the model at the same number of data points as original
        # ??? Should I increase the number of data points ??
        ## Obviously not, since this is only for plotting.
        # self.x = np.linspace(min(self.s.x),max(self.s.x),5*self.max_points)
        self.model_data = test_model.evaluate(self.origin.x)
        self.model = test_model
        
    def fit_data(self):
        """ Attempt to fit data using peak-o-mat Fit function with the 
        generated model. Updates model with fit parameters. """
        
        print "Fitting Data..."
        
        fit = Fit(Spec(self.x,self.active,"currentdata"), self.model)
        result = fit.run()
        
        print "Fit result: ", result[-1]        
        
        # update model, and model plot points
        self.model.update_from_fit(result)
        self.model_data = self.model.evaluate(self.x)
        
    def output_results(self):
        """ Output fit paramters as csv values with errors"""
        
        print "Fit results"
        
        for i in ["pos","fwhm", "area", "amp"]:
            print i
            print self.model.parameters_as_csv(selection=i,witherrors=True)
    
    # ---
    # Helper functions, might make private
    # ---     
    
    def guess_peak_width(self,max_width = 50):
        """ Find an initial guess for the peak with of the data imported, 
        use in peak finding and model buildings and other major functions, 
        probably should call in the constructor
        
        Parameters
        ----------
        max_width : int
            Max width of peaks to search for
        
        Notes
        -------
        Locates the max value in the data
        Finds the peak width associated with this data
        """
        self.data_max = max(self.active)
        self.data_max_pos = np.argmax(self.active)
        self.test_peak_width = self.find_fwhm(self.data_max_pos)
        
        print "Peak width of about ", self.test_peak_width
        
        
    def remove_spikes(self,strength = 0.5):
        """ Attempts to remove spikes in active set using a simple test of 
        the pixels around it. Fractional value of strength needed."""
        print "Removing spikes..."

        mean = lambda x,y: (x+y)/2
        
        y = self.active
        data_max = max(y)
        for i in range(1,len(y)-1):
            if (np.abs( y[i] -  mean(y[i-1],y[i+1]) )/ data_max ) > strength:
                y[i] = mean(y[i-1],y[i+1])
        
        self.active = y
        
    def find_fwhm(self,position):
        """ Find the fwhm of a point using a very simplisitic algorigthm. 
        Could return very large width. """
        left = position
        right = position
        half_max = self.active[position]/2
        
        # change max_width to function of data set
        
        # make sure index does not get out of bounds
        while (self.active[left] > half_max and left > 0 ):
            left = left-1
        while (self.active[right] > half_max and right < (self.data_points-1)):
            right = right + 1
            
        # left = find index to left when height is below half_max
        # right same as above
        # find distance between these two point
        fwhm = self.x[right] - self.x[left]
            
        return fwhm
        
    def filter_high_freq(self, data):
        """ Filter high frequency data using fft 
        
        Parameters
        ----------
        data : np.array
            Data to filter
            
        Returns
        -------
        data : np.array
            Data with high frequcny components removed
        """
        trans = fftpack.fft(data)
        trans[2:] = np.zeros(len(trans)-2)
        return fftpack.ifft(trans)