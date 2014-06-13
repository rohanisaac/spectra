"""
Spectra class
------- -----
Analyze spectral data using combination of numpy, scipy, peak-o-mat and some simple algorithms
@author: Rohan Isaac
"""

from __future__ import division
from scipy import signal
import numpy as np
import sys, os, math

PATH = os.getcwd()+'/'

# better smoothing funciton; should be in np.signal, but import it locally till it appears
from savitzky_golay import savitzky_golay

# peak-o-mat stuff
sys.path.append(PATH + '../peak-o-mat-1.1.9/')
from peak_o_mat.model import Model
from peak_o_mat.spec import Spec
from peak_o_mat.fit import Fit

class Spectra:
    """ Stores spectra data in Spec (peak-o-mat) format
    
    Contains
    --------   
    (Spec objects)
    1. Original Data
    2. Active Data
    3. Background
    
    
    4. Number of peaks found (`num_peaks`)
    5. Model and Model string (`m` and `model_str`)
    
    """
    def __init__(self, filename):
        # import data into spec object
        self.origin = Spec(PATH + filename)
        self.data_points = len(self.origin.y)
        
        # other spec objects        
        self.active = Spec(self.origin.x,self.origin.y,'active')
        self.bg = Spec(self.origin.x,np.zeros(self.data_points),'background')
        
        # basic model attributes
        
        # really only need one spec object, can use y-lists for the rest
        

    def find_background(self,sub_range=None,poly_deg=3,smoothing=5):
        """ Attempts to find the background of the spectra, and updates the `bg` Spec object 
        
        Procedure
        ---------
        1. Smooths data with savtzky_golay
        2. Splits data into subranges
        3. Finds y-value min of subranges, and assigns to mid x-value
        4. Fits n degree polynomial to min values
        5. Returns original x-data and y-values of x-data evaluated with the polynomial
        """
        if sub_range == None:
            sub_range = (self.max_points/20)
            
        # smooth y-data, maybe add poly order as a parameter
        smooth_y = savitzky_golay(self.origin.y,smoothing,3)
    
        # find # of sub-ranges/intervals
        intervals = int(math.ceil(len(self.s.y)/sub_range))
        
        # find min of each subinterval and place it at mid point of 
        bg_x = np.zeros(intervals)
        bg_y = np.zeros(intervals)
        for i in range(intervals):
            pos = i*sub_range
            half_range = int(math.floor(sub_range/2))
            bg_x[i] = self.origin.x[pos+half_range-1]
            bg_y[i] = min(smooth_y[pos:(pos+sub_range)])
            
        #self.slmin = signal.argrelmin(sg.savitzky_golay(y_data,9,5))
        # also find the 10% of min data points in the data set and poly fit them
        # find the space between peaks and fit these
            
        bg_poly = np.polyfit(bg_x,bg_y,poly_deg)    # polynomial coefficeints of fit
        bg = np.poly1d(bg_poly) # evaluate polynomial
        self.bg.y = bg(self.origin.x)
        
    def subtract_background(self):
        """ Subtract background from active spectra """
        self.active = self.active - self.bg

    def find_peaks(self,lower=None,upper=None,threshold=5,limit=20):
        """ Find peaks in actve data set using continous wavelet transformation from `scipy.signal`
        
        Variables modified
        --------- --------
        peak_pos: list of indicies associated with peak positions
        num_peaks: # of peaks found
        
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

        peak_pos = signal.find_peaks_cwt(self.active.y,np.arange(lower,upper))
        
        # remove peaks that are not above the threshold.
        peak_pos = [i for i in peak_pos if (self.active.y[i]/self.data_max) > (threshold/100)]  
        
        # remove peaks that are pretty close together?
        # only use the top limit peaks
        # first sort the peak list by size of y
        # print np.argsort(self.s.y[self.peak_pos])
        
        self.peak_pos = peak_pos[1:limit]
        self.num_peaks = len(self.peak_pos)

    def build_model(self, peak_type='LO', max_width=None):
        """ Builds a peak-o-mat model of peaks in listed by index in `peak_pos`
        
        Peaks can be of the following types: (to setup custom peaks and more, see peak-o-mat docs)
        LO: symmetric lorentzian
        GA: symmetric gaussain
        VO: voigt profile
        PVO: psuedo-voigt profile
        FAN: fano lineshape        
        
        
        Peaks can veof default type lorentzian (LO). Uses some basic algorithms to determine initial parameters for amplitude and fwhm (limit on fwhm to avoid fitting background as peaks. """   
        if max_width==None:
            # peaks should be at most one tenth of the data set
            max_width = self.data_points/10
        
        background = 'CB' # since using a seperate background algorithm, always constant background 

        # start building model string        
        model_str = background
        for i in range(self.num_peaks):
            # build model string of peaks
            model_str = model_str + ' ' + peak_type + str(i+1) # number peaks from 1 not 0
        
        test_model = Model(model_str)
        
        # find initial values for paramters of model (position,height,fwhm) 
        params = [{'const':0.0}] # again background is always zero.
        for i in self.peak_pos:
            amplitude = self.active.y[i]
            position = self.active.x[i]
            
            ## seperate into a seperate function
            # Find fwhm for each peak
            half_max = amplitude/2
            left = i
            right = i
            
            # look for left and right postions of when data is below half max; 
            # also make sure index does not get out of bounds
            while (self.active.y[left]>half_max and left > 0 ):
                left = left - 1
            while (self.active.y[right] > half_max and right < (self.max_points-1)):
                right = right + 1

            # find distance between these two point
            fwhm = self.active.x[right] - self.active.x[left]
            
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
        self.model_x = self.origin.x
        self.model_y = self.test_model.evaluate(self.model_x)
        
        self.model = test_model
        
    def fit_data(self):
        """ Attempt to fit data using peak-o-mat Fit function with the generated model. 
        Updates model with fit parameters. """
        
        fit = Fit(self.active, self.model)
        result = fit.run()
        
        # update model, and model plot points
        self.model.update_from_fit(result)
        self.model_y = self.model.evaluate(self.model_x)
        
    def output_results(self):
        """ Output fit paramters as csv values with errors"""
        for i in ["pos","fwhm", "area", "amp"]:
            print i
            print self.model.parameters_as_csv(selection=i,witherrors=True)
    
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
