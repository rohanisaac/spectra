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
print "Base path is: ", PATH
# better smoothing funciton; should be in np.signal, but import it locally 
# till it appears
from helper_functions import *

# peak-o-mat stuff
sys.path.append(PATH + 'peak-o-mat-1.1.9/')
from peak_o_mat.model import Model
from peak_o_mat.spec import Spec
from peak_o_mat.fit import Fit

class Spectra:
    """ 
    Main class that stores various stages of data processing for a single data
    set (two column x,y data) in peak-o-mat formats. For details see the
    constructor.
    """
    def __init__(self, *args):
        """
        Create an object of spectra 
        
        1 argument : path
            - path to two column text data
        2 arguments : x , y 
            - numpy arrays or lists of x and y data. should be of equal lenght
        
        Usage
        -----
        >>> import spectra as sp
        >>> sp_obj = sp.Specta("/this/is/the/path.txt")
        
        Member Functions
        ----------------
        For details see docs for individual function
        
        find_background()
        subtract_background()
        remove_spikes()
        guess_peak_width(max_width=50)
        find_peaks(limit=30)
        build_model() 
        fit_data()
        output_results()
        
        Plotting
        --------
        
        plt.plot(S.x,S.active,'-')
        plt.plot(S.x,S.bg,'-r')
        plt.plot(S.x,S.model_data,'r-')    
        plt.plot(S.x[S.peak_pos],S.active[S.peak_pos],'oy')
        plt.plot(S.x,S.model_data,'r-')
        
    
        Data Members
        ------------
        
        base : spec object
            Store the active x,y data, loads the file, 
            used for fitting routines
        oy : original y-data
        bg : background y-data 
        md : model y-data
        fd : fitted model y-data
        ox : original x-data
        xc : corrected x-data
        num_points : int
            number of data points in 
        num_peaks : int
            number of 
        peak_pos : list
            x positions of peaks
        m : model object
            peak-o-mat model used in fitting
        model_str : string
            model string used for building peak-o-mat model
        data_max : int
            max of y-data
        data_max_pos : int
            index associated with max data
        
        """
        # import data into spec object
        print "Loading file ... "
        print args[0]
        if len(args) == 1:
            file_name = args[0]
            self.base = Spec(PATH + file_name)
        elif len(args) == 2:
            x, y = args
            self.base = Spec(x,y,'data')
        
        # save original data
        self.ox = self.base.x  
        self.oy = self.base.y
        
        self.num_points = len(self.base.y)

        # make a first guess of peak width
        # also updates max, and max position    
        self.guess_peak_width()   

    def find_background(self, sub_range=None, window_size=5, order=3,
        bg_thresh=0.3):
        """ Attempts to find the background of the spectra, 
        
        Updates
        -------
        bg : array
            background spectrum
        
        Arguments
        ---------
        sub_range : int (default: points/5)
             size of range to split data into for background search
        window_size : int
            [Savitzky Golay] the length of the window. Must be an odd integer
            number.
        order : int
            [Savitzky Golay] the order of the polynomial used in the filtering.
            Must be less then `window_size` - 1.
        bg_thresh : float (default=0.3 i.e 30%)
            fraction that background should be below 
        
        Procedure
        ---------
        1. Smooths y-data with savtzky_golay
        2. Splits data into subranges
        3. Finds y-value min of subranges, and assigns to mid x-value
        4. Fits n degree polynomial to min values
        5. Updates background data with the polynomial evaluated at x-data
        """
        x = self.base.x
        
        print "Finding background ... " 
        
        if sub_range == None:
            sub_range = (self.num_points/5)
            
        # smooth y-data
        smooth_y = savitzky_golay(self.base.y,window_size,order)
    
        # find # of sub-ranges/intervals
        intervals = int(math.ceil(self.num_points/sub_range))
        
        # for each subinterval find min and place at mid point
        bg_x = []
        bg_y = []
        for i in range(intervals):
            # find midpoint
            pos = i*sub_range
            half_range = int(math.floor(sub_range/2))
            
            # add min and corresponding x to list
            bg_x_test = x[pos+half_range-1]
            bg_y_test = min(smooth_y[pos:(pos+sub_range)])
            
            # if y-value is below 30% add to list
            if bg_y_test/self.data_max < .3:
                bg_x.append(bg_x_test)
                bg_y.append(bg_y_test)
            
        # add endpoints of data as well
        bg_x_full = [x[0]] + bg_x + [x[-1]]
        bg_y_full = [smooth_y[0]] + bg_y + [smooth_y[-1]]
        
        # interpolate this data set       
        bg_spline = interpolate.interp1d(bg_x_full, bg_y_full, kind='quadratic')
        self.bg = bg_spline(x)
        
    def subtract_background(self):
        """ Subtract background from active spectra """
        print "Subtracting background ... "
        self.base = self.base - Spec(x,bg_spline(x),"background")

    def find_peaks(self, lower=None, upper=None, threshold=5, limit=20):
        """ Find peaks in active data set using continuous wavelet 
        transformation from `scipy.signal`
        
        Options
        -------
        lower: int
            lower limit of peak size to search for
        upper: int 
            upper limit for same (in index values)
        threshold: float (default=5)
            min percent of max to count as a peak (eg 5 = only peaks above 5 
            percent reported)
        limit: int 
            max limit of peaks to report (sorted by intensity)
        
        Updates
        -------
        peak_pos : list
            indices associated with peak positions
        num_peaks : int
            number of peaks found
            
        Todo
        ----
        remove peaks that are pretty close together?
              
        """
        # smooth y-data
        # smooth_y = savitzky_golay(self.base.y,window_size=25,order=2)
        y = self.base.y
        print "Looking for peaks ... "  
        
        if lower == None:
            lower = self.test_peak_width*1
        if upper == None:
            upper = self.test_peak_width*5

        peak_pos = signal.find_peaks_cwt(y,np.arange(lower,upper),min_snr=2)
        
        print "Found %s peaks at %s" % (len(peak_pos),peak_pos)
        
        # remove peaks that are not above the threshold.
        peak_pos = [i for i in peak_pos if (y[i]/self.data_max) > (threshold/100)]  
        
        print "After filtering out peaks below ", threshold, \
            "percent, we have ", len(peak_pos), " peaks."
        
        # only use the most intense peaks, zip two lists together, 
        # make the y-values as the first item, and sort by it (descending)
        peak_pos = [y for (x,y) in sorted(zip(y[peak_pos],peak_pos),reverse=True)]
        
        self.peak_pos = peak_pos[0:limit]
        self.num_peaks = len(self.peak_pos)
        
        print "Using ", self.num_peaks, " peaks at ", self.peak_pos 

    def build_model(self, peak_type='LO', max_width=None):
        """ Builds a peak-o-mat model of peaks in listed by index in `peak_pos`
        Uses some basic algorithms to determine initial parameters for amplitude
        and fwhm (limit on fwhm to avoid fitting background as peaks).
        
        Parameters
        ----------
        peak_type : string (default='LO')
            Peaks can be of the following types: 
            (to setup custom peaks and more, see peak-o-mat docs)
            
            | 'LO' : symmetric lorentzian
            | 'GA' : symmetric gaussain
            | 'VO' : voigt profile
            | 'PVO' : psuedo-voigt profile
            | 'FAN' : fano lineshape        
        
        max_width : int (default = total points/10)
            max width (in data points) that peak fitted can be
            
        Updates
        -------
        md : ndarray
            y-data of model
        model : model object
            
      
        """   
        
        print "Building model ... "
        
        if max_width==None:
            # peaks should be at most one tenth of the data set
            max_width = self.num_points/10
        
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
            amplitude = self.base.y[i]
            position = self.base.x[i]
            
            ## seperate into a seperate function
            # Find fwhm for each peak
            half_max = amplitude/2
            left = i
            right = i
            
            # look for left and right postions of when data is below half max; 
            # also make sure index does not get out of bounds
            while (self.base.y[left]>half_max and left > 0 ):
                left = left - 1
            while (self.base.y[right] > half_max and right < (self.num_points-1)):
                right = right + 1

            # find distance between these two point
            fwhm = self.base.x[right] - self.base.x[left]
            
            # make sure doesn't blow up
            if fwhm > max_width:
                fwhm = max_width
            
            params.append({'amp': amplitude, 'fwhm': fwhm, 'pos': position })
        
        # set this model 
        test_model.set_parameters(params)
        
        # update class
        self.md = test_model.evaluate(self.base.x)
        self.model = test_model
        
    def fit_data(self):
        """ Attempt to fit data using peak-o-mat Fit function with the 
        generated model. Updates model with fit parameters. """
        
        print "Fitting Data..."
        
        fit = Fit(Spec(self.base.x,self.base.y,"currentdata"), self.model)
        result = fit.run()
        
        print "Fit result: ", result[-1]        
        
        # update model, and model plot points
        self.model.update_from_fit(result)
        self.fd = self.model.evaluate(self.base.x)
        
    def output_results(self):
        """ Output fit paramters as csv values with errors"""
        
        print "Fit results"
        
        for i in ["pos","fwhm", "area", "amp"]:
            print i
            print self.model.parameters_as_csv(selection=i,witherrors=True)
    
    # ---
    # Helper functions
    # ---     
    
    def guess_peak_width(self, max_width = None):
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
        
        Updates
        -------
        data_max : float
            max intensity of y-data
        data_max_pos : int
            index of max data
        test_peak_width :
            guess for peak width of data
        
        """
        if max_width == None:
            max_width = self.num_points/5
        
        self.data_max = max(self.base.y)
        self.data_max_pos = np.argmax(self.base.y)
        self.test_peak_width = self.find_fwhm(self.data_max_pos)
        
        print "Peak width of about ", self.test_peak_width
        
        
    def remove_spikes(self,strength = 0.5):
        """ Attempts to remove spikes in active set using a simple test of 
        the pixels around it. Fractional value of strength needed.
        
        Arguments
        ---------
        strength : float (default: 0.5)
            ratio of data point over average of surrounding points must be to
            count as a spike 
        """
        print "Removing spikes..."
        
        ## !!! Try scipy.signal.medfilt

        mean = lambda x,y: (x+y)/2
        
        y = self.base.y
        data_max = max(y)
        for i in range(1,len(y)-1):
            if (np.abs( y[i] -  mean(y[i-1],y[i+1]) )/ data_max ) > strength:
                y[i] = mean(y[i-1],y[i+1])
        
        self.base.y = y
        
    def find_fwhm(self,position):
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
        half_max = self.base.y[position]/2
        
        # change max_width to function of data set
        
        # make sure index does not get out of bounds
        while (self.base.y[left] > half_max and left > 0 ):
            left = left-1
        while (self.base.y[right] > half_max and right < (self.num_points-1)):
            right = right + 1
            
        # left = find index to left when height is below half_max
        # right same as above
        # find distance between these two point
        fwhm = self.base.x[right] - self.base.x[left]
            
        return fwhm
        
    def filter_high_freq(self, cof=5):
        """ Filter high frequency y-data using 1-D fft 
        
        Parameters
        ----------
        cof : int (default=5)
            Coefficient above which to truncate Fourier series. 
            
        """
        ft = fftpack.fft(self.base.y)
        ft[cof:] = np.zeros(len(trans)-2)
        self.base.y = fftpack.ifft(ft)
		
	def linear_calibrate(self, m, b):
		""" Calibrate the x-data with a linear correction function 
		of the form (new x values) = m (old x-values) + b, 
		
		Parameters
		----------
		m : float
			Slope of linear correction
		b : float
			y-intercept
            
        Updates
        -------
        xc : 
            corrected x-data
		"""
		self.xc = (self.base.x*m) + b 

