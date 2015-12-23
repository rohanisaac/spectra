"""
Spectra class
-------------
Analyze spectral data using combination of numpy, scipy, lmfit and
some simple algorithms

author: Rohan Isaac
"""

# keep only one spec object, lots of y-data
# print in class functions, not in driver

from __future__ import division
import numpy as np
import re
from scipy import signal
from lmfit.models import LorentzianModel,PolynomialModel


class Spectra:
    """
    Primary spectra class that stores various stages of data processing for a
    single data set (two column x,y data) in np.ndarray formats. For details
    see the constructor.
    """
    def __init__(self, *args):
        """
        Create an object of spectra

        1 argument : path
            - path to two column text data
        2 arguments : x , y
            - numpy arrays or lists of x and y data. should be of equal length

        Usage
        -----
        >>> import spectra as sp
        >>> sp_obj = sp.Spectra("/this/is/the/path.txt")
        >>> dat = np.genfromtxt('/path/to/file.txt')
        >>> x_dat = dat[:,0]
        >>> y_dat = dat[:,1]
        >>> sp_obj2 = sp.Spectra(x_dat, y_dat)

        Member Functions
        ----------------
        For details see docs for individual function

        find_background()
        subtract_background()
        remove_spikes()
        guess_peak_width()
        find_peaks()
        build_model()
        fit_data()
        output_results()

        Plotting
        --------

        plt.plot(S.x,S.y,'-') # active data
        plt.plot(S.x,S.bg,'-r') # backround
        plt.plot(S.x,S.md,'r-') # model data
        plt.plot(S.x[S.peak_pos],S.y[S.peak_pos],'oy')
        plt.plot(Sx,S.md,'r-') # fitted y-data


        Data Members
        ------------

        ox : original x-data
        oy : original y-data
        x : active x-data
        y : active y-data
        bg : found background
        num_peaks : number of peaks found
        peak_list : indices of found peaks

        """
        # import data into spec object
        print "Loading file ... "
        print args[0]
        if len(args) == 1:
            file_name = args[0]
            self.x ,self.y = self.getxy(file_name)
        elif len(args) == 2:
            self.x, self.y = args


        self.num_points = len(self.y)
        # make a first guess of peak width
        # also updates max, and max position
        self.guess_peak_width()

    def getxy(self, file_name, headers = False):
        """Extracts x and y data numpy arrays from passed filename.

        Arguments
        ---------
        file_name: string
            full/relative path to file to extract data from

        Returns
        -------
        if headers == True:
            x,y,xlab,ylab (np.array, np.array, str, str)

        else:
            x,y
        where:
            x,y: 1-d numpy arrays containing first and second column of data present in file
            xlab, ylab: header information about columns
        """
        line_pos = 0 # active line number
        start_pos = 0
        end_pos = 0
        dat_read = False # has data been read
        header = ''

        with open(file_name, 'rb') as fil:
            # find header and footer positions
            for lin in fil:
                line_pos += 1
                # if line contains any of the alphabet (except e for exponents, not data)
                if re.search('[a-df-zA-DF-Z]', lin):
                    if not dat_read:
                        # before data has been read, set start of data pos
                        start_pos = line_pos
                        header = lin
                    if dat_read and end_pos==0:
                        # after data had been read and before end position has been set, set end pos
                        end_pos = line_pos

                # if data line and data has not been read
                elif not dat_read:
                    # find seperator
                    if re.search('\t', lin):
                        sep_char = '\t'
                    elif re.search(',', lin):
                        sep_char = ','
                    else:
                        print "Unknown separator character"
                    # now we know what separator is for the data
                    dat_read = True

                # data line and we already know what the separator character is
                else:
                    continue

        # if we didn't find an end position
        if end_pos == 0:
            skip_foot = 0
        else:
        # if we did compute it
            skip_foot = line_pos - end_pos

        xlab,ylab = ('','')
        # find header row if exists
        header_lst = header.split(sep_char)
        #print headerlst
        if len(header_lst) >= 2:
            xlab, ylab = header_lst[0:2]

        # attempt to load into numpy array, see what happens
        fdat = np.genfromtxt(file_name, delimiter=sep_char, skip_header=start_pos, skip_footer = skip_foot )

        if headers == True:
            return fdat[:,0],fdat[:,1],xlab,ylab
        else:
            return fdat[:,0],fdat[:,1]


    def find_background(self, window_size, order):
        """ Attempts to find the background of the spectra,

        Updates
        -------
        bg : array
            background spectrum

        Arguments
        ---------
        window_size : int
            [Savitzky Golay] the length of the window. Must be an odd integer
            number.
        order : int
            [Savitzky Golay] the order of the polynomial used in the filtering.
            Must be less then `window_size` - 1.
        """

        print "Finding background ... "

        # smooth y-data
        self.bg = savitzky_golay(self.y,window_size,order)
        return

    def subtract_background(self):
        """ Subtract background from active spectra """
        print "Subtracting background ... "
        self.oy = self.y
        self.y = self.y - self.bg

    def find_peaks(self, lower=None, upper=None, threshold=5, limit=20):
        """ Find peaks in active data set using continuous wavelet
        transformation

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

        """
        #smooth y-data
        y = signal.savgol_filter(self.y,window_length=25,polyorder=2)
        #y = self.y
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
        """ Builds a lmfit model of peaks in listed by index in `peak_pos`
        Uses some basic algorithms to determine initial parameters for amplitude
        and fwhm (limit on fwhm to avoid fitting background as peaks).

        Parameters
        ----------
        peak_type : string (default='lorentizian')
            Peaks can be of the following types:

            - 'LO' : symmetric lorentzian
            - 'GA' : symmetric gaussain

        max_width : int (default = total points/10)
            max width (in data points) that peak fitted can be

        Updates
        -------
        md : ndarray
            y-data of model
        model : model object


        """
        norm = 0.398942280401
        x = self.x
        y = self.y
        peak_guess = self.x[self.peak_pos]
        print "Building model ... "

        # start with polynomial background
        model = PolynomialModel(2, prefix='bg_')
        pars = model.make_params()

        # add lorentizian peak for all peaks
        for i, peak in enumerate(peak_guess):
	    temp_model = LorentzianModel(prefix='p%s_' % i)
	    pars.update(temp_model.make_params())
	    model += temp_model


        #print model
        #pars = model.make_params()
        print pars
        #print re.sub(',','\n',pars.viewvalues())

        # set inital background as flat line at zeros
        pars['bg_c0'].set(0)
        pars['bg_c1'].set(0)
        pars['bg_c2'].set(0)

        # give values for other peaks
        for i,peak in enumerate(self.peak_pos):
            print x[peak],
            pars['p%s_center' % i].set(x[peak], min=x[peak]-5,max=x[peak]+5)
            pars['p%s_fwhm' % i].set(2,min=0,max=5)
            pars['p%s_amplitude' % i].set(y[peak]*norm, min=0, max=2*max(y))
            print y[peak]

        self.pars = pars
        self.model = model

    def fit_data(self):
        """ Attempt to fit data using lmfit fit function with the
        generated model. Updates model with fit parameters. """

        print "Fitting Data..."
        out = self.model.fit(self.y,self.pars,x=self.x)
        self.out = out

    def output_results(self):
        """ Output fit paramters as summary table"""

        print(out.fit_report())

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

        self.data_max = max(self.y)
        self.data_max_pos = np.argmax(self.y)
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
        return
        #self.base.y = y # don't need because of python lazy copying ??

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
        half_max = self.y[position]/2

        # change max_width to function of data set

        # make sure index does not get out of bounds
        while (self.y[left] > half_max and left > 0 ):
            left = left-1
        while (self.y[right] > half_max and right < (self.num_points-1)):
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
        ft = fftpack.fft(self.base.y)
        ft[cof:] = np.zeros(len(trans)-2)
        self.base.y = fftpack.ifft(ft)

    def calibrate_x(self, m, b):
        """ Applies a linear correction to the x-values """
        # Need to change the active data set
        # Save the old data etc.

        pass
