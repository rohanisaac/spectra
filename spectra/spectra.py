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
import pandas as pd
from scipy import signal
from lmfit.models import LorentzianModel, PolynomialModel


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

        Data Members
        ------------
        x : x-data
        y : y-data
        out.init_fit : model y-data
        out.best_fit : fit y-data
        peak_pos : index of peaks found

        """
        # import data into spec object
        print "Loading file ... "
        if len(args) == 1:
            file_name = args[0]
            print file_name
            self.x, self.y = self.getxy(file_name)
        elif len(args) == 2:
            self.x, self.y = args

        self.num_points = len(self.y)
        # make a first guess of peak width
        # also updates max, and max position
        self.guess_peak_width()

    def getxy(self, file_name, headers=False):
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
            x,y: 1-d numpy arrays containing first and second column of data
            present in file
            xlab, ylab: header information about columns
        """
        line_pos = 0  # active line number
        start_pos = 0
        end_pos = 0
        dat_read = False  # has data been read
        header = ''

        with open(file_name, 'rb') as fil:
            # find header and footer positions
            for lin in fil:
                line_pos += 1
                # if line contains any of the alphabet (except e for exponents,
                # not data)
                if re.search('[a-df-zA-DF-Z]', lin):
                    if not dat_read:
                        # before data has been read, set start of data pos
                        start_pos = line_pos
                        header = lin
                    if dat_read and end_pos == 0:
                        # after data had been read and before end position has
                        # been set, set end pos
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

        xlab, ylab = ('', '')
        # find header row if exists
        header_lst = header.split(sep_char)
        # print headerlst
        if len(header_lst) >= 2:
            xlab, ylab = header_lst[0:2]

        # attempt to load into numpy array, see what happens
        fdat = np.genfromtxt(
            file_name, delimiter=sep_char, skip_header=start_pos,
            skip_footer=skip_foot)

        if headers:
            return fdat[:, 0], fdat[:, 1], xlab, ylab
        else:
            return fdat[:, 0], fdat[:, 1]

    def smooth_data(self, window_size=25, order=2):
        """ Smooths data using savgol_filter """
        self.y_smooth = signal.savgol_filter(self.y, window_size, order)

    def find_background(self):
        """ Attempts to find the background of the spectra,

        Updates
        -------
        bg : array
            background spectrum

        Arguments
        ---------
        """

        print "Finding background ... "
        pass

    def subtract_background(self):
        """ Subtract background from active spectra """
        print "Subtracting background ... "
        pass

    def find_peaks(self, width=None, w_range=5, threshold=5, limit=20,
                   smooth=False):
        """ Find peaks in active data set using continuous wavelet
        transformation

        Options
        -------
        width: float
            estimate of peak size
        w_range: int (default=5)
            number of widths to use for wavelet transformation
            NOTE: more increases computational time (?prob), does not seem to
            affect estimation siginficantly.
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
        print "Looking for peaks ... "

        if smooth:
            try:
                # see if smoothed data already exists
                self.y_smooth
            except:
                # if it doesn't make it with defaults
                self.smooth_data()
            y = self.y_smooth
        else:
            y = self.y

        x = self.x
        xscale = len(x)/(max(x)-min(x))

        if width is None:
            width = self.test_peak_width
        else:
            # if width is given here, use it everywhere else
            self.test_peak_width = width

        lower = width * xscale * 0.75
        upper = width * xscale * 1.25

        peak_pos = signal.find_peaks_cwt(y, np.linspace(lower, upper, w_range))

        print "Found %s peaks at %s" % (len(peak_pos), peak_pos)

        # remove peaks that are not above the threshold.
        peak_pos = [i for i in peak_pos if
                    (y[i]/self.data_max) > (threshold/100)]

        print "After filtering out peaks below ", threshold, \
            "percent, we have ", len(peak_pos), " peaks."

        # only use the most intense peaks, zip two lists together,
        # make the y-values as the first item, and sort by it (descending)
        peak_pos = [y1 for (x1, y1) in sorted(zip(y[peak_pos], peak_pos),
                    reverse=True)]

        self.peak_pos = sorted(peak_pos[0:limit])
        self.num_peaks = len(self.peak_pos)

        print "Using ", self.num_peaks, " peaks at ", self.peak_pos

    def build_model(self, peak_type='LO', max_width=None, bg_ord=2):
        """ Builds a lmfit model of peaks in listed by index in `peak_pos`
        Uses some basic algorithms to determine initial parameters for
        amplitude and fwhm (limit on fwhm to avoid fitting background as peaks)

        Parameters
        ----------
        peak_type : string (default='lorentizian')
            Peaks can be of the following types:

            - 'LO' : symmetric lorentzian
            - 'GA' : symmetric gaussain

        max_width : int (default = total points/10)
            max width (in data points) that peak fitted can be

        bg_ord: int
            order of the background polynomial
            0: constant, 1: linear, ...

        Updates
        -------
        pars : model parameters
        model : model object


        """
        x = self.x
        y = self.y
        pw = self.test_peak_width
        peak_guess = self.x[self.peak_pos]
        print "Building model ... "

        # start with polynomial background
        # second order
        model = PolynomialModel(bg_ord, prefix='bg_')
        pars = model.make_params()

        # add lorentizian peak for all peaks
        for i, peak in enumerate(peak_guess):
            temp_model = LorentzianModel(prefix='p%s_' % i)
            pars.update(temp_model.make_params())
            model += temp_model

        # print model
        # pars = model.make_params()
        # print re.sub(',','\n',pars.viewvalues())

        # set inital background as flat line at zeros
        for i in range(bg_ord+1):
            pars['bg_c%i' % i].set(0)

        # give values for other peaks
        for i, peak in enumerate(self.peak_pos):
            print 'Peak %i: pos %s, height %s' % (i, x[peak], y[peak])
            # could set bounds #, min=x[peak]-5, max=x[peak]+5)
            pars['p%s_center' % i].set(x[peak])
            pars['p%s_sigma' % i].set(pw/2, min=pw*0.25, max=pw*2)
            # here as well #, min=0, max=2*max(y))
            pars['p%s_amplitude' % i].set(y[peak]*(pw/2)*np.pi)

        self.pars = pars
        self.model = model

    def fit_data(self):
        """ Attempt to fit data using lmfit fit function with the
        generated model. Updates model with fit parameters. """

        print "Fitting Data..."
        out = self.model.fit(self.y, self.pars, x=self.x)
        print out.fit_report(show_correl=False)
        self.out = out

    def output_results(self, filename=None, pandas=False):
        """ Return fit paramters and standard error, modified from lmfit
        class. Can output same data to file if passed file path

        Options
        -------
        filename : string (default: None)
            If filename is given write data to that file as csv

        pandas : bool (default: False)
            True: Return data as pandas dataframe
            False: Return data as string (csv)
        """

        params = self.out.params
        dat_out = ''
        for name in list(params.keys()):
            par = params[name]
            dat_out += '%s\t%s\t%s\n' % (name, par.value, par.stderr)
        # print dat_out
        if filename:
            with open(filename, 'w') as f:
                f.write(dat_out)
        if pandas:
            from io import StringIO
            all_data = pd.read_csv(StringIO(unicode(dat_out)),
                                   delimiter='\t',
                                   header=None,
                                   index_col=0,
                                   names=['value', 'stderr'])
            return all_data
        else:
            return dat_out
    # ---
    # Helper functions
    # ---

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

        Updates
        -------
        data_max : float
            max intensity of y-data
        data_max_pos : int
            index of max data
        test_peak_width :
            guess for peak width of data

        """
        if max_width is None:
            max_width = self.num_points/5

        self.data_max = max(self.y)
        self.data_max_pos = np.argmax(self.y)
        self.test_peak_width = self.find_fwhm(self.data_max_pos)

        print "Peak width of about %s (in x-data units)" % self.test_peak_width

    def remove_spikes(self, strength=0.5):
        """ Attempts to remove spikes in active set using a simple test of
        the pixels around it. Fractional value of strength needed.

        Arguments
        ---------
        strength : float (default: 0.5)
            ratio of data point over average of surrounding points must be to
            count as a spike
        """
        print "Removing spikes..."

        # !!! Try scipy.signal.medfilt
        def mean(x, y): return (x+y)/2

        y = self.base.y
        data_max = max(y)
        for i in range(1, len(y)-1):
            if (np.abs(y[i] - mean(y[i-1], y[i+1]))/data_max) > strength:
                y[i] = mean(y[i-1], y[i+1])
        return
        # self.base.y = y # don't need because of python lazy copying ??

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
        half_max = self.y[position]/2

        # change max_width to function of data set

        # make sure index does not get out of bounds
        while (self.y[left] > half_max and left > 0):
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
        pass

    def calibrate_x(self, m, b):
        """ Applies a linear correction to the x-values """
        # Need to change the active data set
        # Save the old data etc.

        pass
