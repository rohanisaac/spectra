from scipy import signal

def smooth_data(self, window_size=25, order=2):
    """ Smooths data using savgol_filter """
    self.y_smooth = signal.savgol_filter(self.y, window_size, order)

def butter_lp_filter(self, cutoff=None, order=2):
    """
    Performs a low pass butterworth filter on the data with cutoff
    frequency that defaults to 2/len(y)

    Arguments
    ---------
    cutoff (float) [default: 2/len(y)]
        cutoff frequency at which to filter the data
    order (int) [default: 2]
        filter order

    Returns
    -------
    (array)
        low pass filtered array same size as data
    """
    if cutoff is None:
        cutoff = 2 / len(self.y)
    B, A = signal.butter(order, cutoff, output='ba')
    return signal.filtfilt(B, A, self.y)


def butter_lowpass(cutoff, fs, order=5):
    nyq = 0.5 * fs
    normal_cutoff = cutoff / nyq
    b, a = signal.butter(order, normal_cutoff, btype='low', analog=False)
    return b, a


def butter_lowpass_filter(data, cutoff, fs, order=5):
    """
    Maximally flat frequency response in the pass band

    Not actuall a zero phase filter, use filtfilt instead
    http://scipy-cookbook.readthedocs.io/items/FiltFilt.html
    https://github.com/scipy/scipy-cookbook/blob/master/ipython/FIRFilter.ipynb
    That or since it is a FIR filter (?check?), the delay will affect all the
    components the same way. The delay of a filter of length M equals (M-1)/2

    https://oceanpython.org/2013/03/11/signal-filtering-butterworth-filter/

    For details see: https://tomroelandts.com/articles/how-to-create-a-simple-low-pass-filter
    http://dsp.stackexchange.com/questions/2864/how-to-write-lowpass-filter-for-sampled-signal-in-python
    http://scipy-cookbook.readthedocs.io/items/FIRFilter.html
    """
    b, a = butter_lowpass(cutoff, fs, order=order)
    y = signal.lfilter(b, a, data)
    return y

## Use digital filter, regularly sampled data
#
def wicker():
    pass


def savgol_filter():
    """Just use scipy"""
    pass

# https://github.com/scipy/scipy-cookbook/blob/master/ipython/ButterworthBandpass.ipynb
def butterworth_bandpass(b, a, data):
    return y

## for fitting sin's gaussian etc:

http://scipy-cookbook.readthedocs.io/items/FittingData.html


## million ways of fitting a line

http://scipy-cookbook.readthedocs.io/items/LinearRegression.html


### plz read

http://scipy-cookbook.readthedocs.io/items/robust_regression.html


### Probs not relevant most cases
https://en.wikipedia.org/wiki/Random_sample_consensus
