import numpy as np
import rampy as rp
from scipy.signal import butter, filtfilt
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt


def subtract_background_lp_filter(cutoff=None, order=2):
    """ Attempts to find the background of the spectra using low
    pass filter

    Arguments
    ---------
    cutoff (float) [default: 2/len(y)]
        cutoff frequency at which to filter the data
    order (int) [default: 2]
        filter order

    Returns
    -------
    bg : array
        background spectrum

    """
    print("Finding background ... ")
    if cutoff is None:
        cutoff = 2 / len(y)
    # use the low pass filter to find the background
    bg = butter_lp_filter(cutoff, order)
    return bg


def find_background_lp_filter(y, thresh=0.1, fs=8096, order=5, cutoff=20, plots=False):
    """
    Threshold of data to consider for background
    Sampling frequency, order, cutoff for butter filter
    """
    y_bg = np.clip(y, 0, thresh * np.max(y))

    # butter(order, cutoff)
    bb, ab = butter(5, 20.0 * 2.0 / fs, btype="low")
    y_fl = filtfilt(bb, ab, y_bg)

    pxx, freqs = mlab.psd(y, Fs=fs)

    if plots:
        fig, ax = plt.subplots()
        ax.loglog(freqs, np.sqrt(pxx), "r-")
        ax.set_xlabel("Freq (Hz)")

    return y_fl


def subtract_background(x, y, method="arPLS", lam=10 ** 6):
    """
    Subtract baseline using defaults
    """
    bir = np.array([[np.min(x), np.max(x)]])
    yc, bg = rp.baseline(x, y, bir, method, lam=lam)

    return yc.T[0]
