
import numpy as np
import os
from . import find_nearest_index

def normalize(x, type='1-norm'):
    """ Returns the normalized vector

    Arguments
    ---------
    x (np.array)
        Input vector to normalize
    type (string)
        Type of normalization to compute
        '1-norm': Area under curve = 1
        '2-norm': Vector of unit length. Larger values are weighted more heavily.
        'inf-norm': Normalize to max value
        'snv': (Standard Normal Variate), sample has unit standard deviation
    Returns
    -------
    (np.array)
        Normalized vector

    References
    ----------
    http://wiki.eigenvector.com/index.php?title=Advanced_Preprocessing:_Sample_Normalization
    """

    if type == '1-norm':
        return x / np.sum(np.abs(x))
    elif type == '2-norm':
        return x / np.sum(np.square(x))
    elif type == 'inf-norm':
        return x / np.max(x)
    elif type == 'snv':
        print("Probably shouldn't be using, doing nothing")


def normalize_msc(dat):
    """ Computes the normalization of a given dataset using Multiplicative Scatter Correction,
    which attempts to correct for scaling effects and offset.

    Uses mean as the reference

    Parameters
    ----------
    dat (np.array)
        2-d array containing vectors to normalize

    Returns
    -------
    (np.array)
        Array of same size as dat normalized to mean

    References
    ----------
    http://wiki.eigenvector.com/index.php?title=Advanced_Preprocessing:_Sample_Normalization
    Citation in page
    """
    # make a blank array to hold the result
    sample, spectrum = dat.shape
    norm = np.zeros(dat.shape)

    # find the mean to use as a reference
    one = np.zeros(spectrum)
    r = np.mean(dat, axis=0)
    rc = r - np.mean(r) * one
    rcl = np.dot(rc, rc)

    for i in range(sample):
        x = dat[i, :]
        xc = x - np.mean(x) * one
        b = np.dot(rc, xc) / rcl
        norm[i, :] = xc / b + r * one

    return norm


def normalize_pq(dat):
    """
    Computes Probablistic quotient normalization on a dataset

    Parameters
    ----------
    dat (np.array)
        2-d array containing vectors to normalize

    Returns
    -------
    (np.array)
        Normalized rray of same size as input dataset

    Notes
    -----
    Probabilistic quotient normalization more robust and more accurate than
    integral normalization and vector length normalization.

    Procedure
    ---------
    1. Perform integral normalization
    2. Calculate reference spectrum (median)
    3. Calculate quotients of spectrum vs reference
    4. Divide by median quotient

    References
    ----------
    [1] Dieterle, F. et. al. Analytical Chemistry 2006 78 (13), 4281-4290
    DOI: 10.1021/ac051632c
    """
    # setup arrays
    sample, spectrum = dat.shape
    int_norm = np.zeros(dat.shape)
    norm = np.zeros(dat.shape)

    # integral normalization
    for i in range(sample):
        int_norm[i, :] = normalize(dat[i, :])

    # find the meadian to use as a reference
    mead_spec = np.median(int_norm, axis=0)

    for i in range(sample):
        # find quotients
        quot = int_norm[i, :] / mead_spec
        mead = np.median(quot)
        norm[i, :] = int_norm[i, :] / mead

    return norm

def normalize_2pt(x, y, xmin, xmax, ymin=0, ymax=1):
    """
    Normalize a spectrum to two points, not very rigorous
    """
    xmini = find_nearest_index(x, xmin)
    xmaxi = find_nearest_index(x, xmax)
    y_floor = y - y[xmini]
    y_norm = y_floor/y_floor[xmaxi]
    return x, y_norm