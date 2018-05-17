"""
Calibrate spectrum with neon data
"""
from .peaks import find_peaks
from .fitting import fit_data, line_fit
from .array_help import find_nearest_tolerance
import numpy as np
import os


def calibrate(x, found_peaks, source="neon", tolerance=1):
    """
    Uses a know set of peaks to calibrate a spectrum

    Parameters
    ----------
    x (float array)
        x-data from a particular source (should be in nm) to be calibrated
    found_peaks (float array)
        peaks that have been found in the y data corresponding to the uncalibrated x data
    source (string)
        file name to get the matching data from
    tolerance (float)
        tolerance when searching for matching peaks
    peak_width (float)
        peak width when auto searching for peaks

    Returns
    -------
    calibrated-x (float array)
        x data that has been calibrated
    """

    source_path = os.path.join(os.path.dirname(os.path.realpath(__file__)), "atomic_lines", "%s.txt" % source)
    source_peaks = np.genfromtxt(source_path, delimiter='\t')

    fpeaks = []
    speaks = []
    for f in found_peaks:
        nearest = find_nearest_tolerance(f, source_peaks, tolerance=tolerance)
        if nearest is not None:
            fpeaks.append(f)
            speaks.append(nearest)

    if len(set(speaks)) < len(speaks): # duplicates in data
        print("Duplicate matches, reduce tolerance")
    
    slope, inter, _ = line_fit(np.array(fpeaks), np.array(speaks), errors=False)
    print("corrected: %sx + %s" % (slope, inter))
    return ((x*slope) + inter)