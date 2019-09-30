"""
Calibrate spectrum with neon data
"""
from .peaks import find_peaks
from .fitting import fit_data, line_fit, poly_fit, fit_data_bg
from .array_help import find_nearest_tolerance
from .convert import rwn2wl, wl2rwn 
from .normalize import normalize
from .read_files import read_horiba
import peakutils as pu
import numpy as np
import matplotlib.pyplot as plt
import os
import statsmodels.api as sm


def neon_peaks(source='neon'):
    """
    Print peaks
    """
    source_path = os.path.join(os.path.dirname(os.path.realpath(__file__)), "atomic_lines", "%s.txt" % source)
    source_peaks = np.genfromtxt(source_path, delimiter='\t')
    return source_peaks


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
            print(f, nearest)
            fpeaks.append(f)
            speaks.append(nearest)

    if len(set(speaks)) < len(speaks): # duplicates in data
        print("Duplicate matches, reduce tolerance")
    
    slope, inter, _ = line_fit(np.array(fpeaks), np.array(speaks), errors=False)
    print("corrected: %sx + %s" % (slope, inter))
    return ((x*slope) + inter)

def find_laser_wavelength(neon_file, laser_file, tolerance=1, offset=0, neon_thres=0.1, plots=False, plot_error_scale=1000):
    """
    1. Load in both files
    2. Find peaks in Neon (list, threshold > 10, min dist > 20)
    3. Find peak in laser (single, >50%)
    4. Fit peaks in Neon (lmfit)
    5. Fit peaks in Laser (lmfit)
    6. Find expected neon peaks that are closest to the measured peaks
    7. Fit (OLS) the expected vs the measured (with errors)
    8. Convert the Measured Laser peak (with error) with (slope/intercept + errors) to True peak
    
    """
    # 1
    print("= Loading peaks =")
    neon = read_horiba(neon_file, x='nm')
    laser = read_horiba(laser_file, x='nm')
    
    # 2, 3
    neon_peaks = pu.indexes(neon['Intensity'], thres=neon_thres, min_dist=20)
    laser_peak = pu.indexes(laser['Intensity'], thres=0.5, min_dist=20)
    print(neon_peaks, laser_peak)
    
    # Error checks
    if len(neon_peaks) < 3:
        print("May be too few peaks for current system. Check data/settings")
    if len(laser_peak) > 1:
        print("Too many laser peaks found. Check data/settings.")
        
    # 4
    print("\n= Fitting neon peaks =")
    # Fit neon data to this model
    outn = fit_data_bg(neon['Wavelength_nm'], neon['Intensity'], neon_peaks, width=1.0, bg_ord=0)

    meas_neon_peaks = []
    meas_neon_peaks_err = []
    for key  in outn.params:
        if key.endswith('x0'):
            val = outn.params[key]
            print(val)
            meas_neon_peaks.append(val.value)
            meas_neon_peaks_err.append(val.stderr)
            
    # 5
    print("\n= Fitting laser peaks =")
    # Fit laser data to this model
    outl = fit_data_bg(laser['Wavelength_nm'], laser['Intensity'], laser_peak, bg_ord=0, width=1.0)

    for key in outl.params:
        if key.endswith('x0'):
            val = outl.params[key]
            print(val)
            meas_laser_peak = val.value
            meas_laser_peak_err = val.stderr
            
    # 6
    print("\n= Finding matching peaks =")
    source = 'neon'

    source_path = os.path.join(os.path.dirname(os.path.realpath(__file__)), "atomic_lines", "%s.txt" % source)
    ref_neon_peaks = np.genfromtxt(source_path, delimiter='\t')

    # usable peaks
    mes_peaks = []
    mes_peaks_err = []
    ref_peaks = []

    # loop through measured peaks and try to find a match
    for p, e in zip(meas_neon_peaks, meas_neon_peaks_err):
        # include offset
        nearest = find_nearest_tolerance(p-offset, ref_neon_peaks, tolerance=tolerance)
        if nearest is not None:
            print(p, nearest)
            mes_peaks.append(p)
            mes_peaks_err.append(e)
            ref_peaks.append(nearest)

    if len(set(ref_peaks)) < len(ref_peaks): # duplicates in data
        print("Duplicate matches, reduce tolerance or add more neon references")
        
    # 7
    print("\n= Fitting data =")
    # Fit with the predicting the reference with the measured
    f2 = poly_fit(np.array(mes_peaks), np.array(ref_peaks), order=1)
    print(f2.summary())
    
    # 8
    print("\n= Making prediction =")
    # Make a prediction based on data

    # put into same form as input
    laser_meas = np.array([1, meas_laser_peak])
    laser_pred = f2.get_prediction(exog=laser_meas, weights=1./(meas_laser_peak_err**2))
    laser_pred_center = laser_pred.predicted_mean[0]
    
    # find the relative error
    # defauts to 95% confidence interval (two-sigma, 0.05); instead using 68% (one-sigma, 0.32)
    laser_pred_error = np.max(np.abs(laser_pred.conf_int(alpha=0.05) - laser_pred_center))
    
    
    # plots
    if plots:
        fig, ax = plt.subplots(1, 2, figsize=(12, 4))
        
        ax[0].plot(neon['Wavelength_nm'], normalize(neon['Intensity']))
        ax[0].plot(laser['Wavelength_nm'], normalize(laser['Intensity']))
        ax[0].plot(neon['Wavelength_nm'][neon_peaks], normalize(neon['Intensity'])[neon_peaks], 'rx', label='Neon peaks')
        ax[0].plot(laser['Wavelength_nm'][laser_peak], normalize(laser['Intensity'])[laser_peak], 'b+', label='Laser peak')
        ax[0].set_xlabel('Wavelength (nm)')
        ax[0].legend()
        
        # scale up errors to see them
        x = np.array(mes_peaks)
        xe = np.array(mes_peaks_err)
        yf = np.array(ref_peaks)
        

        ax[1].errorbar(x, yf, xerr=plot_error_scale*xe, fmt="rx", label='Neon peaks')
        ax[1].plot(x, f2.fittedvalues, 'b', label='Data')
        ax[1].errorbar(meas_laser_peak, laser_pred_center, yerr=laser_pred_error*plot_error_scale, fmt='k+', label='Laser peak estimate')
        ax[1].set_xlabel('Measured values (nm)')
        ax[1].set_ylabel('True values (nm)')
        ax[1].legend()
    
    print("%4.4f ± %4.4f" % (laser_pred_center, laser_pred_error))
    return laser_pred_center, laser_pred_error

def find_best_offset(neon_file, tolerance=20, neon_thres=0.1, offset_range=30, plot=False):
    """
    1. Load in both files
    2. Find peaks in Neon (list, threshold > 10, min dist > 20)
    3. Find peak in laser (single, >50%)
    4. Fit peaks in Neon (lmfit)
    5. Fit peaks in Laser (lmfit)
    6. Find expected neon peaks that are closest to the measured peaks
    7. Fit (OLS) the expected vs the measured (with errors)
    8. Convert the Measured Laser peak (with error) with (slope/intercept + errors) to True peak
    
    """
    # 1
    print("= Loading peaks =")
    neon = read_horiba(neon_file, x='nm')
    
    # 2, 3
    neon_peaks = pu.indexes(neon['Intensity'], thres=neon_thres, min_dist=20)
    print(neon_peaks)
    
    # Error checks
    if len(neon_peaks) < 3:
        print("May be too few peaks for current system. Check data/settings")
        
    # 4
    print("\n= Fitting neon peaks =")
    # Fit neon data to this model
    outn = fit_data_bg(neon['Wavelength_nm'], neon['Intensity'], neon_peaks, width=1.0, bg_ord=0)

    meas_neon_peaks = []
    meas_neon_peaks_err = []
    for key  in outn.params:
        if key.endswith('x0'):
            val = outn.params[key]
            print(val)
            meas_neon_peaks.append(val.value)
            meas_neon_peaks_err.append(val.stderr)
            
    # 6
    print("\n= Searching for best offset value =")
    source = 'neon'

    source_path = os.path.join(os.path.dirname(os.path.realpath(__file__)), "atomic_lines", "%s.txt" % source)
    ref_neon_peaks = np.genfromtxt(source_path, delimiter='\t')

    
    offsets = np.arange(-1*offset_range, offset_range, 0.1)
    totals = []
    #for each offset
    for offset in offsets:
        # loop through measured peaks and try to find a match
        # usable peaks
        total = 0
        for p, e in zip(meas_neon_peaks, meas_neon_peaks_err):
            # include offset
            nearest = find_nearest_tolerance(p-offset, ref_neon_peaks, tolerance=tolerance)
            if nearest is not None:
                total += (nearest-(p-offset))**2
        totals.append(total)
        
    if plot:
        plt.semilogy(offsets, totals)
        plt.xlabel("Offset value")
        plt.ylabel("Total absolute difference (smaller is better)")
    
    best_offset = offsets[np.argmin(np.array(totals))]
    print("%0.1f" % best_offset)
    return best_offset


def find_best_offset2(neon_file, laser_wl=633, source='neon', tolerance=20, neon_thres=0.1, offset_range=30, plot=False):
    """
    1. Load in Neon file, convert wavenumber (1/cm) to wavelength (nm)
    2. Find peaks in Neon (list, threshold > 10, min dist > 20)
    3. Fit peaks in Neon (lmfit)
    6. Find expected neon peaks that are closest to the measured peaks
    7. Fit (OLS) the expected vs the measured (with errors)
    8. Convert the Measured Laser peak (with error) with (slope/intercept + errors) to True peak
    
    """
    # 1
    print("= Loading peaks =")
    neon = read_horiba(neon_file)
    neon_y = neon['Intensity']
    neon_x = neon['Relative_Wavenumber']
    neon_x_nm = rwn2wl(neon_x, laser_wl)
    
    # 2
    neon_peaks = pu.indexes(neon_y, thres=neon_thres, min_dist=20)
    print(neon_peaks)
    
    # Error checks
    if len(neon_peaks) < 3:
        print("May be too few peaks for current system. Check data/settings")
        
    # 4
    print("\n= Fitting neon peaks =")
    # Fit neon data to this model
    outn = fit_data_bg(neon_x_nm, neon_y, neon_peaks, width=1.0, bg_ord=0)

    meas_neon_peaks = []
    meas_neon_peaks_err = []
    for key  in outn.params:
        if key.endswith('x0'):
            val = outn.params[key]
            print(val)
            meas_neon_peaks.append(val.value)
            meas_neon_peaks_err.append(val.stderr)
            
    # 6
    print("\n= Searching for best offset value =")
    source_path = os.path.join(os.path.dirname(os.path.realpath(__file__)), "atomic_lines", "%s.txt" % source)
    ref_neon_peaks = np.genfromtxt(source_path, delimiter='\t')

    offsets = np.arange(-1*offset_range, offset_range, 0.1)
    totals = []
    #for each offset
    for offset in offsets:
        # loop through measured peaks and try to find a match
        # usable peaks
        total = 0
        for p, e in zip(meas_neon_peaks, meas_neon_peaks_err):
            # include offset
            nearest = find_nearest_tolerance(p-offset, ref_neon_peaks, tolerance=tolerance)
            if nearest is not None:
                total += (nearest-(p-offset))**2
        totals.append(total)
        
    if plot:
        plt.semilogy(offsets, totals)
        plt.xlabel("Offset value")
        plt.ylabel("Total absolute difference (smaller is better)")
    
    best_offset = offsets[np.argmin(np.array(totals))]
    print("%0.1f" % best_offset)
    return best_offset

def calibrate_x_data2(neon_file, laser_wl, tolerance=1, offset=0, neon_thres=0.1, plots=False, plot_error_scale=1000):
    """
    1. Load in neon file (assume data file has the same x-data)
    2. Find peaks in Neon (list, threshold > 10, min dist > 20)
    3. Fit peaks in Neon (lmfit)
    4. Finding matching peaks
    5. Fit (OLS) the expected vs the measured (with errors)
    6. Generating new x-values
    
    """
    # 1
    print("= Loading peaks =")

    neon = read_horiba(neon_file)
    neon_y = neon['Intensity']
    x_old = neon['Relative_Wavenumber']
    x_old_nm = rwn2wl(x_old, laser_wl)
    
    # 2
    neon_peaks = pu.indexes(neon_y, thres=neon_thres, min_dist=20)
    print(neon_peaks)
    
    # Error checks
    if len(neon_peaks) < 3:
        print("May be too few peaks for current system. Check data/settings")
        
    # 3
    print("\n= Fitting neon peaks =")
    # Fit neon data to this model
    outn = fit_data_bg(x_old_nm, neon_y, neon_peaks, width=1.0, bg_ord=0)

    meas_neon_peaks = []
    meas_neon_peaks_err = []
    for key  in outn.params:
        if key.endswith('x0'):
            val = outn.params[key]
            print(val)
            meas_neon_peaks.append(val.value)
            meas_neon_peaks_err.append(val.stderr)
            
    # 4
    print("\n= Finding matching peaks =")
    source = 'neon'

    source_path = os.path.join(os.path.dirname(os.path.realpath(__file__)), "atomic_lines", "%s.txt" % source)
    ref_neon_peaks = np.genfromtxt(source_path, delimiter='\t')

    # usable peaks
    mes_peaks = []
    mes_peaks_err = []
    ref_peaks = []

    # loop through measured peaks and try to find a match
    for p, e in zip(meas_neon_peaks, meas_neon_peaks_err):
        # include offset
        nearest = find_nearest_tolerance(p-offset, ref_neon_peaks, tolerance=tolerance)
        if nearest is not None:
            print(p, nearest)
            mes_peaks.append(p)
            mes_peaks_err.append(e)
            ref_peaks.append(nearest)

    if len(set(ref_peaks)) < len(ref_peaks): # duplicates in data
        print("Duplicate matches, reduce tolerance or add more neon references")
        
    # 5
    print("\n= Fitting data =")
    # Fit with the predicting the reference with the measured
    f2 = poly_fit(np.array(mes_peaks), np.array(ref_peaks), order=1)
    print(f2.summary())
    
    # 6
    print("\n= Generating new x-values =")

    if len(neon_peaks) == 1:
        # No fitting, just an offset
        x_new_nm = x_old_nm + (ref_peaks[0] - mes_peaks[0])
        x_new = wl2rwn(x_new_nm, laser_wl)
    else:
        x_new_nm = f2.get_prediction(exog=sm.add_constant(x_old_nm))
        x_new = wl2rwn(x_new_nm.predicted_mean, laser_wl)


    print("Difference at the start -> Old: %3.3f; New: %3.3f" % (x_old[0], x_new[0]))
    print("Difference at the end   -> Old: %3.3f; New: %3.3f" % (x_old[-1], x_new[-1]))

    return x_new