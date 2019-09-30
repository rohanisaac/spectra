spectra
=======

Various routines for analysis and preprocessing of spectral data. 

Prerequisites
-------------
1. numpy: basic data storage and manipulation
2. scipy: signal processing
3. lmfit: non-linear fitting
4. matplotlib: plotting
5. pandas: tables


Running
-------
```
# add to $PYTHONPATH or add to system path
import sys
sys.path.append('spectra/')  # relative or absolute path to dir

import spectra as sp

# call appropriate functions
sp.get_xy('')
```


Contents
--------

### array_help
* find_nearest
* find_nearest_index
* find_nearest_tolerance
* copy_range
* copy_range_array
* copy_range_yarray
* sort_by_list
* sort_array_column
* split_columns

### background
*See (https://github.com/charlesll/rampy)[rampy] for better functions

### calibrate
* neon_peaks
* calibrate
* calibrate_x_data2
* find_laser_wavelength
* find_best_offset
* find_best_offset2

### convert

* wl2wn
* wl2rwn
* wn2wl
* rwn2wn
* rwn2wl
* absorption
* nm2ev
* nm2ev_xy
* nm2ev_xyz

### file_io

* path
* rpath
* assure_path_exists
* make_fname
* write2col
* getxy
* clean_file
* load_folder
* quick_load_xy
* list_all_files
* write_json
* read_json

### filters
* smooth_data
* butter_lp_filter
* butter_lowpass
* butter_lowpass_filter
* wicker
* savgol_filter
* butterworth_bandpass
* resample

### fitting
* fit_peaks
* split_and_fit
* fit_data
* fit_data_bg
* output_results
* fit_peak_table
* batch_fit_single_peak
* line_fit
* exponential_fit_offset
* exponential_fit
* poly_fit
* build_model
* set_parameters
* build_model_d
* build_model_dd

### misc
* crop
* generate_spectrum
* activity_to_intensity
* find_common
* remove_absorption_jumps

### normalize
* normalize
* normalize_msc
* normalize_pq
* normalize_2pt
* normalize_fs

### peaks
* find_peaks
* gaussian
* lorentzian
* voigt
* guess_peak_width
* find_fwhm
* lorentzian_d
* lorentzian_dd
* gaussian_d
* gaussian_dd

### plot
* plot_xy
* fit_plot_single
* plot_peak_fit
* plot_components

### read_files
* data_details
* read_cary
* read_craic
* read_nicolet
* read_horiba
* read_renishaw


Todo
----

1. ~~Background search~~ Better background search
2. Rebinning spectra (is this useful?)
3. Normalizing
4. Dealing with a batch of spectra
5. Better data extraction function
6. Working with different peak shapes
7. ~~Error estimates~~
