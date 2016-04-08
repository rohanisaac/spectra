spectra
=======

Basic spectral analysis and automation using lmfit. Designed for optical spectroscopic data such as UV-Vis, IR, Raman, Photoluminescence, but can be used for a variety of other data types such as X-ray diffraction, NMR, Mass spec data.

Performs some basic data cleanup (background fitting, noise removal). Can conduct basic peak finding using continuous wavelet transformation and generate a list of parameters (peak position, FWHM, amplitude) to send to `lmfit` as initial model parameters. Fit results including uncertainties are retrieved from `lmfit` and output.

Can plot various stages of the process including raw data, processed data, background estimate, fit model.

Prerequisites
-------------

1. numpy : basic data storage
2. scipy : signal processing
3. lmfit : fitting routines
4. matplotlib : for plotting


Running
-------

```python
import sys
sys.path.append('spectra/') # if not on path
import spectra

# Loading data
S = spectra.Specta("/this/is/the/path.txt")

# Plotting
import matplotlib.pyplot as plt
plt.plot(S.x,S.y,'r-')

# Finding background
S.reset() # reset y-data to original data (not really needed here)
S.find_background(cutoff=1/4500, order=3)
plt.plot(S.x,S.y,'r-',S.x,S.bg,'b-')

# Subtract background
S.subtract_background()
plt.plot(S.x,S.y,'r-')

# Find peaks
S.find_peaks()
plt.plot(S.x,S.y,'b-',S.x[S.peak_pos],S.y[S.peak_pos],'oy')

# Build model and fit dat
S.build_model()
S.fit_data()

# plot data and fit
plt.plot(S.x,S.y,'-',S.x,S.out.best_fit,'b-')
plt.plot(S.x,S.y,'-',S.x,S.fd,'r-')

# output results
print S.output_results()
```


Member Functions
----------------
For more details see docstring of individual function

- find_background()
	+ Attempt to find background data by a filtering data
	+ Sets `bg` numpy array
- subtract_background()
	+ Subtracts `bg` from active data
- guess_peak_width(max_width=50)
	+ attempt to guess an approximate scale for the width of the peaks to improve peak finding
- find_peaks(limit=30, lower=2, upper=10)
	+ search for peaks in the data set using continuous wavelet transformation
- build_model(peaktype)
	+ build a model function to fit the data, using the peak profile specified, defaults to lorentzian. Builds a model string used by peak-o-mat
- fit_data()
	+ Uses the peak-o-mat module to perform the fitting routine using the model on the data
- output_results()
	+ Output the fit results from peak-o-mat object


Plotting
--------

```python
plt.plot(S.x,S.y,'-') # active data
plt.plot(S.x,S.bg,'-r') # backround
plt.plot(S.x,S.md,'r-') # model data
plt.plot(S.x[S.peak_pos],S.y[S.peak_pos],'oy') # peak positions
plt.plot(S.x,S.md,'r-') # fitted y-data
```


Data Members
------------

x : active x-data
y : active y-data
y_bak : original y-data
x_bak : original x-data
bg : background y-data
y_smooth : smoothed y-data
num_points : length of data
num_peaks : number of peaks found
peak_pos (list) : x positions of peaks
model : lmfit model used in fitting
pars : parameters to optimize with lmfit
out : fitted model
test_peak_width : peak width used for finding peaks

Todo
----

1. ~~Background search~~ Better background search
2. Rebinning spectra (is this useful?)
3. Normalizing
4. Dealing with a batch of spectra
5. Better data extraction function
6. Working with different peak shapes
7. ~~Error estimates~~
