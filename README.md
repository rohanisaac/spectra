spectra
=======

Basic spectral analysis using lmfit. Designed for optical spectroscopic data such as UV-Vis, IR, Raman, Photoluminescence, but can be used for a variety of other data types such as X-ray diffraction, NMR, Mass spec data.

Performs some basic data cleanup (background fitting, noise removal). Can conduct basic peak finding using continuous wavelet transformation and generate a list of parameters (peak position, FWHM, amplitude) to send to `lmfit` as initial model parameters. Fit results including uncertainties are retrieved from `lmfit` and output.

Can plot various stages of the process including raw data, processed data, background estimate, fit model.

**Update: moved from peak-o-mat to lmfit**

**NOTE: Module does not actively run**

Prerequisites
-------------

1. numpy : basic data storage
2. scipy : signal processing
3. lmfit : fitting routines
4. matplotlib : for plotting


Running
-------

Importing module

```python
import sys
sys.path.append('spectra/') # if not on path
import spectra
```

Loading data

```python
S = spectra.Specta("/this/is/the/path.txt")
```

Plotting

```python
import matplotlib.pyplot as plt
plt.plot(S.ox,S.oy,'r-') # original data
plt.plot(S.x,S.y,'b-') # active data
```

Find peaks

```python
S.find_peaks(lower=7, upper=99, limit=8)
plt.plot(S.x,S.y,'b-',S.x[S.peak_pos],S.y[S.peak_pos],'oy')
```

Build model and fit dat

```python
S.build_model()
S.fit_data()
plt.plot(S.x,S.y,'-',S.x,S.fd,'r-') # plot data and fit
S.output_results()
```

Other functions

```python
S.remove_spikes()
```

Member Functions
----------------
For more details see docstring of individual function

- find_background()
	+ Attempt to find background data by a filtering data
	+ Sets `bg` numpy array
- subtract_background()
	+ Subtracts `bg` from active data
- remove_spikes()
	+ Attempts to very sharp noise from data
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

Untested script
---------------

```python
import matplotlib.pyplot as plt
from spectra import Spectra

S = Spectra('samples/SAMPLE.CSV')
S.find_background()

# Plot data and bg
plt.figure(1)
plt.plot(S.x,S.active,'-')
plt.plot(S.x,S.bg,'-r')

S.subtract_background()
S.remove_spikes()
S.guess_peak_width(max_width=50)
S.find_peaks(limit=30)
S.build_model()

# plot spectra, model, found_peaks
plt.figure(2)
plt.plot(S.x,S.active,'-')
plt.plot(S.x,S.model_data,'r-')
plt.plot(S.x[S.peak_pos],S.active[S.peak_pos],'oy')
S.fit_data()

# plot spectra,fit
plt.figure(3)
plt.plot(S.x,S.active,'g-',alpha=0.3)
plt.plot(S.x,S.model_data,'r-')

S.output_results()
```

Todo
----

1. ~~Background search~~ Better background search
2. Rebinning spectra (is this useful?)
3. Normalizing
4. Dealing with a batch of spectra
5. Better data extraction function
6. Working with different peak shapes
7. ~~Error estimates~~
