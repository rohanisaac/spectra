spectra
=======

Basic spectral analysis using peak-o-mat (http://lorentz.sourceforge.net/). Designed for optical spectroscopic data such as UV-Vis, IR, Raman, Photoluminescence, but can be used for a variety of other data types such as X-ray diffraction, NMR, Mass spec data.

Prerequisites
-------------

**peak-o-mat (http://lorentz.sourceforge.net/):** Set path in spectra/spectra.py

Basic Usage 
-----------

	# load data
	import spectra
	s_obj = spectra.Specta("/this/is/the/path.txt")

Class functions
---------------

find_background()
	Attempt to find and fit a polynomial function to the data set
	Sets `bg` numpy array
subtract_background()
	Subtracts `bg` from active data
remove_spikes()
	Attempts to very sharp noise from data
guess_peak_width(max_width=50)
	attempt to guess an approximate scale for the width of the peaks to improve peak finding
find_peaks(limit=30, lower=2, upper=10)
	search for peaks in the data set using continuous wavelet transformation 
build_model(peaktype)
	build a model function to fit the data, using the peak profile specified, defaults to lorentzian. Builds a model string used by peak-o-mat
fit_data()
	Uses the peak-o-mat module to perform the fitting routine using the model on the data
output_results()
	Output the fit results from peak-o-mat object

Class data
----------

origin
	A spec (peak-o-mat) format object to store the original x,y data, load the file, as well as to use for the fitting routines
The following y-data (the x-data remains constant for the other operations)
bg
	attempt poly background fit
o
	Original Data
active
	Active Data
bg
	Background
numpeaks
	Number of peaks found
m
	peak-o-mat model with model string `model_str`


Detailed Usage
--------------

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

