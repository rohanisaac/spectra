"""
Driving the functionalized version of spectra
"""
from __future__ import division

import matplotlib.pyplot as plt
import definitions as df

#from scipy import signal, interpolate, fftpack
#import numpy as np
import sys, os

# peak-o-mat stuff
PATH = os.getcwd()+'/'
sys.path.append(PATH + '../peak-o-mat-1.1.9/')

from peak_o_mat.spec import Spec
# from peak_o_mat.fit import Fit
# from peak_o_mat.model import Model

# use peak-o-mat to import model
spec_obj = Spec(PATH+'samples/PMDA-KBr-1-75x1.CSV') 

# use local variables for data
x = spec_obj.x
y = spec_obj.y

# Plot data and bg 
plt.figure(1)
plt.plot(x,y,'-')

bg = df.find_background(x,y)
plt.plot(x,bg,'-r')
    