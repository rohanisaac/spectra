"""
Driving the functionalized version of spectra
"""
from __future__ import division

import matplotlib.pyplot as plt
from definitions import *
from helper_functions import *

#from scipy import signal, interpolate, fftpack
#import numpy as np
import sys, os

# peak-o-mat stuff
PATH = os.getcwd()+'/'
sys.path.append(PATH + '../peak-o-mat-1.1.9/')

from peak_o_mat.spec import Spec
# from peak_o_mat.fit import Fit
# from peak_o_mat.model import Model

S = Spec(PATH+'samples/PMDA-KBr-1-75x1.CSV') 

# Plot data and bg    
plt.figure(1)
plt.plot(S.x,S.y,'-')

bg = find_background(S.x,S.y)
plt.plot(S.x,bg,'-r')
    