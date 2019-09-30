import matplotlib.pyplot as plt
from .file_io import getxy
from lmfit.lineshapes import lorentzian, gaussian, pvoigt
import numpy as np

def plotxy(filename):
    x, y = getxy(filename)
    plt.plot(x, y)
    plt.title(filename)

def fit_plot_single(range):
    pass

def plot_peak_fit(out):
    """
    Plot the results of a peak fitting routine

    Requires: 
    out: lmfit fit output
    """
    fig, ax = plt.subplots(2, 1, sharex=True, figsize=(12, 6), gridspec_kw={'height_ratios': [1, 4]})
    xd = out.userkws[out.model.independent_vars[0]]
    xnew = np.linspace(xd[0], xd[-1], 1000)
    ax[0].plot(xd, out.residual, 'k-', ms=3, lw=1, alpha=0.5)
    ax[1].plot(xd, out.data, 'b+-', ms=3, lw=1, alpha=0.8)
    ax[1].plot(xnew, out.eval(x=xnew), 'k--', lw=1.5)
    for i, p in enumerate(out.peak_pos):
        #print(out.params['p%s_'])
        ax[1].plot(xnew, lorentzian(xnew, \
                                    amplitude = out.params['p%s_amplitude' % i], \
                                    center = out.params['p%s_center' % i], \
                                    sigma = out.params['p%s_sigma' % i]), lw=0.5)
    return fig


def plot_components(out):
    independent_var = out.model.independent_vars[0]
    x = out.userkws[independent_var]
    
    fig, ax = plt.subplots(figsize=(12, 6))
    comp = out.eval_components(x=x)
    
    for mname, mval in comp.items():
        ax.plot(x, mval, label=mname[:-1], lw=1, alpha=0.5)

    fig.subplots_adjust(wspace=0, hspace=0)
        
    ax.scatter(x, out.data, s=2, alpha=0.5, label='data')
    ax.plot(x, out.eval(x=x), ls='--', lw=2, alpha=0.5, c='k', label='fit')
    fig.legend()
    return fig