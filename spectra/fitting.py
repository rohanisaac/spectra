from __future__ import division
import numpy as np
from lmfit.models import PolynomialModel
from lmfit import Model
from peaks import gaussian, lorentzian, voigt, guess_peak_width


def fit_data(x, y, peak_pos, peak_type='LO', max_width=None, bg_ord=2):
    """ Builds a lmfit model of peaks in listed by index in `peak_pos`
    Uses some basic algorithms to determine initial parameters for
    amplitude and fwhm (limit on fwhm to avoid fitting background as peaks)

    Parameters
    ----------
    peak_type : string (default='lorentizian')
        Peaks can be of the following types:

        - 'LO' : symmetric lorentzian
        - 'GA' : symmetric gaussain
        - 'VO' : symmetric pseudo voigt

    max_width : int (default = total points/10)
        max width (in data points) that peak fitted can be

    bg_ord: int
        order of the background polynomial
        0: constant, 1: linear, ...

    Returns
    -------
    pars : model parameters
    model : model object


    """
    # need to define peak width finding
    pw = guess_peak_width(x, y)
    peak_guess = x[peak_pos]

    # start with polynomial background
    model = PolynomialModel(bg_ord, prefix='bg_')
    pars = model.make_params()

    if peak_type == 'LO':
        peak_function = lorentzian
    elif peak_type == 'GA':
        peak_function = gaussian
    elif peak_type == 'VO':
        peak_function = voigt

    # add lorentizian peak for all peaks
    for i, peak in enumerate(peak_guess):
        temp_model = Model(peak_function, prefix='p%s_' % i)
        pars.update(temp_model.make_params())
        model += temp_model

    # set inital background as flat line at zeros
    for i in range(bg_ord + 1):
        pars['bg_c%i' % i].set(0)

    # give values for other peaks
    for i, peak in enumerate(peak_pos):
        # could set bounds #, min=x[peak]-5, max=x[peak]+5)
        pars['p%s_x0' % i].set(x[peak])
        pars['p%s_fwhm' % i].set(pw / 2, min=pw * 0.25, max=pw * 2)
        # here as well #, min=0, max=2*max(y))
        pars['p%s_amp' % i].set(y[peak])

    out = model.fit(y, pars, x=x)
    return out

def output_results(fit_model, filename=None, pandas=False):
    """ Return fit paramters and standard error, modified from lmfit
    class. Can output same data to file if passed file path

    Arguments
    ---------
    filename : string (default: None)
        If filename is given write data to that file as csv

    pandas : bool (default: False)
        True: Return data as pandas dataframe
        False: Return data as string (csv)
    """
    params = fit_model.params
    dat_out = []
    for name in params.keys():
        par = params[name]
        dat_out.append((name, par.value, par.stderr))

    output_np = np.array(dat_out)

    if filename or pandas:
        try:
            import pandas as pd
            out_df = pd.DataFrame(output_np,
                                  columns=['parameter', 'value', 'stderr']
                                  ).set_index('parameter')
            if filename:
                out_df.to_csv(filename)
                return
            else:
                return out_df
        except ImportError:
            print "Need pandas to output to file or pandas dataframe"
    else:
        return output_np

