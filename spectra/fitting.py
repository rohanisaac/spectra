"""
Various fitting functions, mainly based on lmfit
"""
import numpy as np
from lmfit.models import PolynomialModel
import statsmodels.api as sm
from lmfit import Model
from lmfit.models import ConstantModel, LorentzianModel, GaussianModel, PseudoVoigtModel
from .peaks import (
    gaussian,
    lorentzian,
    voigt,
    guess_peak_width,
    lorentzian_d,
    lorentzian_dd,
    gaussian_d,
    gaussian_dd,
)
from .array_help import copy_range_array, find_nearest_index, copy_range


def fit_peaks(
    x,
    y,
    peak_pos,
    bg="constant",
    sigma_guess=2,
    center_pm=20,
    sigma_min=0.5,
    amplitude_max_m=3.0,
    bg_pm=100,
):
    mod = ConstantModel()

    for i, p in enumerate(peak_pos):
        mod += LorentzianModel(prefix="p%s_" % i)

    pars = mod.make_params()

    for i, p in enumerate(peak_pos):
        pars["p%s_center" % i].set(p, min=p - center_pm, max=p + center_pm)
        pars["p%s_sigma" % i].set(sigma_guess, min=sigma_min)
        # pars['p%s_amplitude' % i].set(10**2, min=0.0)
        pars["p%s_amplitude" % i].set(
            amplitude_max_m * y[find_nearest_index(x, p)], min=0.0
        )

    pars["c"].set(0, min=-1 * bg_pm, max=bg_pm)

    out = mod.fit(y, pars, x=x, method="leastsq")
    out.peak_pos = peak_pos
    return out


def make_region_array(split_pos):
    """
    Make a set of regions based on split positions (including endpoints)
    """
    regions = []
    init = split_pos[0]

    for i in split_pos[1:]:
        regions.append([init, i])
        init = i

    return regions


def make_regions(peaks, n_groups=5, plot=False, x_low=0, x_high=2000):
    """
    Groups peak in the peak list by k-means clustering and thus split the fitting region
    """

    # find the groups by k means (make deterministic/repeatable)
    from sklearn.cluster import KMeans
    from collections import OrderedDict

    X = np.array(peaks).reshape(-1, 1)
    kmeans = KMeans(n_clusters=n_groups, random_state=0).fit(X)
    groups = kmeans.predict(X)

    peaks_pos_list = []
    for i in range(n_groups):
        peaks_pos_list.append(np.where(groups == i)[0])

    # sort the groups by min_peak
    min_pos = []
    for pl in peaks_pos_list:
        min_pos.append(pl[0])

    peak_pos_sor = [x for _, x in sorted(zip(min_pos, peaks_pos_list))]

    # make the divisions as the average between the max of the the group and the min of the next group
    seperators = []
    for i in range(n_groups - 1):
        seperators.append(
            (peaks[peak_pos_sor[i][-1]] + peaks[peak_pos_sor[i + 1][0]]) / 2
        )

    ns = [x_low] + seperators + [x_high]

    if plot:
        prop_cycle = plt.rcParams["axes.prop_cycle"]
        colors = prop_cycle.by_key()["color"]

        fig, ax = plt.subplots(figsize=(15, 3))
        # color code by peaks
        for i, c in zip(peaks, groups):
            ax.axvline(i, c=colors[c], lw=2, label="group %i" % c)
        # label centers
        for i in kmeans.cluster_centers_:
            ax.axvline(i, c="k", ls="--", lw=1, alpha=1)

        # label separators
        for i in seperators:
            ax.axvline(i, c="k", alpha=0.4, lw=5)

        # ax.legend()

        handles, labels = ax.get_legend_handles_labels()
        by_label = OrderedDict(zip(labels, handles))
        ax.legend(by_label.values(), by_label.keys(), ncol=5, loc="upper center")

    return make_region_array(ns)


""" def works
def fit_peaks(x, y, peak_pos):
    mod = ConstantModel()

    for i, p in enumerate(peak_pos):
        mod += LorentzianModel(prefix='p%s_' % i)
    
    pars = mod.make_params()

    for i, p in enumerate(peak_pos):
        pars['p%s_center' % i].set(p, min=p-20, max=p+20)
        pars['p%s_sigma' % i].set(2, min=0.5)
        #pars['p%s_amplitude' % i].set(10**2, min=0.0)
        pars['p%s_amplitude' % i].set(3*y[sp.find_nearest_index(x, p)], min=0.0)

    pars['c'].set(0, min=-100, max=100)
    
    out = mod.fit(y, pars, x=x, method='leastsq')
    out.peak_pos = peak_pos
    return out
"""


def split_and_fit(x, y, sub_ranges, all_peaks, index=False):
    """
    Split data into regions and fit subsets
    """
    out_list = []
    for xmin, xmax in sub_ranges:
        print("Fitting between range {0} and {1}".format(xmin, xmax))
        xt, yt = copy_range(x, y, xmin, xmax)
        if index:
            sub_set = [x[p] for p in all_peaks if x[p] > xmin and x[p] < xmax]
        else:
            sub_set = [p for p in all_peaks if p > xmin and p < xmax]
        print("Peaks in this range are at: {0}".format(sub_set))
        if len(sub_set) == 0:
            print("No peaks to fit, skipping")
            continue
        out = fit_peaks(xt, yt, sub_set)
        out_list.append(out)
    return out_list


def fit_data(x, y, peak_pos, peak_type="LO", width=None):
    """ 
    Builds a lmfit model of peaks in listed by index in `peak_pos` with lineshape
    given by `peak_type`. 

    Parameters
    ----------
    x : array
    y : array
    peak_pos : array of peak postion indices
    
    peak_type : string
        Peaks can be of the following types:
        - 'LO' : symmetric lorentzian (default)
        - 'GA' : symmetric gaussain
        - 'VO' : symmetric pseudo voigt

    width : float
    amp : float

    Returns
    -------
    out : Fitted lmfit model

    """

    # get appropriate line shape function

    if peak_type == "LO":
        peak_function = lorentzian
    elif peak_type == "GA":
        peak_function = gaussian
    elif peak_type == "VO":
        peak_function = voigt

    # build model
    model = Model(peak_function, prefix="p0_")
    for i, p in enumerate(peak_pos[1:]):
        model += Model(peak_function, prefix="p%s_" % str(i + 1))

    pars = model.make_params()

    # initialize params
    for i, p in enumerate(peak_pos):
        pars["p%s_x0" % i].set(x[p])
        pars["p%s_fwhm" % i].set(width, min=0)
        pars["p%s_amp" % i].set(y[p], min=0)

    out = model.fit(y, pars, x=x)
    return out


def fit_data_bg(x, y, peak_pos, peak_type="LO", width=None, bg_ord=0):
    """ 
    Builds a lmfit model of peaks in listed by index in `peak_pos`

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
    out: 
        fitted model

    """
    # need to define peak width finding
    if width is None:
        width = guess_peak_width(x, y)

    # start with polynomial background
    model = PolynomialModel(bg_ord, prefix="bg_")
    pars = model.make_params()

    if peak_type == "LO":
        peak_function = lorentzian
    elif peak_type == "GA":
        peak_function = gaussian
    elif peak_type == "VO":
        peak_function = voigt

    # add peak type for all peaks
    for i, peak in enumerate(peak_pos):
        temp_model = Model(peak_function, prefix="p%s_" % i)
        pars.update(temp_model.make_params())
        model += temp_model

    # set initial background as flat line at zeros
    for i in range(bg_ord + 1):
        pars["bg_c%i" % i].set(0)

    # give values for other peaks, keeping width and height positive
    for i, peak in enumerate(peak_pos):
        pars["p%s_x0" % i].set(x[peak])
        pars["p%s_fwhm" % i].set(width, min=0)
        pars["p%s_amp" % i].set(y[peak], min=0)

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
    for name in list(params.keys()):
        par = params[name]
        dat_out.append((name, par.value, par.stderr))

    output_np = np.array(dat_out)

    if filename or pandas:
        try:
            import pandas as pd

            out_df = pd.DataFrame(
                output_np, columns=["parameter", "value", "stderr"]
            ).set_index("parameter")
            if filename:
                out_df.to_csv(filename)
                return
            else:
                return out_df
        except ImportError:
            print("Need pandas to output to file or pandas dataframe")
    else:
        return output_np


def fit_peak_table(out_list):
    outf = "num,center,c_err,height,h_err,fwhm,f_err\n"
    peak_num = 1
    for o in out_list:
        for i in range(len(o.peak_pos)):
            outf += str(peak_num) + ","
            outf += str(o.params["p%s_center" % i].value) + ","
            outf += str(o.params["p%s_center" % i].stderr) + ","
            outf += str(o.params["p%s_height" % i].value) + ","
            outf += str(o.params["p%s_height" % i].stderr) + ","
            outf += str(o.params["p%s_fwhm" % i].value) + ","
            outf += str(o.params["p%s_fwhm" % i].stderr) + "\n"
            peak_num += 1
    return outf


def batch_fit_single_peak(x, y, xlow, xhigh):
    """
    Batch fit a single peak
    """
    r_xfull, r_y = copy_range_array(x, y, xlow, xhigh)
    r_x = r_xfull[0, :]

    cen_l = []
    amp_l = []
    fwh_l = []

    for i, y in enumerate(r_y):
        peak_pos = [np.argmax(y)]
        fit_out = fit_data_bg(r_x, y, peak_pos)
        data_out = output_results(fit_out)
        cen_l.append(data_out[1][1:])
        amp_l.append(data_out[2][1:])
        fwh_l.append(data_out[3][1:])

    return (
        np.array(cen_l, dtype=np.float64),
        np.array(amp_l, dtype=np.float64),
        np.array(fwh_l, dtype=np.float64),
    )


def line_fit(x, y, errors=True):
    """
    Simple helper function that speeds up line fitting. Uses lmfit
    
    Parameters
    ----------
    x, y (float)
        Sample length x, y arrays of data to be fitted

    Returns
    -------
    Returns slope and intercept, with uncertainties (use uncertainties package if availabe)
    Also returns fit object out, can be dropped
    """
    from lmfit.models import LinearModel

    mod = LinearModel()
    par = mod.guess(y, x=x)
    out = mod.fit(y, par, x=x)

    s = out.params["slope"]
    i = out.params["intercept"]

    if errors:
        try:
            from uncertainties import ufloat

            return ufloat(s.value, s.stderr), ufloat(i.value, i.stderr), out
        except:
            return s.value, s.stderr, i.value, i.stderr, out
    else:
        return s.value, i.value, out


def exponential_fit_offset(
    x, y, amp_guess=1, decay_guess=1, offset_guess=0, errors=True
):
    """
    Simple helper function that speeds up single exponetial fit with offset. Uses lmfit
    
    Parameters
    ----------
    x, y (float)
        Sample length x, y arrays of data to be fitted

    Returns
    -------
    Returns slope and intercept, with uncertainties (use uncertainties package if availabe)

    !!!!!
    not tested or working
    !!!!!
    """
    from lmfit.models import ExpressionModel

    mod = ExpressionModel("offset + amp * exp(-x/decay)")
    par = mod.make_params(amp=amp_guess, decay=decay_guess, offset_guess=offset_guess)
    out = mod.fit(y, params=par, x=x)

    s = out.params["slope"]
    i = out.params["intercept"]

    if errors:
        try:
            from uncertainties import ufloat

            return ufloat(s.value, s.stderr), ufloat(i.value, i.stderr)
        except:
            return s.value, s.stderr, i.value, i.stderr
    else:
        return s.value, i.value


def exponential_fit(x, y, errors=True):
    from lmfit.models import ExponentialModel

    mod = ExponentialModel()
    par = mod.guess(data=y, x=x)
    out = mod.fit(y, params=par, x=x)

    a = out.params["amplitude"]
    d = out.params["decay"]

    return a, d


def poly_fit(x, y, order=1, err=None):
    """
    Use the statsmodel api to perform a OLS fit of the data to a polynomial of order (defualt=1, linear)

    If error are given, perform a WLS instead. Return the fitted model 
    """
    X = np.column_stack((np.power(x, i) for i in range(1, order + 1)))
    X = sm.add_constant(X)
    if err is None:
        mod = sm.OLS(y, X)
    else:
        mod = sm.WLS(y, X, weights=1.0 / (err ** 2))
    return mod.fit()


def build_model(bg_ord=1, peak_num=1, peak_type="lorentzian"):
    model = PolynomialModel(bg_ord, prefix="bg_")

    peak_str = peak_type.lower()[:2]

    if peak_str == "lo":
        peak_function = lorentzian
    elif peak_str == "ga":
        peak_function = gaussian
    elif peak_str == "vo" or peak_type == "ps":
        peak_function = voigt
    else:
        print(
            "'peak_type' should be one of 'lorentzian', 'gaussian' or 'pseudo-voight'"
        )
        return

    for i in range(peak_num):
        model += Model(peak_function, prefix="p%s_" % i)

    params = model.make_params()
    return model, params


def set_parameters(
    par, peak_centers, center_s=np.inf, fwhm=1, fwhm_max=np.inf, amp=1, amp_max=np.inf
):
    for p in par:
        if p.startswith("bg"):
            par[p].set(0e0, min=-1e9, max=1e9)
        if p.endswith("fwhm"):
            par[p].set(fwhm, min=0.1, max=fwhm_max)
        if p.endswith("amp"):
            par[p].set(amp, min=0.1, max=amp_max)

    peak_center_names = sorted([p for p in par if p.endswith("x0")])
    if len(peak_center_names) != len(peak_centers):
        print("Need centers for each peak")

    for pc, pcn in zip(peak_centers, peak_center_names):
        par[pcn].set(pc, min=(pc - center_s), max=(pc + center_s))

    return par


def build_model_d(bg_ord=1, peak_num=1, peak_type="lorentzian"):
    model = PolynomialModel(bg_ord, prefix="bg_")

    peak_str = peak_type.lower()[:2]

    if peak_str == "lo":
        peak_function = lorentzian_d
    elif peak_str == "ga":
        peak_function = gaussian_d
    elif peak_str == "vo" or peak_type == "ps":
        peak_function = voigt
    else:
        print(
            "'peak_type' should be one of 'lorentzian', 'gaussian' or 'pseudo-voight'"
        )
        return

    for i in range(peak_num):
        model += Model(peak_function, prefix="p%s_" % i)

    params = model.make_params()
    return model, params


def build_model_dd(bg_ord=1, peak_num=1, peak_type="lorentzian"):
    model = PolynomialModel(bg_ord, prefix="bg_")

    peak_str = peak_type.lower()[:2]

    if peak_str == "lo":
        peak_function = lorentzian_dd
    elif peak_str == "ga":
        peak_function = gaussian_dd
    elif peak_str == "vo" or peak_type == "ps":
        peak_function = voigt
    else:
        print(
            "'peak_type' should be one of 'lorentzian', 'gaussian' or 'pseudo-voight'"
        )
        return

    for i in range(peak_num):
        model += Model(peak_function, prefix="p%s_" % i)

    params = model.make_params()
    return model, params


def peak_table(out_list, pandas=False):
    """Loop over all peaks fitted over a set of models and print all the information"""
    pn = 0
    table = "center_v,center_std,height_v,height_std,fwhm_v,fwhm_std\n"
    for o in out_list:
        psn = 0
        p = o.params
        while "p%s_center" % psn in p:
            c = p["p%s_center" % psn]
            h = p["p%s_height" % psn]
            f = p["p%s_fwhm" % psn]
            table += "{},{},{},{},{},{}\n".format(
                c.value, c.stderr, h.value, h.stderr, f.value, f.stderr
            )
            psn += 1
    if pandas:
        from io import StringIO
        import pandas as pd

        return pd.read_csv(StringIO(table))
    else:
        return table
