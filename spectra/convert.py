import numpy as np
from scipy import constants as sc


def wl2wn(wl):
    """
    Converts wavelength (nm) to wavenumber (1/cm)
    """
    return (10.0 ** 7) / wl


def wl2rwn(wl, wlr):
    """
    Converts wavelength (nm) to relative wavenumber (1/cm) with a reference wavelength
    """
    return wl2wn(wlr) - wl2wn(wl)


def wn2wl(wn):
    """
    Converts wavenumber (1/cm) to wavelength (nm)
    """
    return (10.0 ** 7) / wn


def rwn2wn(wn, wl):
    """
    Convert relative wavenumber to absolute using nm shift
    """
    return wl2wn(wl) - wn


def rwn2wl(wn, wl):
    """
    Converts relative wavenumber to wavelength using a wavelength reference
    """
    return wn2wl(rwn2wn(wn, wl))


def absorption(transmission):
    """
    Convert transmission data from a percent scale (0-100) to an absorption values
    """
    return 2 - np.log10(transmission)


def nm2ev(x):
    """
    Convert nm to eV
    """
    c1 = sc.h * sc.c / (sc.eV * sc.nano)  # hc = 1240 eV nm
    return c1 / x


def nm2ev_xy(x, y):
    """
    Properly convert data from wavelength on the x-scale from nm to eV
    
    DOI:10.1021/jz401508t
    """
    c1 = sc.h * sc.c / (sc.eV * sc.nano)  # hc = 1240 eV nm
    x_c = c1 / x  # E = hc/lambda in eV
    y_c = y / np.power(x_c, 2)  # f(E) = f(lambda) hc/E^2
    return x_c, y_c


def nm2ev_xyz(x, y, z):
    """
    Properly convert 2d data from nm to eV
    
    Parameters
    ----------
    
    x: some axis
    y: nm
    z: intensity
    
    Returns
    -------
    x: same
    y: eV
    z: corrected intensity
    """
    print(x.shape, y.shape, z.shape)
    c1 = sc.h * sc.c / (sc.eV * sc.nano)  # hc = 1240 eV nm
    y_c = c1 / y  # E = hc/lambda in eV
    z_c = z.T / np.power(y_c, 2)  # f(E) = f(lambda) hc/E^2
    return x, y_c, z_c.T
