import numpy as np
from scipy import constants as sc

def wl2wn(wl):
    """
    Converts wavelength (nm) to wavenumber (1/cm)
    """
    return (10.0**7) / wl


def wn2wl(wn):
    """
    Converts wavenumber (1/cm) to wavelength (nm)
    """
    return (10.0**7) / wn

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
    return (2 - np.log10(transmission))

def nm2ev(x, y):
    """
    Properly convert data from wavelength on the x-scale from nm to eV
    
    DOI:10.1021/jz401508t
    """
    c1 = sc.h * sc.c / ( sc.eV * sc.nano) # hc = 1240 eV nm
    x_c = c1 / x # E = hc/lambda in eV
    y_c = y / np.power(x_c, 2) # f(E) = f(lambda) hc/E^2
    return x_c, y_c