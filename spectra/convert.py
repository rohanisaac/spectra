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
