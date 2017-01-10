def wl2wn(wl):
    """Converts wavelength (nm) to wavenumber (1/cm)"""
    return (10.0**7) / wl


def wn2wl(wn):
    """Converts wavenumber (1/cm) to wavelength (nm)"""
    return (10.0**7) / wn
