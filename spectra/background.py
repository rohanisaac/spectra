def find_background(self, cutoff=None, order=2):
    """ Attempts to find the background of the spectra using low
    pass filter

    Arguments
    ---------
    cutoff (float) [default: 2/len(y)]
        cutoff frequency at which to filter the data
    order (int) [default: 2]
        filter order

    Returns
    -------
    bg : array
        background spectrum

    """
    print "Finding background ... "
    if cutoff is None:
        cutoff = 2 / len(self.y)
    # use the low pass filter to find the background
    self.bg = self.butter_lp_filter(cutoff, order)
    return self.bg

def subtract_background(self):
    """ Subtract background from active spectra """
    print "Subtracting background ... "
    self.y = self.y - self.bg


