
import numpy as np
from scipy.constants import h, c, k
from .peaks import lorentzian

A1 = (2 * np.pi**2 / 45) * (h / c)
A2 = -1 * h * c / k

def crop(self, xmin, xmax):
    """
    Crops data using x values

    Args:
        xmin: min x value
        ymax: max y value
    """
    r1 = np.argmin(abs(self.x - xmin))
    r2 = np.argmin(abs(self.x - xmax))
    # print r1,r2
    if r1 < r2:
        self.x, self.y = self.x[r1:r2], self.y[r1:r2]
    elif r1 > r2:
        self.x, self.y = self.x[r2:r1], self.y[r2:r1]
    else:
        print("Error, no subrange")


def generate_spectrum(wavelengths, peak_pos, peak_amp, width=5):
    """Generate a simulated spectrum using Lorentzian functions of width over
    the with peak positions and peak amplitudes"""
    y = np.zeros(len(wavelengths))
    y_blank = np.zeros(len(wavelengths))
    for p, a in zip(peak_pos, peak_amp):
        y += lorentzian(y_blank, a, p, width)
    return y


def activity_to_intensity(S, v0, v, T):
    """
    From Polavarapu et al. 1990
    """
    return A1 * (v0 - v)**4 * S / (v * (1 - np.exp(A2 * v / T)))


def find_common(sub_list, full_list, peak_range):
    """ Find all occurances of sub_list positions within full_list within
    a range of peak_range

    Usage
    -----
    res1, res2 = find_common(list1, list2, width)

    """
    com_sub = []
    com_full = []
    for sub in sub_list:
        for full in full_list:
            if sub < (full + peak_range) and sub > (full - peak_range):
                com_sub.append(sub)
                com_full.append(full)
    return com_sub, com_full
