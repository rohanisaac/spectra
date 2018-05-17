"""
NOTE: filename not named array.py
    Potential problems with opening notebooks, array.py is a file that is already in use
    (Not a general issue, only affects when started in the same directory)
"""

import numpy as np

def find_nearest(array, value):
    """
    Returns the position(s) of the point(s) in the list closes to the
    passed point
    """
    idx = (np.abs(array-value)).argmin()
    return array[idx]


def find_nearest_index(array, value):
    """
    Returns the position(s) of the point(s) in the list closes to the
    passed point
    """
    return (np.abs(array-value)).argmin()


# copy range from x,y data based on x values
def copy_range(x, y, xmin, xmax):
    """
    Copy a range of x, y values based on x values
    """
    r1 = np.argmin(abs(x - xmin))
    r2 = np.argmin(abs(x - xmax))
    # print r1,r2
    if r1 < r2:
        return x[r1:r2], y[r1:r2]
    elif r1 > r2:
        return x[r2:r1], y[r2:r1]
    else:
        print("Error, no subrange")

def copy_range_array(x, y, xmin, xmax):
    """
    Copy range from x,y data based on x values, but here x and y are both arrays
    """
    x0 = x[0, :]
    r1 = np.argmin(abs(x0 - xmin))
    r2 = np.argmin(abs(x0 - xmax))
    # print r1,r2
    if r1 < r2:
        return x[:, r1:r2], y[:, r1:r2]
    elif r1 > r2:
        return x[:, r2:r1], y[:, r2:r1]
    else:
        print("Error, no subrange")

def sort_by_list(names, x, y):
    """
    Sort x and y array by the list values

    The first index of the arrays x and y should be the same length as names

    Returns
    -------
    names, x, y
        sorted version of all three arrays
    """
    # note slx and sly are lists
    snames, slx, sly = zip(*sorted(zip(names, list(x), list(y))))
    return snames, np.array(slx), np.array(sly)


def find_nearest_tolerance(value, array, tolerance=10):
    """
    Returns nearest value in the list if it falls within a tolerance,
    none otherwise
    """
    idx = (np.abs(array-value)).argmin()
    if np.abs(array[idx]-value) <= tolerance:
        return array[idx]
    else:
        return None