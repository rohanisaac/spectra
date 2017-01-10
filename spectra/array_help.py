"""
NOTE: Not renamed array due to potential problems seen opening notebooks (only affects when started in the same directory, not an issue in general
"""
from __future__ import division
import numpy as np

def find_nearest(array, value):
    idx = (np.abs(array-value)).argmin()
    return array[idx]


def find_nearest_index(array, value):
    """Returns the position(s) of the point(s) in the list closes to the
    passed point"""
    return (np.abs(array-value)).argmin()


# copy range from x,y data based on x values
def copy_range(x, y, xmin, xmax):
    r1 = np.argmin(abs(x - xmin))
    r2 = np.argmin(abs(x - xmax))
    # print r1,r2
    if r1 < r2:
        return x[r1:r2], y[r1:r2]
    elif r1 > r2:
        return x[r2:r1], y[r2:r1]
    else:
        print "Error, no subrange"
