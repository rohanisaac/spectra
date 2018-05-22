"""
File loaders for common text based file formats. Importing them into a standard numpy array with column headers. 
No other details are saved
"""
import numpy as np
# 
# Helper functions
#

def data_details(f):
    """
    Return column names and data size of array of a function returning a single numpy array
    """
    def wrapped(*args, **kwargs):
        dat = f(*args, **kwargs)
        print(dat.dtype.names, dat.shape)
        return dat
    return wrapped

# Loaders

@data_details
def read_cary(file):
    """
    Read text file created by CARY Scan software and return a single numpy array with the data and
    column labels
    
    Parameters
    ----------
    filename: str
        Full path to the file to open (relative/absolute)
    
    """
    with open(file) as f:
        for i, l in enumerate(f):
            # extract relevant column headers
            if i==0:
                samples = l.split(',')
            elif i==1:
                columns = l.split(',')
            # stop reading when hit a blank line
            elif l in ['\n', '\r\n']:
                break
    ncols = len(columns) - 1
    return np.genfromtxt(file, delimiter=',', skip_header=2, max_rows=(i-2), usecols=range(ncols), names=columns[:ncols])

@data_details
def read_craic(file):
    """
    Read text file created by CRAIC microphotospectrometer software and return a single numpy array with the data and
    column labels
    
    Parameters
    ----------
    filename: str
        Full path to the file to open (relative/absolute)
    
    """
    with open(file) as f:
        for i, l in enumerate(f):
            # extract relevant column headers
            if i==2:
                columns = l.split(',')
            # stop reading when hit a blank line
            elif l in ['\n', '\r\n']:
                break
    return np.genfromtxt(file, delimiter=',', skip_header=9, names=['Wavelength', columns[0]])

@data_details
def read_nicolet(file):
    """
    Read text file created by Nicolet Nexus FTIR software and return a single numpy array with the data and
    column labels
    
    Parameters
    ----------
    filename: str
        Full path to the file to open (relative/absolute)
    
    """
    return np.genfromtxt(file, delimiter=',', names=['Wavenumber', 'Absorbance'])

@data_details
def read_horiba(file):
    """
    Read text file created by Horiba LabSpec software and return a single numpy array with the data and
    column labels
    
    Parameters
    ----------
    filename: str
        Full path to the file to open (relative/absolute)
    
    """
    return np.genfromtxt(file, delimiter='\t', names=['Relative_Wavenumber', 'Intensity'])