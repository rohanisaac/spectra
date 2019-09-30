import os
import re
import numpy as np
import json

def path(*args):
    """
    Return path relative to home

    Parameters
    ----------
    Comma separated directories
    """
    return os.path.join(os.path.expanduser('~'), *args)

def rpath(*args):
    """
    Return path relative to current working directory

    Parameters
    ----------
    Comma separated directories
    """
    return os.path.join(os.getcwd(), *args)

def assure_path_exists(path):
    """
    Create directory structure of file if it doens't already exist
    """
    dir = os.path.dirname(path)
    if not os.path.exists(dir):
        os.makedirs(dir)
    return os.path.abspath(path)

    
def make_fname(fol_name, fil_name, extension):
    """
    Make a path from a folder name, filename and extension
    If the filename already has an extension it will be stripped off
    If the folder name doesn't exist, it will be created
    """
    path = os.path.join(os.getcwd(), fol_name, '%s.%s' % (os.path.splitext(fil_name)[0], extension))
    direct = os.path.dirname(path)
    if not os.path.exists(direct):
        os.makedirs(direct)
    return path

def write2col(fname, x, y, delim='\t'):
    """
    Write two column data to file
    """
    if not os.path.exists(os.path.dirname(fname)):
        os.makedirs(os.path.dirname(fname))
    with open(fname, 'w') as f:
        for xi, yi in zip(x, y):
            f.write("%s%s%s\n" % (xi, delim, yi))
    return

def getxy(file_name, headers=False):
    """Extracts x and y data numpy arrays from passed filename.

    Arguments
    ---------
    file_name: string
        full/relative path to file to extract data from

    Returns
    -------
    if headers == True:
        x,y,xlab,ylab (np.array, np.array, str, str)

    else:
        x,y
    where:
        x,y: 1-d numpy arrays containing first and second column of data
        present in file xlab, ylab: header information about columns
    """

    line_pos = 0  # active line number
    start_pos = 0
    end_pos = 0
    dat_read = False  # has data been read
    header = ''
    sep_char = ','  # default to csv

    with open(file_name, 'rb') as fil:
        # find header and footer positions
        for lin in fil:
            line_pos += 1
            # if line contains any of the alphabet (except e for exponents,
            # not data)
            if re.search(b'[a-df-zA-DF-Z]', lin):
                if not dat_read:
                    # before data has been read, set start of data pos
                    start_pos = line_pos
                    header = lin
                if dat_read and end_pos == 0:
                    # after data had been read and before end position has
                    # been set, set end pos
                    end_pos = line_pos

            # if data line and data has not been read
            elif not dat_read:
                # find seperator
                if re.search(b'\t', lin):
                    sep_char = '\t'
                elif re.search(b',', lin):
                    sep_char = ','
                else:
                    print("Unknown separator character")
                # now we know what separator is for the data
                dat_read = True

            # data line and we already know what the separator character is
            else:
                continue

    # if we didn't find an end position
    if end_pos == 0:
        skip_foot = 0
    else:
        # if we did compute it
        skip_foot = line_pos - end_pos

    xlab, ylab = ('', '')
    # find header row if exists
    header_lst = header.split(sep_char)
    # print headerlst
    if len(header_lst) >= 2:
        xlab, ylab = header_lst[0:2]

    # attempt to load into numpy array, see what happens
    fdat = np.genfromtxt(
        file_name, delimiter=sep_char, skip_header=start_pos,
        skip_footer=skip_foot)

    if headers:
        return fdat[:, 0], fdat[:, 1], xlab, ylab
    else:
        return fdat[:, 0], fdat[:, 1]

""" There has to be a regular expression that does the same thing that the
code above does, just need to find it out and test it on a full set of data
"""


def clean_file(infile, outfile):
    """
    Cleans out header information from file and writes it to another file

    Arguments
    ---------
    infile (str)
        Path to input file

    outfile (str)
        Path to output file (folder should be created)
    """
    x, y = getxy(infile)
    write2col(outfile, x, y)

def load_folder(path, extension='csv', reader=getxy, regex='[0-9]+'):
    """
    Load all the files from a folder assuming they are a number of files
    file1.ext
    file2.ext
    ...
    filen.ext

    Parameters
    ----------
    path (str)
        full path to the folder
    extension (str)
        file extension to filter by in the folder
    reader (function)
        Function to read the files with

    Returns
    -------
    y-values (2d array)
        (nfiles x ndata) array containing all the yvalues
    x-values (2d array)
        (nfiles x ndata) array containing all the xvalues
    filenames (list)
        array containing all the filenames corresponding to the each file
        processed in the order they were processed
    """
    # start with arrays b/c we don't know size
    x_values = []
    y_values = []
    filenames = []
    sort_key = []

    for f in os.listdir(path):
        if f.lower().endswith(extension):
            dat = reader(os.path.join(path, f))
            if isinstance(dat, tuple) and len(dat) == 2:
                # If they return two elements, must be x and y
                x, y = dat
            elif isinstance(dat, np.ndarray):
                if dat.dtype.names is not None:
                    # if array has header
                    header = dat.dtype.names
                    x = dat[header[0]]
                    y = dat[header[1]]
                else:
                    x = dat[:, 0]
                    y = dat[:, 1]
            else:
                print("Cannot load file %s" % f)
                continue
            x_values.append(x)
            y_values.append(y)
            filenames.append(f)
            sort_key.append(int(re.search(regex, f).group(0)))
    sk, sx, sy, sf = zip(*sorted(zip(sort_key, x_values, y_values, filenames)))

    return np.array(sx), np.array(sy), sf, sk

def quick_load_xy(path, delimiter=",", skip_header=1):
    """
    Read a file using numpy.genfromtext to x, y
    """
    data = np.genfromtxt(path, delimiter=delimiter, skip_header=skip_header)
    return data[:, 0], data[:, 1]

def list_all_files(base_dir, extension='csv'):
    """
    Walk through a directory and list all files that end with the extension in the directory tree
    """
    files = []
    for root, dirs, files in os.walk(base_dir):
        for name in files:
            if name.lower().endswith(extension):
                fname = os.path.join(root, name)
                files.append(fname)
    return files


def write_json(data_file, arb_data):
    with open(data_file, 'w') as dump:
        dump.write(json.dumps(arb_data, indent=0))


def read_json(data_file):
    source = open(data_file, 'r').read()
    data = json.loads(source)
    return data