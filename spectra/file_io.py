import os
import re
import numpy as np

def assure_path_exists(path):
    """
    Create directory structure of file if it doens't already exist
    """
    dir = os.path.dirname(path)
    if not os.path.exists(dir):
        os.makedirs(dir)

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
            if re.search('[a-df-zA-DF-Z]', lin):
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
                if re.search('\t', lin):
                    sep_char = '\t'
                elif re.search(',', lin):
                    sep_char = ','
                else:
                    print "Unknown separator character"
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

def batch_load(folder, extension=None, delimiter=None):
    """
    Load a batch of files in a folder into a 2-d numpy array

    Parameters
    ----------
    folder (string)
        folder path containing the text based data files, one data set per file
        will not look in subdirectories
    extension (string, optional)
        a three letter extension of the text file. Defaults to all txt or csv
        files in the directory
    delimiter (character)
        Single character corresponding to the datatype separating the columns in
        the text file. Use delimiter='fwf' for fixed width text format

    Returns
    -------
    x (1-d array)
        x values (common for all the files, taken from first file read)
    data (2-d array)
        array of shape=(num_files, data_points) containing all the y data
    f_names (list)
        list containing the filenames of the files, same order as data
    """
    # generate a list of names
    f_names = []
    for f in os.listdir(folder):
        f_lower = f.to_lower()
        if f_lower.endswith('csv') or f_lower.endswith('txt'):
            f_names.apped(f)
    x, y = getxy(os.path.join(folder, f_name[0]))

    # allocate space for data
    data = np.zeros(len(f_name), len(x))
    for i, f in enumerate(f_name):
        _, y = getxy(os.path.join(folder, f))
        data[i, :] = y

    return x, data, f_names
