import matplotlib.pyplot as plt
from .file_io import getxy

def plotxy(filename):
    x, y = getxy(filename)
    plt.plot(x, y)
    plt.title(filename)

def fit_plot_single(range):
    pass
