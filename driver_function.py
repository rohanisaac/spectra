"""
Analyze spectral data using combination of numpy, scipy, peak-o-mat and some simple algorithms
@author: Rohan Isaac
"""

import matplotlib.pyplot as plt
#import guiqwt.pyplot as  # alternative plot
from spectra import Spectra

# Driver function    
def main():
    plt.close('all')
    print "Starting script ... "
    
    print "Loading file ... "    
    #S = Spectra('samples/neon72.txt')
    #S = Spectra('samples/80K-426c.txt')
    S = Spectra('samples/PMDA-KBr-1-75x1.CSV')
    print S.data_max, S.data_max_pos, S.test_peak_width

    print "Finding background ... "    
    S.find_background()    
    # Plot data and bg    
    plt.figure(1)
    plt.plot(S.s.x,S.s.y,'-')
    plt.plot(S.bg.x,S.bg.y,'-r')
    plt.plot(S.s.x[S.slmin],S.s.y[S.slmin],'o')
    
    print "Subtracting background ... "
    S.subtract_background()
    
    y2 = S.remove_spikes()
    plt.figure(4)
    plt.plot(S.s.x,y2,'-')

    """print "Looking for peaks ... "  
    S.find_peaks()
    print "Found " + str(S.num_peaks) + " peaks."
    print S.peak_pos
    
    print "Building model ... "
    S.build_model()    
    """
    # plot spectra, model, found_peaks
    plt.figure(2)
    plt.plot(S.s.x,S.s.y,'-')
    """    
    plt.plot(S.x,S.y,'b-')    
    plt.plot(S.s.x[S.peak_pos],S.s.y[S.peak_pos],'or')

    print "Fitting Data"
    S.fit_data()
    print S.res[-1]
    
    # plot spectra,fit
    plt.figure(3)
    plt.plot(S.s.x,S.s.y,'ko',alpha=0.3)
    plt.plot(S.x,S.y,'r-')

    S.output_results()
    """
if __name__ == "__main__":
    main()
    