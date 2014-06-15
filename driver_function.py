"""
Analyze spectral data using combination of numpy, scipy, peak-o-mat and some simple algorithms
@author: Rohan Isaac
"""

import matplotlib.pyplot as plt
#import guiqwt.pyplot as  # alternative plot
from spectra import Spectra

a =10
# Driver function    
#making it a script so it is easier to debug
plt.close('all')
print "Starting script ... "

#S = Spectra('samples/neon72.txt')
#S = Spectra('samples/80K-426c.txt')
S = Spectra('samples/PMDA-KBr-1-75x1.CSV')
   
S.find_background()  

# Plot data and bg    
plt.figure(1)
plt.plot(S.x,S.active,'-')
plt.plot(S.x,S.bg,'-r')

S.subtract_background()

plt.figure(4)
plt.plot(S.x,S.active,'b-')    
S.remove_spikes()
plt.plot(S.x,S.active,'r-')

S.guess_peak_width(max_width=50)
S.find_peaks(limit=30)
S.build_model()    

# plot spectra, model, found_peaks
plt.figure(2)
plt.plot(S.x,S.active,'-')
plt.plot(S.x,S.model_data,'r-')    
plt.plot(S.x[S.peak_pos],S.active[S.peak_pos],'oy')


S.fit_data()

# plot spectra,fit
plt.figure(3)
plt.plot(S.x,S.active,'g-',alpha=0.3)
plt.plot(S.x,S.model_data,'r-')

S.output_results()
