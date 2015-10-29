import sys
sys.path.append('..')
import spectra as sp
import matplotlib.pyplot as plt


s = sp.Spectra('../data/PhillipsXRD.txt')
#print s.x, s.y
#f,ax = plt.subplots()
plt.plot(s.x,s.y)
plt.show()
s.guess_peak_width()
s.find_peaks()

#f,ax = plt.subplots()
plt.plot(s.x[s.peak_pos], s.y[s.peak_pos] ,'ro')
plt.show()