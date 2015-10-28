import sys
sys.path.append('..')
import spectra as sp
import matplotlib.pyplot as plt

s = sp.Spectra('../data/PhillipsXRD.txt')

print s.x, s.y
plt.figure()
plt.plot(s.x,s.y)
plt.show()