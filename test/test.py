import sys
import os
import matplotlib.pyplot as plt
sp_path = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
sys.path.append(sp_path)
import spectra as sp

s = sp.Spectra(os.path.join(sp_path, 'data', 'PhillipsXRD.txt'))
# print s.x, s.y
# f,ax = plt.subplots()
s.find_peaks(limit=8)
s.build_model(bg_ord=2)
s.fit_data()
print s.output_results()

# f,ax = plt.subplots()
plt.plot(s.x, s.y, 'b-',  # input data
         s.x[s.peak_pos], s.y[s.peak_pos], 'ro',  # peaks
         s.x, s.out.init_fit, 'r--',  # model
         s.x, s.out.best_fit, 'g-',  # fit data)
         )
plt.legend(['Data', 'Peaks', 'Model', 'Fit'])
plt.show()
