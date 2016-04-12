import numpy as np
import matplotlib.pyplot as plt
import pickle

run = 360 # AB+EF
pad = 2
img = 1
run_bkgd = 373
data1 = pickle.load(open("r%3d_p%1d_i%1d_rb%3d.pickle" % (run, pad, img, run_bkgd), 'rb'))

run = 375 # AB
data2 = pickle.load(open("r%3d_p%1d_i%1d_rb%3d.pickle" % (run, pad, img, run_bkgd), 'rb'))

run = 377 # EF 
data3 = pickle.load(open("r%3d_p%1d_i%1d_rb%3d.pickle" % (run, pad, img, run_bkgd), 'rb'))

plt.plot(data1['twotheta']/np.pi*180, data1['fr'], 'b-',\
         label='AB+EF')
plt.plot(data2['twotheta']/np.pi*180, data2['fr'], 'r-',\
         label='AB')
plt.plot(data3['twotheta']/np.pi*180, data3['fr'], 'g-',\
         label='EF')

plt.xlabel(r"$2\theta$ (deg)")
plt.legend()
plt.savefig("compare_abef.pdf")
plt.show()
