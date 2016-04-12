import numpy as np
import matplotlib.pyplot as plt
import pickle

run = 360
pad = 2
img = 0
if img == 0:
    run_bkgd = 334
else:
    run_bkgd = 373
data1 = pickle.load(open("r%3d_p%1d_i%1d_noBkgdSubtract.pickle" % (run, pad, img), 'rb'))
data2 = pickle.load(open("r%3d_p%1d_i%1d_noBkgdSubtract.pickle" % (run_bkgd, pad, 0), 'rb'))
data3 = pickle.load(open("r%3d_p%1d_i%1d_rb%3d.pickle" % (run, pad, img, run_bkgd), 'rb'))

plt.plot(data1['twotheta']/np.pi*180, data1['fr'], 'b-',\
         label='original')
plt.plot(data2['twotheta']/np.pi*180, data2['fr'], 'r-',\
         label='background')
plt.plot(data3['twotheta']/np.pi*180, data3['fr'], 'g-',\
         label='orig. - backg.')

plt.xlabel(r"$2\theta$ (deg)")
plt.legend()
plt.show()
plt.savefig("compare_abef.pdf")
