import os, sys
import pickle
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import linregress

test = True 
fit = False 
compare = False 

test = int(sys.argv[1])

run = 375
pad = 0
img = 0
run_bkgd = 334
data1 = pickle.load(open("r%3d_p%1d_i%1d_rb%3d.pickle" % (run, pad, img, run_bkgd), 'rb'))
sin1 = np.sin(0.5*data1['twotheta'])
fr1 = data1['fr']

pad = 2
img = 1
run_bkgd = 373
data2 = pickle.load(open("r%3d_p%1d_i%1d_rb%3d.pickle" % (run, pad, img, run_bkgd), 'rb'))
sin2 = np.sin(0.5*data2['twotheta'])
fr2 = data2['fr']
if compare: # compare with pad 2 before shot
    data3 = pickle.load(open("r%3d_p%1d_i%1d_rb%3d.pickle" % (run, pad, 0, 334), 'rb'))
    sin3 = np.sin(0.5*data3['twotheta'])
    fr3 = data3['fr']

if fit:
    ## fit background (assume exponential)
    start_ind = 250
    end_ind = 400
    slope, intercept, r_val, p_val, stderr = linregress(
             np.tan(data1['twotheta'][start_ind:end_ind]), np.log(fr2[start_ind:end_ind]))
    fr2_bkgd = np.exp(slope * np.tan(data1['twotheta']) + intercept)

fig = plt.figure()

if test:
    ax = fig.add_subplot(1,2,1)
    ax.plot(sin1, fr1, 'b-')

    ax = fig.add_subplot(1,2,2)
    if fit:
        ax.plot(sin2, fr2, 'b-')
        ax.plot(sin2, fr2_bkgd, 'r--')
        ax.plot(sin2, fr2 - fr2_bkgd, 'r-')
    else:
        ax.plot(fr2, 'b-')


else:
    # first find the peaks
    p11 = sin1[160 + np.argmax(fr1[160:220])]
    p12 = sin1[250 + np.argmax(fr1[250:300])]
    p21_range = (40, 80)
    p21 = sin2[p21_range[0] + np.argmax(fr2[p21_range[0]:p21_range[1]])]
    p22_range = (100, 160)
    p22 = sin2[p22_range[0] + np.argmax(fr2[p22_range[0]:p22_range[1]])]

    ratio = 0.5*((p21/p11) + (p22/p12))
    #ratio = p22/p12
    ax = fig.add_subplot(1,1,1)
    ax.plot(sin1*ratio, fr1, 'b-', label='before')
    if fit:
        ax.plot(sin2, fr2-fr2_bkgd, 'r-', label='after')
    else:
        ax.plot(sin2, fr2, 'r-', label='after')

    ymin, ymax = ax.get_ylim()
    ax.vlines([p11*ratio, p12*ratio], ymin, ymax, \
              colors='k', linestyles='--')

    if compare: # compare with pad2 before shot
        ax.plot(sin3, fr3/np.nanmax(fr3)*np.nanmax(fr1), 'g-', label='pad 2 before')

    ax.set_xlim(0.3, 0.5)
    ax.set_ylim(ymin, ymax)
    ax.set_xlabel(r"$\sin\theta$")
    ax.set_title(r"$\sin\theta/\sin\theta_0$=%.2f" % ratio)
    ax.legend()

plt.show()
if not compare:
    fig.savefig("r%3d_compression.pdf" % run)
