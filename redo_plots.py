import os, sys
import numpy as np
import matplotlib.pyplot as plt
import pickle
from process_2 import process, findind, fitfunc

data = pickle.load(open("data_2.pickle", 'rb'))
if data['pad'] != 2:
    print "Wrong pad! Exiting..."
    sys.exit(1)

pad = 2
run_bkgd = 401
do_debug = False
bkgdSubtract = True

for run in data.keys():
    if run == 'pad':
        continue
    
    rundata = data[run]
    img = rundata['img'] 
    params = rundata['params']

    imData_polar, rr, tt, twotheta, fr = process(pad, run, img, run_bkgd, \
                                         do_debug, bkgdSubtract)
    twotheta_deg = twotheta/np.pi*180

    fig = plt.figure(figsize=(7.5, 15))
    ax1 = fig.add_subplot(2,1,1)
    ax1.imshow(imData_polar, \
            extent=(twotheta_deg[0], twotheta_deg[-1], tt[0]/np.pi*180, tt[-1]/np.pi*180),\
            aspect='auto', origin='lower')
    ax1.set_ylabel(r"$\phi$ (deg)")

    ax2 = fig.add_subplot(2,1,2, sharex=ax1)
    ax2.plot(twotheta_deg, fr, label='original')

    # plot fit curves
    fit_s = findind(340, rr)
    fit_e = findind(650, rr)
    ax2.plot(twotheta_deg[fit_s:fit_e], \
             fitfunc(rr[fit_s:fit_e], *params), \
             label='fit')    # fit curve
    ax2.plot(twotheta_deg[fit_s:fit_e], \
             params[0]*np.exp(-params[1]*rr[fit_s:fit_e])+params[2], \
             label='background')
            # background
    for i in range(3, len(params), 3):
        ax2.plot(twotheta_deg[fit_s:fit_e], \
                params[i]*np.exp(-0.5*((rr[fit_s:fit_e]-params[i+1])/params[i+2])**2),
                label='peak %d' % (i/3)) # peaks

    ax2.set_xlabel(r"$2\theta$ (deg)")
    ax2.legend()
    ax2.set_xlim(twotheta_deg[0], twotheta_deg[-1])

    ax1.set_title("Shot #%d" % run, fontsize=20)
    fig.savefig("r%03d_i%01d_rb%03d_redo.pdf" % (run, img, run_bkgd))
    fig.savefig("r%03d_i%01d_rb%03d_redo.png" % (run, img, run_bkgd))
