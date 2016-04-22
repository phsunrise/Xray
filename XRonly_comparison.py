import sys, os
from getopt import getopt
import pickle
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm

from helpers import findind, get_data
from adjust_subpads import adjust_subpads

if __name__=="__main__":
    # default values    
    pad = 0
    run = 359 
    run_bkgd = 401
    img_bkgd = 0

    # parse command line arguments 
    opts, args = getopt(sys.argv[1:], "p:r:b:")
    for opt, arg in opts:
        if opt == '-p':
            pad = int(arg)
        elif opt == '-r':
            run = int(arg)
        elif opt == '-b':
            run_bkgd = int(arg)
 
    # plot all images
    fig = plt.figure(figsize=(10,20))
    ax1 = fig.add_subplot(3,1,1)
    ax2 = fig.add_subplot(3,1,2)
    ax3 = fig.add_subplot(3,1,3)
    
    # first get the total number of images for colormap
    Nimg = 0
    while True:
        try:
            get_data(pad, run, Nimg)
        except IOError:
            break
        Nimg += 1

    stddev = []
    for img in xrange(Nimg):
        imData = get_data(pad, run, img)
        figname = "r%04d_i%02d_rb%04d_p%d_XRonly" % (run, img, run_bkgd, pad)
        (imData_polar, onPad, rr, tt, \
             twotheta_deg, fr, offset) = adjust_subpads(imData, pad, \
                                           fixpad=0, do_debug=False, \
                                           figname=figname, offset=0.)
        tt_deg = tt / np.pi * 180.
        twotheta_start = findind(40.0, twotheta_deg)
        twotheta_end = findind(48.0, twotheta_deg)
        if img == 0:
            ft0 = (np.sum(imData_polar[:, twotheta_start:twotheta_end], axis=1)
                 / np.sum(onPad[:, twotheta_start:twotheta_end], axis=1))
            continue

        ft = (np.sum(imData_polar[:, twotheta_start:twotheta_end], axis=1)
            / np.sum(onPad[:, twotheta_start:twotheta_end], axis=1))
        ft = ft-ft0
        ax1.plot(tt_deg, ft+200.*img, c=cm.jet(img*1./Nimg), label="img%d" % img)
        # do fft, but need to take care of NaN and inf
        mask = np.isfinite(ft)
        ft_filtered = ft[mask] 
        stddev.append(np.std(ft_filtered))
        ft_fft = np.abs(np.fft.fft(ft_filtered))**2
        ft_fft = ft_fft[1:len(ft_fft)/2] # only keeping the first half
        
        ax2.semilogy(ft_fft, c=cm.jet(img*1./Nimg), label="img%d" % img)

    #ax1.legend()
    ax1.set_xlim(tt_deg[mask][0], tt_deg[mask][-1])
    ax1.set_xlabel(r"$\phi$ (deg)")
    ax2.set_xlim(0, len(ft_fft)-1)
    ax2.set_xticklabels([])
    ax2.set_yticklabels([])
    ax3.plot(stddev, 'ro--')
    ax3.set_xticklabels([])
    fig.savefig("r%04d_rb%04d_p%d_XRonly.pdf" % (run, run_bkgd, pad))
    plt.show()
