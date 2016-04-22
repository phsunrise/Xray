import sys, os
from getopt import getopt
import pickle
import numpy as np
import matplotlib
matplotlib.use('pdf')
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
 
    runs = np.genfromtxt("XRonly_runnumber.txt").astype(int)
    for run in runs[::-1]:
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

        ft_sum = []
        for img in xrange(Nimg):
            imData = get_data(pad, run, img)
            figname = "r%04d_i%02d_rb%04d_p%d_XRonly" % (run, img, run_bkgd, pad)
            (imData_polar, onPad, rr, tt, \
                 twotheta_deg, fr, offset) = adjust_subpads(imData, pad, \
                                               fixpad=0, do_debug=False, \
                                               figname=figname, offset=0.)
                 
            tt_deg = tt / np.pi * 180.
            rr_start = findind(380, rr)
            rr_end = findind(400, rr)
            if img == 0:
                ft0 = (np.sum(imData_polar[:, rr_start:rr_end], axis=1)
                     / np.sum(onPad[:, rr_start:rr_end], axis=1))
                continue
            #ax1.imshow(imData_polar, \
            #    extent=[rr[0], rr[-1], tt_deg[0], tt_deg[-1]], aspect='auto')
            #ax2.imshow(imData_polar[:, rr_start:rr_end])
            #plt.show()
            #sys.exit(0)

            ft = (np.sum(imData_polar[:, rr_start:rr_end], axis=1)
                / np.sum(onPad[:, rr_start:rr_end], axis=1))
            ft = ft-ft0

            ax1.plot(tt_deg, ft+200.*img, c=cm.jet(img*1./Nimg), label="img%d" % img)
            # do fft, but need to take care of NaN and inf
            mask = np.isfinite(ft)
            ft_filtered = ft[mask] 
            ft_sum.append(np.sum(ft_filtered))
            ft_fft = np.abs(np.fft.fft(ft_filtered))**2
            ft_fft = ft_fft[1:len(ft_fft)/2] # only keeping the first half
            
            ax2.semilogy(ft_fft, c=cm.jet(img*1./Nimg), label="img%d" % img)

        #ax1.legend()
        ax1.set_xlim(tt_deg[mask][0], tt_deg[mask][-1])
        ax1.set_xlabel(r"$\phi$ (deg)")
        ax1.set_title("Shot %d" % run)

        ax2.set_xlim(0, len(ft_fft)-1)
        ax2.set_xticklabels([])
        ax2.set_yticklabels([])

        ax3.plot(ft_sum, 'ro--')
        ax3.set_xticklabels([])
        ax3.set_title("Sum")

        plt.tight_layout()
        fig.savefig("r%04d_rb%04d_p%d_XRonly.pdf" % (run, run_bkgd, pad))
        print "Done run %d" % run
        plt.close()
