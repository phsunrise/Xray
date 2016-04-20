import sys, os
from getopt import getopt
import pickle
import numpy as np
import matplotlib.pyplot as plt
from helpers import findind, get_data, get_pads, get_calib
from adjust_subpads import adjust_subpads

if __name__=="__main__":
    # default values    
    pad = 0
    run = 336
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
    fig = plt.figure(figsize=(10,10))
    ax1 = fig.add_subplot(1,1,1)
    img = 0
    while True:
        try:
            imData = get_data(pad, run, img)
        except IOError:
            break

        figname = "r%04d_i%02d_rb%04d_p%d_offsets" % (run, img, run_bkgd, pad)
        (imData_polar, onPad, rr, tt, \
             twotheta_deg, fr, offset) = adjust_subpads(imData, pad, \
                                           fixpad=0, do_debug=False, \
                                           figname=figname)
        ax1.semilogy(twotheta_deg, fr, label="img%d" % img)

        img += 1

    ax1.legend()
    plt.show()
