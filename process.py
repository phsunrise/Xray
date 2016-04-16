import sys, os
from getopt import getopt
import pickle
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib.widgets import SpanSelector, RectangleSelector 
from scipy.optimize import curve_fit
from helpers import findind, get_data, get_pads, get_calib
from adjust_subpads import adjust_subpads

def fitfunc(x, *params):
# parameter format:
# background (exponential): params[0]*np.exp(-params[1]*x)+params[2]
# starting from params[3], gaussian peaks: 
#     params[i]*np.exp(-0.5*((x-params[i+1])/params[i+2])**2)
    y = params[0]*np.exp(-params[1]*x)+params[2]
    for i in range(3, len(params), 3):
        y = y + params[i]*np.exp(-0.5*((x-params[i+1])/params[i+2])**2)
    return y

def process(pad, run, img, run_bkgd, do_debug, bkgdSubtract):
    return 0

if __name__=="__main__":
    # default values    
    pad = 2
    run = 360
    img = 1
    run_bkgd = 401
    img_bkgd = 0
    do_debug = False 
    bkgdSubtract = True

    # parse command line arguments 
    opts, args = getopt(sys.argv[1:], "p:r:i:b:d")
    for opt, arg in opts:
        if opt == '-p':
            pad = int(arg)
        elif opt == '-d':
            do_debug = True
        elif opt == '-r':
            run = int(arg)
        elif opt == '-i':
            img = int(arg)
        elif opt == '-b':
            run_bkgd = int(arg)
            bkgdSubtract = True
    
    # process (adjust pad offset)
    imData = get_data(pad, run, img) 
    imData_bkgd = get_data(pad, run_bkgd, img_bkgd) 
    ny, nx = imData.shape

    if bkgdSubtract:
        imData = imData - imData_bkgd

    fixpad = 0
    imData_polar, rr, tt, twotheta_deg, fr, offset = adjust_subpads(imData, pad, fixpad, do_debug)
    
    # plot data after process
    fig = plt.figure(figsize=(7.5, 15))
    ax1 = fig.add_subplot(2,1,1)
    ax1.imshow(imData_polar, \
            extent=(rr[0], rr[-1], tt[0]/np.pi*180, tt[-1]/np.pi*180),\
            aspect='auto', origin='lower')
    ax1.set_xlabel(r"$r$")
    ax1.set_ylabel(r"$\phi$ (deg)")

    ax2 = fig.add_subplot(2,1,2, sharex=ax1)
    ax2.plot(rr, fr, label='orig.')

    if do_debug:
        ax2.legend()
        plt.show()
        print "Debug mode: exiting before fit"
        sys.exit(0)

    ## let user select range for background fit
    ## first select fit range
    fig_fit = plt.figure(figsize=(20,10))
    ax1_fit = fig_fit.add_subplot(1,1,1)
    ax1_fit.plot(twotheta_deg, fr)
    ax1_fit.set_xlabel(r"$2\theta$ (deg)")
    ax1_fit.set_title("Please select fit range, then close window.")
    fit_range = [0,0]
    def onselect1(thmin, thmax):
        global fit_range
        indmin, indmax = np.searchsorted(twotheta_deg, (thmin, thmax))
        indmax = min(len(twotheta_deg)-1, indmax)
        fit_range = [indmin, indmax]

    span1 = SpanSelector(ax1_fit, onselect1, 'horizontal', useblit=True,\
                        rectprops=dict(alpha=0.5, facecolor='r'), \
                        span_stays=True)
    plt.show()

    ## next select peak range(s) to exclude
    while True:
        fig_fit = plt.figure(figsize=(20,10))
        ax1_fit = fig_fit.add_subplot(1,1,1)
        ax1_fit.plot(twotheta_deg[fit_range[0]:fit_range[1]], \
                     fr[fit_range[0]:fit_range[1]])
        ax1_fit.set_xlabel(r"$2\theta$ (deg)")
        ax1_fit.set_title("Please select peak range(s) to exclude, then close window.")
        peaks_range = []
        def onselect2(thmin, thmax):
            global fit_range
            indmin, indmax = np.searchsorted(twotheta_deg, (thmin, thmax))
            indmax = min(len(twotheta_deg)-1, indmax)
            peaks_range.append([indmin, indmax])
            ax1_fit.axvspan(twotheta_deg[peaks_range[-1][0]], twotheta_deg[peaks_range[-1][1]], \
                            facecolor='g', alpha=0.5)

        span2 = SpanSelector(ax1_fit, onselect2, 'horizontal', useblit=True,\
                            rectprops=dict(alpha=0.5, facecolor='r'), \
                            span_stays=False)
        plt.show()
        
        val = raw_input("Continue? (Enter n to redo this step: )")
        if val not in ['n', 'N']:
            break 

    # first fit background
    fit_mask = np.zeros(len(fr))
    fit_mask[fit_range[0]:fit_range[1]] = 1
    for i in xrange(len(peaks_range)):
        fit_mask[peaks_range[i][0]:peaks_range[i][1]] = 0
    fitmin, fitmax = np.nonzero(fit_mask)[0][0], np.nonzero(fit_mask)[0][-1]
    bkgd_fit = lambda x,a,b,c: a*np.exp(-b*x)+c
    b0 = -np.log(fr[fitmin]/fr[fitmax])/(rr[fitmin]-rr[fitmax])
    a0 = fr[fitmin]/np.exp(-b0*rr[fitmin])
    c0 = 0.
    fit_mask = fit_mask.astype(bool)
    popt, pcov = curve_fit(bkgd_fit, rr[fit_mask], fr[fit_mask], p0=[a0, b0, c0])
    a0, b0, c0 = popt
    print "fit:", popt, np.sqrt(np.diag(pcov))

    # then fit peaks + background
    # Let user choose approximate peak positions (using rectangular regions)
    while True:
        params = [a0, b0, c0]
        fig_fit = plt.figure(figsize=(20,10))
        ax_fit = fig_fit.add_subplot(1,1,1)
        ax_fit.plot(twotheta_deg[fitmin:fitmax], fr[fitmin:fitmax], 'b-', \
                 label='orig.')
        ax_fit.plot(twotheta_deg[fitmin:fitmax], a0*np.exp(-b0*rr[fitmin:fitmax])+c0, \
                 'r-', label='bkgd fit')
        ax_fit.plot(twotheta_deg[fitmin:fitmax], \
                 fr[fitmin:fitmax]-a0*np.exp(-b0*rr[fitmin:fitmax])+c0, \
                 'g-', label='minus bkgd')
        ax_fit.set_xlabel(r"$2\theta$ (deg)")
        ax_fit.legend()
        def onselect3(eclick, erelease):
            x1, y1 = eclick.xdata, eclick.ydata 
            x2, y2 = erelease.xdata, erelease.ydata 
            m = 0.5*(x1+x2)
            s = 0.5*0.5*abs(x2-x1) # this is approximate
            A = abs(y2-y1)
            print "mean =", m, "sigma =", s, "amp =", A
            params.extend([A,m,s])
            ax_fit.add_patch(Rectangle([min(x1,x2), min(y1,y2)], \
                                       abs(x2-x1), abs(y2-y1), alpha=0.5,\
                                       color='g'))
        rect = RectangleSelector(ax_fit, onselect3, drawtype='box')
        plt.show()

        try:
            params, params_cov = curve_fit(fitfunc, twotheta_deg[fitmin:fitmax], \
                                     fr[fitmin:fitmax], p0=params)
        except RuntimeError:
            print "Fit failed! Please choose the peaks again..."
            continue

        # show the fit functions
        plt.plot(twotheta_deg[fitmin:fitmax], fr[fitmin:fitmax], \
                 label='orig.')
        plt.plot(twotheta_deg[fitmin:fitmax], \
                 fitfunc(twotheta_deg[fitmin:fitmax], *params), \
                 label='fit')
        plt.plot(twotheta_deg[fitmin:fitmax], \
                 params[0]*np.exp(-params[1]*twotheta_deg[fitmin:fitmax])+params[2], \
                 label='bkgd')
        for i in range(3, len(params), 3):
            plt.plot(twotheta_deg[fitmin:fitmax], \
                    params[i]*np.exp(-0.5*((twotheta_deg[fitmin:fitmax]-params[i+1])/params[i+2])**2),
                    label='peak%d' % (i/3))

        plt.xlabel(r"$2\theta$ (deg)")
        plt.legend()
        plt.show()

        val = raw_input("Continue? Enter n to redo this step: ")
        if val not in ['n', 'N']:
            break

    # save data
    do_save = raw_input("Save data? ")
    if do_save in ['y', 'yes']:
        if not os.path.isfile("data_2.pickle"):
            data = {'pad': 2}
            pickle.dump(data, open("data_2.pickle", 'wb'))
        data = pickle.load(open("data_2.pickle", 'rb'))
        thisrun = {'img': img, 'run_bkgd': run_bkgd, 'twotheta': twotheta,\
                   'fr': fr, 'bkgdSubtract': bkgdSubtract}
        thisrun['Npeaks'] = len(params)/3 - 1
        thisrun['params'] = params
        thisrun['params_err'] = np.sqrt(np.diag(params_cov))
        data[run] = thisrun
        pickle.dump(data, open("data_2.pickle", 'wb'))
