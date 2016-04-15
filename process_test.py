import sys, os
from getopt import getopt
import pickle
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib.widgets import SpanSelector 
from scipy.optimize import curve_fit
from helpers import findind, get_data, get_pads, get_calib

def fitfunc(x, *params):
# parameter format:
# background: params[0]*np.exp(-params[1]*x)+params[2]
# starting from params[3], gaussian peaks: 
#     params[i]*np.exp(-0.5*((x-params[i+1])/params[i+2])**2)
    y = params[0]*np.exp(-params[1]*x)+params[2]
    for i in range(3, len(params), 3):
        y = y + params[i]*np.exp(-0.5*((x-params[i+1])/params[i+2])**2)
    return y

def process(pad, run, img, run_bkgd, do_debug, bkgdSubtract):
    # Read data file into an array
    imData = get_data(pad, run, img) 
    imData_bkgd = get_data(pad, run_bkgd, 0) 
    ny, nx = imData.shape

    # get position of the four pads
    pads = get_pads(pad) 

    if bkgdSubtract:
        imData = imData - imData_bkgd
    
    ## read calibration data
    x0, y0, D, intercept = get_calib(pad)

    ## convert to polar coordinates
    # find the extent by looking at the corners
    y_corners = np.array([0, 0, ny-1, ny-1])
    x_corners = np.array([0, nx-1, 0, nx-1])
    rmin = min(np.sqrt((x_corners-x0)**2 + (y_corners-y0)**2))
    rmax = max(np.sqrt((x_corners-x0)**2 + (y_corners-y0)**2))
    nr = 500
    rr = np.linspace(rmin, rmax, nr, endpoint=False) # r coordinates
    twotheta = np.arctan((rr-intercept)/D)
    tmin = min(np.arctan((y_corners-y0) / (x_corners-x0)))
    tmax = max(np.arctan((y_corners-y0) / (x_corners-x0)))
    nt = 500
    tt = np.linspace(tmin, tmax, nt, endpoint=False) # t coordinates

    imData_polar = np.zeros((4, nt, nr))
    onPad = np.zeros((4, nt, nr))
    for y in range(ny):
        for x in range(nx):
            pad = -1
            for i in xrange(4):
                if pads[i,0]<x and x<pads[i,2] and pads[i,1]<y and y<pads[i,3]:
                    pad = i
            if pad >= 0:
                r = np.sqrt((x-x0)**2 + (y-y0)**2)
                t = np.arctan((y-y0) / (x-x0))
                imData_polar[pad,
                             findind(t,tt), \
                             findind(r,rr)] += imData[y,x]
                onPad[pad,
                      findind(t,tt), \
                      findind(r,rr)] += 1 

    if do_debug:
        # Create a figure
        fig = plt.figure()

        ## Display the image
        ax1 = fig.add_subplot(2,2,1)
        ax1.imshow(imData, origin='lower')
        for i in xrange(4):
            ax1.add_patch(Rectangle((pads[i,0], pads[i,1]), \
                        pads[i,2]-pads[i,0], pads[i,3]-pads[i,1],\
                        fill=False, ls='--', color='r'))

        ## display image in polar coordiantes before adjusting for offset
        ax2 = fig.add_subplot(2,2,2)
        ax2.imshow(np.sum(imData_polar, axis=0), \
                extent=(rmin, rmax, tmin/np.pi*180, tmax/np.pi*180),\
                aspect='auto', origin='lower')

        ### print out pad ranges as index in rr
        #print "pad0:"
        #print findind(np.sqrt((pads_x[0]-x0)**2+(pads_y[1]-y0)**2), rr)
        #print findind(np.sqrt((pads_x[1]-x0)**2+(pads_y[0]-y0)**2), rr)
        #print "pad1:"
        #print findind(np.sqrt((pads_x[2]-x0)**2+(pads_y[1]-y0)**2), rr)
        #print findind(np.sqrt((pads_x[3]-x0)**2+(pads_y[0]-y0)**2), rr)
        #print "pad2:"
        #print findind(np.sqrt((pads_x[0]-x0)**2+(pads_y[3]-y0)**2), rr)
        #print findind(np.sqrt((pads_x[1]-x0)**2+(pads_y[2]-y0)**2), rr)
        #print "pad3:"
        #print findind(np.sqrt((pads_x[2]-x0)**2+(pads_y[3]-y0)**2), rr)
        #print findind(np.sqrt((pads_x[3]-x0)**2+(pads_y[2]-y0)**2), rr)


    twotheta_deg = twotheta/np.pi*180

    fr = np.zeros((4, nr))
    for i_pad in xrange(4):
        fr[i_pad] = np.sum(imData_polar[i_pad], axis=0)/np.sum(onPad[i_pad], axis=0)

    if do_debug:
        # display curves before adjusting for offset
        ax3 = fig.add_subplot(2,2,3)
        
        for i_pad in xrange(4):
            ax3.plot(twotheta_deg, fr[i_pad], label='%d'%i_pad)
        ax3.set_xlabel(r"$2\theta$")
        ax3.legend()

    ## calculate offsets 
    # offset1 = pad2 - pad3
    offset1 = np.mean((fr[2][140:220]-fr[3][140:220]))
    print "offset1:", offset1
    # offset2 = pad3 - pad0
    offset2 = np.mean((fr[3][200:350]-fr[0][200:350]))
    print "offset2:", offset2
    # offset3 = pad0 - pad1
    offset3 = np.mean((fr[0][265:310]-fr[1][265:310]))
    print "offset3:", offset3
    # offset4 = pad1[-10:] - 0.
    offset4 = np.mean(fr[1][487:496])
    print "offset4:", offset4
    # plot fr after adjusting for offset
    fr[1] -= offset4
    fr[0] -= offset3 + offset4
    fr[3] -= offset2 + offset3 + offset4
    fr[2] -= offset1 + offset2 + offset3 + offset4

    if do_debug:
        ax4 = fig.add_subplot(2,2,4)
        for i_pad in xrange(4):
            ax4.plot(twotheta_deg, fr[i_pad], label='%d'%i_pad)
        ax4.set_xlabel(r"$2\theta$")
        ax4.legend()
        fig.savefig("r%03d_i%01d_rb%03d_debug.pdf" % (run, img, run_bkgd))

    # now redo coordinate change with offsets
    imData_polar = np.zeros((nt, nr))
    onPad = np.zeros((nt, nr))
    for y in range(ny):
        for x in range(nx):
            pad = -1
            for i in xrange(4):
                if pads[i,0]<x and x<pads[i,2] and pads[i,1]<y and y<pads[i,3]:
                    pad = i

            if pad >= 0:
                r = np.sqrt((x-x0)**2 + (y-y0)**2)
                t = np.arctan((y-y0) / (x-x0))
                imData_polar[findind(t,tt), \
                             findind(r,rr)] += imData[y,x]
                onPad[findind(t,tt), \
                      findind(r,rr)] += 1 

    fr = np.sum(imData_polar, axis=0)/np.sum(onPad, axis=0)

    # return values: data in polar coordinates; 2*theta array; r array;
    #                fr array (average for each r)
    return imData_polar, rr, tt, twotheta, fr



if __name__=="__main__":
    # default values    
    pad = 2
    run = 360
    img = 1
    if img == 0:
        run_bkgd = 334
    else:
        run_bkgd = 373
        #run_bkgd = 390
    #run_bkgd = 334
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
    imData_polar, rr, tt, twotheta, fr = process(pad, run, img, run_bkgd, \
                                         do_debug, bkgdSubtract)
    twotheta_deg = twotheta/np.pi*180
    
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
        
        val = raw_input("Was the range okay? ")
        if val in ['y', 'Y', 'yes']:
            break 

    # first fit background
    fit_mask = np.zeros(len(fr))
    fit_mask[fit_range[0]:fit_range[1]] = 1
    for i in xrange(len(peaks_range)):
        fit_mask[peaks_range[i][0]:peaks_range[i][1]] = 0
    maskmin, maskmax = np.nonzero(fit_mask)[0][0], np.nonzero(fit_mask)[0][-1]
    bkgd_fit = lambda x,a,b,c: a*np.exp(-b*x)+c
    b0 = -np.log(fr[maskmin]/fr[maskmax])/(rr[maskmin]-rr[maskmax])
    a0 = fr[maskmin]/np.exp(-b0*rr[maskmin])
    c0 = 0.
    fit_mask = fit_mask.astype(bool)
    popt, pcov = curve_fit(bkgd_fit, rr[fit_mask], fr[fit_mask], p0=[a0, b0, c0])
    a0, b0, c0 = popt
    print "fit:", popt, np.sqrt(np.diag(pcov))

    plt.plot(twotheta_deg[maskmin:maskmax], fr[maskmin:maskmax], 'b-', \
             label='orig.')
    plt.plot(twotheta_deg[maskmin:maskmax], a0*np.exp(-b0*rr[maskmin:maskmax])+c0, \
             'r-', label='bkgd fit')
    plt.xlabel(r"$2\theta$ (deg)")
    plt.show()
    sys.exit()

    # now fit peaks + background
    # Let user choose approximate peak positions (using rectangular regions)
    peaks = []
    def onselect(eclick, erelease):
        x1, x2 = 
        m = (ereleast.xdata-eclick.xdata

    # setting initial values
    m1 = 390 # mean
    m1_ind = findind(m1, rr)
    s1 = 10 # standard deviation
    A1 = .5*fr[m1_ind]-a0*np.exp(-b0*rr[m1_ind])-c0 # amplitude
    m2 = 420 # mean
    m2_ind = findind(m2, rr)
    s2 = 10 # standard deviation
    A2 = .5*fr[m1_ind]-a0*np.exp(-b0*rr[m2_ind])-c0 # amplitude
    m3 = 520 # mean
    m3_ind = findind(m3, rr)
    s3 = 10 # standard deviation
    A3 = fr[m3_ind]-a0*np.exp(-b0*rr[m3_ind])-c0 # amplitude
    # see fitfunc for definition of parameter list
    if Npeaks == 3:
        params0 = [a0, b0, c0, A1, m1, s1, A2, m2, s2, A3, m3, s3]
    elif Npeaks == 2: 
        #params0 = [a0, b0, c0, A1, m1, s1, A2, m2, s2]
        params0 = [a0, b0, c0, A1, m1, s1, A3, m3, s3]
    elif Npeaks == 1:
        params0 = [a0, b0, c0, A1, m1, s1]
    fit_s = fit_s1 
    fit_e = fit_e2
    params, pcov = curve_fit(fitfunc, rr[fit_s:fit_e], fr[fit_s:fit_e], p0=params0)
    print params
    ax2.plot(rr[fit_s:fit_e], fitfunc(rr[fit_s:fit_e], *params), label='fit')

    ## plot each curve
    ax2.plot(rr[fit_s:fit_e], \
             params[0]*np.exp(-params[1]*rr[fit_s:fit_e])+params[2], \
             label='background')
            # background
    for i in range(3, len(params), 3):
        ax2.plot(rr[fit_s:fit_e], \
                params[i]*np.exp(-0.5*((rr[fit_s:fit_e]-params[i+1])/params[i+2])**2),
                label='peak%d' % (i/3)) # peaks

    ax2.set_xlabel(r"$r$")
    ax2.legend()
    fig.savefig("r%03d_i%01d_rb%03d.pdf" % (run, img, run_bkgd))
    plt.show()

    # save data
    do_save = raw_input("Save data? ")
    if do_save in ['y', 'yes']:
        if not os.path.isfile("data_2.pickle"):
            data = {'pad': 2}
            pickle.dump(data, open("data_2.pickle", 'wb'))
        data = pickle.load(open("data_2.pickle", 'rb'))
        thisrun = {'img': img, 'run_bkgd': run_bkgd, 'twotheta': twotheta,\
                   'fr': fr, 'bkgdSubtract': bkgdSubtract}
        thisrun['Npeaks'] = Npeaks
        thisrun['params'] = params
        data[run] = thisrun
        pickle.dump(data, open("data_2.pickle", 'wb'))
