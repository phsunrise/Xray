import sys, os
from getopt import getopt
import pickle
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit, leastsq
import PIL.Image
from fitcircle import fitcircle

def findind(x, x_array): # assume x_array is evenly spaced
    if len(x_array) <= 1:
        return 0
    dx = x_array[1]-x_array[0] 
    xx = int((x-x_array[0]) / dx)
    if xx < 0:
        return 0
    elif xx >= len(x_array):
        return len(x_array)-1
    else:
        return xx 

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
    imFile = PIL.Image.open(
        "LA61_MgO_analysis/MgO_CSPAD/image_cspad%02d_r%04d_e%05d.tif" 
        % (pad, run, img))
    imFile_bkgd = PIL.Image.open(
        "LA61_MgO_analysis/MgO_CSPAD/image_cspad%02d_r%04d_e%05d.tif" 
        % (pad, run_bkgd, 0)) # assume this to be the background
    nx, ny = imFile.size
    imData = np.array(imFile.getdata()).reshape((ny, nx)).astype(float)
    imData_bkgd = np.array(imFile_bkgd.getdata()).reshape((ny, nx)).astype(float)

    # these two arrays define the position of the four pads
    pads_x = [12, 190, 197, 385] 
    pads_y = [0, 184, 212, 395]
    #fig = plt.figure()
    ## pad 1
    #ax = fig.add_subplot(2,2,1)
    #pad = imData[pads_y[0]:pads_y[1], pads_x[0]:pads_x[1]]
    #pad_bkgd = imData_bkgd[pads_y[0]:pads_y[1], pads_x[0]:pads_x[1]]
    #ax.plot(np.sum(pad, axis=0)+np.sum(pad_bkgd[:, 10])-np.sum(pad[:, 10]), 'r-')
    #ax.plot(np.sum(pad_bkgd, axis=0), 'b-')
    ## pad 2 
    #ax = fig.add_subplot(2,2,2)
    #pad = imData[pads_y[0]:pads_y[1], pads_x[2]:pads_x[3]]
    #pad_bkgd = imData_bkgd[pads_y[0]:pads_y[1], pads_x[2]:pads_x[3]]
    #ax.plot(np.sum(pad, axis=0)+np.sum(pad_bkgd[:, 10])-np.sum(pad[:, 10]), 'r-')
    #ax.plot(np.sum(pad_bkgd, axis=0), 'b-')
    ## pad 3 
    #ax = fig.add_subplot(2,2,3)
    #pad = imData[pads_y[2]:pads_y[3], pads_x[0]:pads_x[1]]
    #pad_bkgd = imData_bkgd[pads_y[2]:pads_y[3], pads_x[0]:pads_x[1]]
    #ax.plot(np.sum(pad, axis=0), 'r-')
    #ax.plot(np.sum(pad_bkgd, axis=0), 'b-')
    ## pad 4 
    #ax = fig.add_subplot(2,2,4)
    #pad = imData[pads_y[2]:pads_y[3], pads_x[2]:pads_x[3]]
    #pad_bkgd = imData_bkgd[pads_y[2]:pads_y[3], pads_x[2]:pads_x[3]]
    #ax.plot(np.sum(pad, axis=0), 'r-')
    #ax.plot(np.sum(pad_bkgd, axis=0), 'b-')
    #
    #plt.show()

    if bkgdSubtract:
        imData = imData - imData_bkgd

    ## subtract average background
    #imData = imData - np.mean(imData[0:150, 0:150])
    #imData[imData < 0.] = 0.
    
    ## read calibration data
    params = pickle.load(open("calibration_2.pickle", 'rb'))
    x0 = params['x0']
    y0 = params['y0']
    D = params['D']
    intercept = params['intercept']

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
            if pads_x[0]<x and x<pads_x[1]:
                if pads_y[0]<y and y<pads_y[1]:
                    pad = 0
                elif pads_y[2]<y and y<pads_y[3]:
                    pad = 2
            elif pads_x[2]<x and x<pads_x[3]:
                if pads_y[0]<y and y<pads_y[1]:
                    pad = 1 
                elif pads_y[2]<y and y<pads_y[3]:
                    pad = 3
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
        xmin, xmax = ax1.get_xlim()
        ymin, ymax = ax1.get_ylim()
        ax1.vlines(pads_x, ymin, ymax, colors='r', linestyles='--')
        ax1.hlines(pads_y, xmin, xmax, colors='r', linestyles='--')
        ax1.set_xlim(xmin, xmax)
        ax1.set_ylim(ymin, ymax)

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
            if pads_x[0]<x and x<pads_x[1]:
                if pads_y[0]<y and y<pads_y[1]:
                    pad = 0
                    imData[y,x] -= offset3 + offset4
                elif pads_y[2]<y and y<pads_y[3]:
                    pad = 2
                    imData[y,x] -= offset1 + offset2 + offset3 + offset4
            elif pads_x[2]<x and x<pads_x[3]:
                if pads_y[0]<y and y<pads_y[1]:
                    pad = 1 
                    imData[y,x] -= offset4
                elif pads_y[2]<y and y<pads_y[3]:
                    pad = 3
                    imData[y,x] -= offset2 + offset3 + offset4
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
    Npeaks = 2
    do_debug = False 
    bkgdSubtract = True

    # parse command line arguments 
    opts, args = getopt(sys.argv[1:], "r:i:b:p:d")
    for opt, arg in opts:
        if opt == '-d':
            do_debug = True
        elif opt == '-r':
            run = int(arg)
        elif opt == '-i':
            img = int(arg)
        elif opt == '-b':
            run_bkgd = int(arg)
            bkgdSubtract = True
        elif opt == '-p':
            Npeaks = int(arg)
    
    # process (adjust pad offset)
    imData_polar, rr, tt, twotheta, fr = process(pad, run, img, run_bkgd, \
                                         do_debug, bkgdSubtract)
    
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
        plt.show()
        print "Debug mode: exiting before fit"
        sys.exit(0)

    ## fit data
    if Npeaks == 3:
        fit_s1 = findind(340, rr)
        fit_e1 = findind(380, rr)
        fit_s2 = findind(570, rr)
        fit_e2 = findind(650, rr)
    elif Npeaks in [1, 2]:
        fit_s1 = findind(340, rr)
        fit_e1 = findind(380, rr)
        fit_s2 = findind(450, rr)
        fit_e2 = findind(600, rr)
    # first fit background
    fr_fit = lambda x,a,b,c: a*np.exp(-b*x)+c
    b0 = -np.log(fr[fit_e2]/fr[fit_s1])/(rr[fit_e2]-rr[fit_s1])
    a0 = fr[fit_s1]/np.exp(-b0*rr[fit_s1])
    c0 = 0.
    popt, pcov = curve_fit(fr_fit, np.concatenate((rr[fit_s1:fit_e1], rr[fit_s2:fit_e2])),\
                           np.concatenate((fr[fit_s1:fit_e1], fr[fit_s2:fit_e2])), \
                           p0=[a0, b0, c0])
    a0, b0, c0 = popt
    print "fit:", popt, np.sqrt(np.diag(pcov))
    bkgd_start = 6
    bkgd_end = 500
    #ax6.plot(rr[bkgd_start:bkgd_end], a0*np.exp(-b0*rr[bkgd_start:bkgd_end])+c0, \
    #         label='bkgd fit')
    #ax6.plot(rr[bkgd_start:bkgd_end], (fr-a0*np.exp(-b0*rr)-c0)[bkgd_start:bkgd_end], \
    #         'r-', label='bkgd fit')

    # now fit peaks + background
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
