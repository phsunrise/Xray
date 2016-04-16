import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from helpers import findind, get_pads, get_calib

'''
This module adjusts for offsets on the subpads so that 
they agree on overlapping radial ranges 
'''

def adjust_subpads(imData, pad, fixpad=0, do_debug=False, figname='debug'):
# fixpad is the subpad to be fixed; other pads will be adjusted to this pad 
    ny, nx = imData.shape
    pads = get_pads(pad)
    Npads = pads.shape[0]
    x0, y0, D, rr, tt, twotheta_deg = get_calib(pad)
    nr = len(rr)
    nt = len(tt)

    ## convert data to polar coordinates
    imData_polar = np.zeros((Npads, nt, nr))
    onPad = np.zeros((Npads, nt, nr))
    for y in range(ny):
        for x in range(nx):
            i_pad = -1
            for i in xrange(Npads):
                if pads[i,0]<x and x<pads[i,2] and pads[i,1]<y and y<pads[i,3]:
                    i_pad = i
            if i_pad >= 0:
                r = np.sqrt((x-x0)**2 + (y-y0)**2)
                t = np.arctan((y-y0) / (x-x0))
                imData_polar[i_pad,
                             findind(t,tt), \
                             findind(r,rr)] += imData[y,x]
                onPad[i_pad,
                      findind(t,tt), \
                      findind(r,rr)] += 1 

    # Create a figure
    fig = plt.figure()

    ## Display the image
    ax1 = fig.add_subplot(2,2,1)
    ax1.imshow(imData, origin='lower')
    for i in xrange(Npads):
        ax1.add_patch(Rectangle((pads[i,0], pads[i,1]), \
                    pads[i,2]-pads[i,0], pads[i,3]-pads[i,1],\
                    fill=False, ls='--', color='r'))

    ## display image in polar coordiantes before adjusting for offset
    ax2 = fig.add_subplot(2,2,2)
    ax2.imshow(np.sum(imData_polar, axis=0), \
            extent=(rr[0], rr[-1], tt[0]/np.pi*180, tt[-1]/np.pi*180),\
            aspect='auto', origin='lower')

    fr = np.zeros((Npads, nr))
    for i_pad in xrange(Npads):
        fr[i_pad] = np.sum(imData_polar[i_pad], axis=0)/np.sum(onPad[i_pad], axis=0)

    # display curves before adjusting for offset
    ax3 = fig.add_subplot(2,2,3)
    
    for i_pad in xrange(Npads):
        ax3.plot(twotheta_deg, fr[i_pad], label='pad%d'%(i_pad+1))
    ax3.set_xlabel(r"$2\theta$ (deg)")
    ax3.legend()

    ## adjusting for offsets
    ## Method: starting from fixpad, find the subpad that has the most
    ##      overlap with it, and adjust that pad. In the adjusting
    ##      process we find again the pad with largest overlap
    onPad_r = np.any(onPad, axis=1)
        # this 2D array now records if a given r is covered by each pad
    offset = np.zeros(Npads)
        # record the offset value for each pad
    adjusted = np.zeros(Npads).astype(bool)
        # record if each pad has been adjusted
    adjusted[fixpad] = True
    current_pad = fixpad
    while not np.all(adjusted):
        # first find the unadjusted subpad that has the largest overlap
        # with current_pad
        overlap = 0
        adj_pad = fixpad 
        for i_pad in xrange(Npads):
            if adjusted[i_pad]:
                continue
            if np.sum(np.logical_and(onPad_r[current_pad], onPad_r[i_pad])) > overlap:
                adj_pad = i_pad
                overlap = np.sum(np.logical_and(onPad_r[current_pad], onPad_r[i_pad]))
        # next set current_pad to this subpad and find the subpad 
        # that has the largest overlap
        # with adj_pad
        overlap = 0
        current_pad = adj_pad
        adj_pad = fixpad 
        for i_pad in xrange(Npads):
            if not adjusted[i_pad]:
                continue
            if np.sum(np.logical_and(onPad_r[current_pad], onPad_r[i_pad])) > overlap:
                adj_pad = i_pad
                overlap = np.sum(np.logical_and(onPad_r[current_pad], onPad_r[i_pad]))
        # finally, adjust current_pad to match adj_pad
        adj_min = np.nonzero(np.logical_and(onPad_r[current_pad], onPad_r[adj_pad]))[0][0] + 5
        adj_max = np.nonzero(np.logical_and(onPad_r[current_pad], onPad_r[adj_pad]))[0][-1] - 5
        if adj_max <= adj_min:
            raise RuntimeError("Not enough overlap!")
        offset[current_pad] = np.mean(fr[adj_pad][adj_min:adj_max]-fr[current_pad][adj_min:adj_max])
        fr[current_pad] += offset[current_pad]
        adjusted[current_pad] = True
        print "Adjusted pad %d against pad %d" % (current_pad+1, adj_pad+1)
        
    # plot fr after adjusting for offset
    ax4 = fig.add_subplot(2,2,4)
    for i_pad in xrange(4):
        ax4.plot(twotheta_deg, fr[i_pad], label='pad%d'%(i_pad+1))
    ax4.set_xlabel(r"$2\theta$")
    ax4.legend()

    fig.savefig(figname+".pdf")
    if do_debug:
        plt.show()

    # now redo coordinate change with offsets
    imData_polar = np.zeros((nt, nr))
    onPad = np.zeros((nt, nr))
    for y in range(ny):
        for x in range(nx):
            i_pad = -1
            for i in xrange(4):
                if pads[i,0]<x and x<pads[i,2] and pads[i,1]<y and y<pads[i,3]:
                    i_pad = i

            if i_pad >= 0:
                r = np.sqrt((x-x0)**2 + (y-y0)**2)
                t = np.arctan((y-y0) / (x-x0))
                imData_polar[findind(t,tt), \
                             findind(r,rr)] += imData[y,x]+offset[i_pad]
                onPad[findind(t,tt), \
                      findind(r,rr)] += 1 

    fr = np.sum(imData_polar, axis=0)/np.sum(onPad, axis=0)

    # return values: data in polar coordinates; 2*theta array; r array;
    #                fr array (average for each r)
    return imData_polar, rr, tt, twotheta_deg, fr, offset
