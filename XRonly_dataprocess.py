import sys, os
import pickle
import numpy as np

from helpers import findind, get_data, get_pads, get_calib

if __name__=="__main__":
    pad = 0
    r_start = 380
    r_end = 400
    imData_bkgd = get_data(pad, 401, 0)
    # read calibration values and subpad locations
    x0, y0, D, rr, tt, twotheta_deg = get_calib(pad)
    pads = get_pads(pad)
    Npads = pads.shape[0]

    # read XRonly run numbers and loop over
    runs = np.genfromtxt("XRonly_runnumber.txt").astype(int)
    data = {}
    for run in runs:
        data[run] = {}
        img = 0
        while True:
            try:
                imData = get_data(pad, run, img)-imData_bkgd
            except IOError:
                break
                
            ny, nx = imData.shape
            ft = np.zeros(len(tt))
            ft_count = np.zeros(len(tt))
            for y in xrange(ny):
                for x in xrange(nx):
                    i_pad = -1
                    for i in xrange(Npads):
                        if pads[i,0]<x and x<pads[i,2] and pads[i,1]<y and y<pads[i,3]:
                            i_pad = i
                    if i_pad >= 0:
                        r = np.sqrt((x-x0)**2 + (y-y0)**2)
                        if r_start <= r and r <= r_end: 
                            t = np.arctan((y-y0) / (x-x0))
                            ft[findind(t,tt)] += imData[y,x]
                            ft_count[findind(t,tt)] += 1
            ft = ft.astype(float) / ft_count
            data[run][img] = ft

            tt_deg = tt / np.pi * 180.

            print "\tCompleted img %d" % img
            img += 1
       
        ## END img loop

        print "Completed run %d" % run
    ## END run loop 

    data['tt'] = tt
    data['pad'] = pad
    pickle.dump(data, open("XRonly_p%d_data.pickle" % pad, 'wb'))
