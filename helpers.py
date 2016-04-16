import pickle
import os, sys
import numpy as np
import PIL.Image

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

def get_datadir():
    try:
        return pickle.load(open("parameters.pickle", 'r'))['data_dir']
    except KeyError:
        print "Error: Need to specify data_dir in \"parameters.pickle\""
        sys.exit(0)

def get_pads(pad):
    try:
        return pickle.load(open("parameters.pickle", 'r'))[pad]['pads']
    except KeyError:
        print "Error: Need to specify pads in \"parameters.pickle\""
        sys.exit(0)

def get_data(pad, run, img):
    data_dir = get_datadir()
    imFile = PIL.Image.open(
        data_dir+"/image_cspad%02d_r%04d_e%05d.tif" 
        % (pad, run, img))
    nx, ny = imFile.size
    imData = np.array(imFile.getdata()).reshape((ny, nx)).astype(float)

    return imData
    

## read calibration data
def get_calib(pad):
    params = pickle.load(open("calibration.pickle", 'rb'))[pad]
    x0 = params['x0']
    y0 = params['y0']
    D = params['D']
    rr = params['r_array']
    tt = params['t_array']
    twotheta_deg = params['twotheta_deg']

    return x0, y0, D, rr, tt, twotheta_deg
