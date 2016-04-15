import sys, os
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import PIL.Image
import pickle

from fitcircle import fitconccircle
from roipoly import roipoly
from helpers import findind, get_data

# Read data file into an array
pad = 0
run = 389
img = 1 
run_bkgd = 401
if pad == 2:
    twotheta_deg = np.array([47.9, 57.4, 63.1]) # array for peaks
elif pad == 0:
    twotheta_deg = np.array([32.2, 34.4, 36.7, 47.9])
tantwotheta = np.tan(twotheta_deg/180.*np.pi)
imData = get_data(pad, run, img)
ny, nx = imData.shape
imData_bkgd = get_data(pad, run_bkgd, 0) 
imData = imData - imData_bkgd

# use ROI to select region
N_circ = len(twotheta_deg) 
ROI = []
mask = []
for i_circ in xrange(N_circ):
    fig = plt.figure(figsize=(10,10))
    ax = fig.add_subplot(1,1,1)
    ax.imshow(imData, origin='lower')
    for i in xrange(i_circ):
        ROI[i].displayROI()
    ax.set_title("Select region %d: right or double click to finish" % (i_circ+1))

    ROI.append(roipoly(roicolor='r'))
    mask.append(np.logical_and(ROI[i_circ].getMask(imData), \
                imData>0.5*np.max(imData[ROI[i_circ].getMask(imData)]))
               )

# fit to circle, weight by intensity
xx, yy = np.meshgrid(np.arange(nx), np.arange(ny))
X = np.array([])
Y = np.array([])
R = np.array([])
weights = np.array([])
for i_circ in xrange(N_circ):
    X = np.concatenate((X, xx[mask[i_circ]]))
    Y = np.concatenate((Y, yy[mask[i_circ]]))
    R = np.concatenate((R, np.ones(np.sum(mask[i_circ]))*tantwotheta[i_circ]))
    weights = np.concatenate((weights, imData[mask[i_circ]]))
    
x0, y0, r0 = fitconccircle(X, Y, R, weights)
print "x0, y0, r0 =", x0, y0, r0
# draw the image and the circles
fig = plt.figure(figsize=(10,10))
ax1 = fig.add_subplot(2,2,1)
ax1.imshow(imData, origin='lower')
xmin, xmax = ax1.get_xlim()
ymin, ymax = ax1.get_ylim()
for i_circ in xrange(N_circ):
    r = r0*tantwotheta[i_circ]
    cir_x = np.linspace(x0-r+0.5, x0+r-0.5, 200)
    cir_y = y0 - np.sqrt(r**2 - (cir_x-x0)**2)
    ax1.plot(cir_x, cir_y, 'r-')
ax1.set_xlim(xmin, xmax)
ax1.set_ylim(ymin, ymax)

## convert to polar coordinates
ax2 = fig.add_subplot(2,2,2)
# find the extent by looking at the corners
y_corners = np.array([0, 0, ny-1, ny-1])
x_corners = np.array([0, nx-1, 0, nx-1])
rmin = min(np.sqrt((x_corners-x0)**2 + (y_corners-y0)**2))
rmax = max(np.sqrt((x_corners-x0)**2 + (y_corners-y0)**2))
nr = 500 
rr = np.linspace(rmin, rmax, nr, endpoint=False) # r coordinates
tmin = min(np.arctan((y_corners-y0) / (x_corners-x0)))
tmax = max(np.arctan((y_corners-y0) / (x_corners-x0)))
nt = 500
tt = np.linspace(tmin, tmax, nt, endpoint=False) # t coordinates

imData_polar = np.zeros((nt, nr))
for y in range(ny):
    for x in range(nx):
        r = np.sqrt((x-x0)**2 + (y-y0)**2)
        t = np.arctan((y-y0) / (x-x0))
        imData_polar[findind(t,tt), \
                     findind(r,rr)] = imData[y,x]

ax2.imshow(imData_polar, \
        extent=(rmin, rmax, tmin/np.pi*180, tmax/np.pi*180),\
        aspect='auto', origin='lower')
ax2.set_xlabel(r"$r$")
ax2.set_ylabel(r"$\phi$ (deg)")

ax4 = fig.add_subplot(2,2,4, sharex=ax2)
fr = np.sum(imData_polar, axis=0)
ax4.plot(rr, fr, 'b-')
ax4.set_xlabel(r"$r$")

for i in xrange(N_circ):
    ax4.axvline(x=tantwotheta[i]*r0, ls='--')

# plot against 2theta
ax3 = fig.add_subplot(2,2,3)
ax3.plot(np.arctan(rr/r0)/np.pi*180, fr, 'b-')
for i in xrange(N_circ):
    ax3.axvline(x=twotheta_deg[i], ls='--')
ax3.set_xlabel(r"$2\theta$ (deg)")

plt.show()

## save data to file
val = raw_input("Save data to file? ")
if val in ['y', 'Y', 'yes', 'Yes', 's', 'S']:
    params = {'pad': pad, 'x0': x0, 'y0': y0, 'D': r0, \
              'r_array': rr, 't_array': tt}
    data = pickle.load(open("calibration", 'rb'))
    data[pad] = params
    pickle.dump(params, open("calibration.pickle", 'wb'))
