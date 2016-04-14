import sys, os
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import PIL.Image
import pickle

from fitcircle import fitcircle
from roipoly import roipoly
from helpers import findind, get_data, get_pads

# Read data file into an array
pad = 2
run = 389
img = 1 
run_bkgd = 401
imData = get_data(pad, run, img)
ny, nx = imData.shape
imData_bkgd = get_data(pad, run_bkgd, 0) 
imData = imData - imData_bkgd

# Display the image
fig = plt.figure(figsize=(10,10))
ax = fig.add_subplot(1,1,1)
ax.imshow(imData, origin='lower')
ax.set_title("Select polygon region: right or double click to finish")

# use ROI to select region
ROI = roipoly(roicolor='r')
#plt.imshow(imData, origin='lower')
#ROI.displayROI()
mask = np.logical_and(ROI.getMask(imData), imData>0.3*np.max(imData))
#plt.imshow(mask, origin='lower', cmap='Greys')

## fit to circle, weight by intensity
xx, yy = np.meshgrid(np.arange(nx), np.arange(ny))
x0, y0, r0 = fitcircle(xx[mask], yy[mask], imData[mask])
print "x0, y0, r0 =", x0, y0, r0
# draw the image and the circle
cir_x = np.linspace(x0-r0+0.5, x0+r0-0.5, 500)
cir_y = y0 - np.sqrt(r0**2 - (cir_x-x0)**2)

fig = plt.figure()
ax1 = fig.add_subplot(2,2,1)
ax1.imshow(imData, origin='lower')
xmin, xmax = ax1.get_xlim()
ymin, ymax = ax1.get_ylim()

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
ax4.plot(fr, 'b-')
ax4.set_xlabel("r index")
plt.show()

# find the three peaks
rpeaks_ind = np.array([
                100 + np.argmax(fr[100:200]),
                250 + np.argmax(fr[250:350]),
                400 + np.argmax(fr[400:500])
                ])
rpeaks = rr[rpeaks_ind]
ymin, ymax = ax.get_ylim()
ax.vlines(rpeaks_ind, ymin, ymax, linestyles='dashed')
ax.set_ylim(ymin, ymax)

## calibration
twotheta = np.array([47.9, 57.4, 63.1])
tantwotheta = np.tan(twotheta/180.*np.pi)
# use linear regression to find d in r = d*tan(2*theta)
from scipy.stats import linregress
D, intercept, r_val, p_val, std = linregress(tantwotheta, rpeaks)
print "D =", D, "intercept =", intercept
print "r^2 =", r_val**2

# plot calibrated values
ax = fig.add_subplot(2,3,5)

#ax.plot(tantwotheta, rpeaks, 'ro')
#x = np.linspace(np.tan(30*np.pi/180), np.tan(50*np.pi/180), 100)
#y = d * x + intercept
#ax.plot(x, y, 'b-')

ax.plot(np.arctan((rr-intercept)*1./D)/np.pi*180, fr, 'b-')
ax.set_xlabel(r"$2\theta$ (deg)")
# plot vertical lines
ymin, ymax = ax.get_ylim()
ax.vlines(twotheta, ymin, ymax, linestyles='dashed')
ax.set_ylim(ymin, ymax)

plt.show()

## save data to file
val = raw_input("Save data to file?")
if val in ['y', 'Y', 'yes', 'Yes', 's', 'S']:
    params = {'pad': pad, 'x0': x0, 'y0': y0, 'D': D, 'intercept': intercept}
    data = pickle.load(open("calibration", 'rb'))
    data[pad] = params
    pickle.dump(params, open("calibration.pickle", 'wb'))
