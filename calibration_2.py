import sys, os
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import PIL.Image
from fitcircle import fitcircle
from roipoly import roipoly
from process_2 import findind

# Read data file into an array
pad = 2
run = 389
img = 1 
imFile = PIL.Image.open(
    "LA61_MgO_analysis/MgO_CSPAD/image_cspad%02d_r%04d_e%05d.tif" 
    % (pad, run, img))
imFile_bkgd = PIL.Image.open(
    "LA61_MgO_analysis/MgO_CSPAD/image_cspad%02d_r%04d_e%05d.tif" 
    % (pad, 401, 0)) # assume this to be the background
nx, ny = imFile.size
imData = np.array(imFile.getdata()).reshape((ny, nx)).astype(float)
imData_bkgd = np.array(imFile_bkgd.getdata()).reshape((ny, nx)).astype(float)
imData = imData - imData_bkgd

# pyplot setup 
fig = plt.figure()
ax = fig.add_subplot(1,1,1)

# Display the image
ax.imshow(imData, origin='lower')
xmin, xmax = ax.get_xlim()
ymin, ymax = ax.get_ylim()

# use ROI to select region
ROI = roipoly(roicolor='r')
ax.imshow

xx, yy = np.meshgrid(np.arange(nx), np.arange(ny))
# mask for data points to be fitted
mask = np.logical_and.reduce((yy <= k1*xx+b1, yy >= k2*xx+b2, \
            imData > 0.1 * np.max(imData)))
#from pprint import pprint
#pprint(zip(xx[mask], yy[mask], imData[mask]))

# delineate approximate circle region
ax.imshow(mask, origin='lower')
ax.plot([xmin, xmax], [k1*xmin+b1, k1*xmax+b1], 'r--') 
ax.plot([xmin, xmax], [k2*xmin+b2, k2*xmax+b2], 'r--') 

## fit to circle, weight by intensity
x0, y0, r0 = fitcircle(xx[mask], yy[mask], [1.]*np.sum(mask), -150., 700., 480.)
print "x0, y0, r0 =", x0, y0, r0
# draw the circle
cir_x = np.linspace(x0-r0, x0+r0, 500)
cir_y = y0 - np.sqrt(r0**2 - (cir_x-x0)**2)
ax.plot(cir_x, cir_y, 'r-')

ax.set_xlim(xmin, xmax)
ax.set_ylim(ymin, ymax)

## convert to polar coordinates
ax = fig.add_subplot(2,3,3)
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

ax.imshow(imData_polar, \
        extent=(rmin, rmax, tmin/np.pi*180, tmax/np.pi*180),\
        aspect='auto', origin='lower')
ax.set_xlabel(r"$r$")
ax.set_ylabel(r"$\phi$ (deg)")

ax = fig.add_subplot(2,3,4)
fr = np.sum(imData_polar, axis=0)
ax.plot(fr, 'b-')
ax.set_xlabel("r index")

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
import pickle
params = {'pad': pad, 'x0': x0, 'y0': y0, 'D': D, 'intercept': intercept}
pickle.dump(params, open("calibration_2.pickle", 'wb'))
