import sys, os
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import PIL.Image
from fitcircle import fitcircle

def findind(x, xmin=0., xmax=0., nx=1):
    dx = (xmax-xmin) * 1. / nx
    xx = int((x-xmin) / dx)
    if xx < 0:
        return 0
    elif xx >= nx:
        return nx-1
    else:
        return int((x-xmin) / dx)


# Read data file into an array
pad = 0 
run = 343
img = 0 
if len(sys.argv) >= 2:
    run = int(sys.argv[1])
if len(sys.argv) >= 3:
    img = int(sys.argv[2])
imFile = PIL.Image.open(
    "LA61_MgO_analysis/MgO_CSPAD/image_cspad%02d_r%04d_e%05d.tif" 
    % (pad, run, img))
run_bkgd = 334
imFile_bkgd = PIL.Image.open(
    "LA61_MgO_analysis/MgO_CSPAD/image_cspad%02d_r%04d_e%05d.tif" 
    % (pad, run_bkgd, 0)) # assume this to be the background
nx, ny = imFile.size
imData = np.array(imFile.getdata()).reshape((ny, nx)).astype(float)
imData_bkgd = np.array(imFile_bkgd.getdata()).reshape((ny, nx)).astype(float)
imData = imData - imData_bkgd
#imData[imData < 0.] = 0.

# Create a figure
fig = plt.figure()

# Create the axes to plot data
ax = fig.add_subplot(2,2,1)

# Display the image
ax.imshow(imData, origin='lower')
xmin, xmax = ax.get_xlim()
ymin, ymax = ax.get_ylim()
pads_x = [0, 191, 200, 385] # these two arrays define the position of the four pads
pads_y = [0, 184, 212, 395]
ax.vlines(pads_x, ymin, ymax, colors='r', linestyles='--')
ax.hlines(pads_y, xmin, xmax, colors='r', linestyles='--')
ax.set_xlim(xmin, xmax)
ax.set_ylim(ymin, ymax)

#ax = fig.add_subplot(2,2,2)
#xx, yy = np.meshgrid(np.arange(nx), np.arange(ny))
## mask for data points to be fitted
#mask = np.logical_and.reduce((-2./3*xx+yy >= 280, -2./3*xx+yy <= 320, \
#            imData > 0.1 * np.max(imData)))
##from pprint import pprint
##pprint(zip(xx[mask], yy[mask], imData[mask]))
#
## delineate approximate circle region
#ax.imshow(mask, origin='lower')
#ax.plot([xmin, xmax], [2./3*xmin+280, 2./3*xmax+280], 'r--') 
#    # line1: y = (2./3) * x + 280
#ax.plot([xmin, xmax], [2./3*xmin+320, 2./3*xmax+320], 'r--') 
#    # line2: y = (2./3) * x + 320
#
### fit to circle, weight by intensity
#x0, y0, r0 = fitcircle(xx[mask], yy[mask], [1.]*np.sum(mask), -150., 700., 420.)
#print "x0, y0, r0 =", x0, y0, r0
## draw the circle
#cir_x = np.linspace(x0-r0, x0+r0, 500)
#cir_y = y0 - np.sqrt(r0**2 - (cir_x-x0)**2)
#ax.plot(cir_x, cir_y, 'r-')
#
#ax.set_xlim(xmin, xmax)
#ax.set_ylim(ymin, ymax)

## read calibration data
import pickle
params = pickle.load(open("calibration_0.pickle", 'rb'))
x0 = params['x0']
y0 = params['y0']
D = params['D']
intercept = params['intercept']
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
onPad = np.zeros((nt, nr))
for y in range(ny):
    for x in range(nx):
        if  (((pads_x[0]<x and x<pads_x[1]) or (pads_x[2]<x and x<pads_x[3]))
         and ((pads_y[0]<y and y<pads_y[1]) or (pads_y[2]<y and y<pads_y[3]))):
            r = np.sqrt((x-x0)**2 + (y-y0)**2)
            t = np.arctan((y-y0) / (x-x0))
            imData_polar[findind(t,tmin,tmax,nt), \
                         findind(r,rmin,rmax,nr)] = imData[y,x]
            onPad[findind(t,tmin,tmax,nt), \
                  findind(r,rmin,rmax,nr)] = 1 

ax2.imshow(imData_polar, \
        extent=(rmin, rmax, tmin/np.pi*180, tmax/np.pi*180),\
        aspect='auto', origin='lower')
ax2.set_xlabel(r"$r$")
ax2.set_ylabel(r"$\phi$ (deg)")

ax = fig.add_subplot(2,2,3)
fr = np.sum(imData_polar, axis=0)/np.sum(onPad, axis=0)*2*np.pi*rr*(rr[1]-rr[0])
twotheta = np.arctan((rr-intercept)/D)
twotheta_deg = twotheta/np.pi*180

ax.plot(twotheta_deg, fr, 'b-')
ax.set_xlabel(r"$2\theta$")

ax4 = fig.add_subplot(2,2,4, sharex=ax2)
ax4.plot(rr, fr, 'b-')
ax4.set_xlabel(r"$r$")
plt.show()

# save data
data = {'pad': pad, 'run': run, 'img': img, 'run_bkgd': run_bkgd,
        'twotheta': twotheta, 'fr': fr}
pickle.dump(data, open("r%3d_p%1d_i%1d_rb%3d.pickle" % (run, pad, img, run_bkgd), 'wb'))

## find the peaks
#rpeaks_ind = np.array([100+np.argmax(fr[100:130]), \
#                       130+np.argmax(fr[130:150]), \
#                       150+np.argmax(fr[150:180]), \
#                       290+np.argmax(fr[290:320])])
#rpeaks = rr[rpeaks_ind]
#ymin, ymax = ax.get_ylim()
#ax.vlines(rpeaks_ind, ymin, ymax, linestyles='dashed')
#ax.set_ylim(ymin, ymax)
