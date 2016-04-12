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
pad = 2
run = 373 
img = 0
run_bkgd = 401 # dark current run
bkgdSubtract = True 
if len(sys.argv) >= 2:
    run = int(sys.argv[1])
if len(sys.argv) >= 3:
    img = int(sys.argv[2])
if len(sys.argv) >= 4:
    run_bkgd = int(sys.argv[3])
    bkgdSubtract = True

imFile = PIL.Image.open(
    "LA61_MgO_analysis/MgO_CSPAD/image_cspad%02d_r%04d_e%05d.tif" 
    % (pad, run, img))
imFile_bkgd = PIL.Image.open(
    "LA61_MgO_analysis/MgO_CSPAD/image_cspad%02d_r%04d_e%05d.tif" 
    % (pad, run_bkgd, 0)) # assume this to be the background
nx, ny = imFile.size
imData = np.array(imFile.getdata()).reshape((ny, nx)).astype(float)
imData_bkgd = np.array(imFile_bkgd.getdata()).reshape((ny, nx)).astype(float)

pads_x = [12, 190, 197, 385] # these two arrays define the position of the four pads
pads_y = [0, 184, 212, 395]
## printing each subpad
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

# subtract average background
#imData = imData - np.mean(imData[0:150, 0:150])
#imData[imData < 0.] = 0.

# Create a figure
fig = plt.figure()

# Create the axes to plot data
ax1 = fig.add_subplot(2,3,1)

## Display the image
ax1.imshow(imData, origin='lower')
xmin, xmax = ax1.get_xlim()
ymin, ymax = ax1.get_ylim()
ax1.vlines(pads_x, ymin, ymax, colors='r', linestyles='--')
ax1.hlines(pads_y, xmin, xmax, colors='r', linestyles='--')
ax1.set_xlim(xmin, xmax)
ax1.set_ylim(ymin, ymax)

## read calibration data
import pickle
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
                         findind(t,tmin,tmax,nt), \
                         findind(r,rmin,rmax,nr)] = imData[y,x]
            onPad[pad,
                  findind(t,tmin,tmax,nt), \
                  findind(r,rmin,rmax,nr)] = 1 
ax2 = fig.add_subplot(2,3,2)
ax2.imshow(np.sum(imData_polar, axis=0), \
        extent=(rmin, rmax, tmin/np.pi*180, tmax/np.pi*180),\
        aspect='auto', origin='lower')

## print out pad ranges as index in rr
print "pad0:"
print findind(np.sqrt((pads_x[0]-x0)**2+(pads_y[1]-y0)**2), rmin, rmax, nr)
print findind(np.sqrt((pads_x[1]-x0)**2+(pads_y[0]-y0)**2), rmin, rmax, nr)
print "pad1:"
print findind(np.sqrt((pads_x[2]-x0)**2+(pads_y[1]-y0)**2), rmin, rmax, nr)
print findind(np.sqrt((pads_x[3]-x0)**2+(pads_y[0]-y0)**2), rmin, rmax, nr)
print "pad2:"
print findind(np.sqrt((pads_x[0]-x0)**2+(pads_y[3]-y0)**2), rmin, rmax, nr)
print findind(np.sqrt((pads_x[1]-x0)**2+(pads_y[2]-y0)**2), rmin, rmax, nr)
print "pad3:"
print findind(np.sqrt((pads_x[2]-x0)**2+(pads_y[3]-y0)**2), rmin, rmax, nr)
print findind(np.sqrt((pads_x[3]-x0)**2+(pads_y[2]-y0)**2), rmin, rmax, nr)


ax4 = fig.add_subplot(2,3,4)
twotheta = np.arctan((rr-intercept)/D)
twotheta_deg = twotheta/np.pi*180

fr = np.zeros((4, nr))
for i_pad in xrange(4):
    fr[i_pad] = np.sum(imData_polar[i_pad], axis=0)/np.sum(onPad[i_pad], axis=0)
    ax4.plot(twotheta_deg, fr[i_pad], label='%d'%i_pad)
ax4.set_xlabel(r"$2\theta$")
ax4.legend()

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
ax5 = fig.add_subplot(2,3,5)
for i_pad in xrange(4):
    ax5.plot(twotheta_deg, fr[i_pad], label='%d'%i_pad)
ax5.set_xlabel(r"$2\theta$")
ax5.legend()

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
            imData_polar[findind(t,tmin,tmax,nt), \
                         findind(r,rmin,rmax,nr)] = imData[y,x]
            onPad[findind(t,tmin,tmax,nt), \
                  findind(r,rmin,rmax,nr)] = 1 
ax3 = fig.add_subplot(2,3,3)
ax3.imshow(imData_polar, \
        extent=(rmin, rmax, tmin/np.pi*180, tmax/np.pi*180),\
        aspect='auto', origin='lower')
ax3.set_xlabel(r"$r$")
ax3.set_ylabel(r"$\phi$ (deg)")

fr = np.sum(imData_polar, axis=0)/np.sum(onPad, axis=0)
ax6 = fig.add_subplot(2,3,6, sharex=ax3)
#ax6.semilogy(rr, fr, 'ro', label='orig.')
# moving average
N_ma = 10
window = np.ones(N_ma)*1./N_ma
ax6.semilogy(rr, fr, 'b-', label='orig.')
fr_fit = lambda x,a,b,c: a*np.exp(-b*x)+c
fit_start = 20
fit_end = 400
c0 = 0
b0 = -np.log(fr[fit_end]*1./fr[fit_start])/(rr[fit_end]-rr[fit_start])
a0 = fr[fit_start]/np.exp(-b0*rr[fit_start])
print "a0, b0, c0:", a0, b0, c0
popt, pcov = curve_fit(fr_fit, rr[fit_start:fit_end], fr[fit_start:fit_end], \
                       p0=[a0, b0, c0])
print "fit:", popt, np.sqrt(np.diag(pcov)) 
ax6.semilogy(rr[fit_start:fit_end], \
             popt[0]*np.exp(-popt[1]*rr[fit_start:fit_end])+popt[2], \
             'r-', label='fit')
ax6.set_xlabel(r"$r$")
ax6.legend()
plt.show()

# save data
data = {'pad': pad, 'run': run, 'img': img, 'run_bkgd': run_bkgd,
        'twotheta': twotheta, 'fr': fr}
if bkgdSubtract:
    pickle.dump(data, open("r%3d_p%1d_i%1d_rb%3d.pickle" % (run, pad, img, run_bkgd), 'wb'))
else:
    pickle.dump(data, open("r%3d_p%1d_i%1d_noBkgdSubtract.pickle" % (run, pad, img), 'wb'))


## find the peaks
#rpeaks_ind = np.array([100+np.argmax(fr[100:130]), \
#                       130+np.argmax(fr[130:150]), \
#                       150+np.argmax(fr[150:180]), \
#                       290+np.argmax(fr[290:320])])
#rpeaks = rr[rpeaks_ind]
#ymin, ymax = ax.get_ylim()
#ax.vlines(rpeaks_ind, ymin, ymax, linestyles='dashed')
#ax.set_ylim(ymin, ymax)
