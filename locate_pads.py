import sys, os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import RectangleSelector
from matplotlib import patches
import PIL.Image
import pickle

## initialization
pads = np.zeros((4,4))
pads_format = "Format: numpy 4*4 array, each line is (x_LL, y_LL, x_UR, y_UR),\
               and the pads are LL, LR, UL, UR"
i_pad = 0

def sortpads(pads):
## this funciton sorts the "pads" array to meet the format
    print pads
    pads1 = np.copy(pads)
    pads2 = np.copy(pads)
    # first sort each pad
    for i in xrange(4):
        if pads1[i,0] > pads1[i,2]:
            pads2[i,0], pads2[i,2] = pads1[i,2], pads1[i,0]
        if pads1[i,1] > pads1[i,3]:
            pads2[i,1], pads2[i,3] = pads1[i,3], pads1[i,1]
    pads1 = np.copy(pads2)

    # now sort y direction
    args = np.argsort(pads1[:, 1])
    pads2 = np.array([pads1[args[0]], pads1[args[1]],
                      pads1[args[2]], pads1[args[3]]])
    pads1 = np.copy(pads2)
    print pads1
    
    # now sort x direction
    if pads1[0,0] > pads1[1,0]:
        pads2[0], pads2[1] = pads1[1], pads1[0]
    if pads1[2,0] > pads1[3,0]:
        pads2[2], pads2[3] = pads1[3], pads1[2]
    pads1 = np.copy(pads2)
    print pads1

    return pads1

def showpads(pads, ax):
    pads1 = sortpads(pads)
    for i in xrange(4):
        ax.add_patch(patches.Rectangle((pads1[i,0], pads1[i,1]),\
                 pads1[i,2]-pads1[i,0], pads1[i,3]-pads1[i,1],\
                 fill=False))
    plt.show()


def line_select_callback(eclick, erelease):
    # eclick and erelease are the press and release events
    global pads, i_pad
    x1, y1 = eclick.xdata, eclick.ydata
    x2, y2 = erelease.xdata, erelease.ydata
    print("Pad %d selected: (%3.2f, %3.2f) --> (%3.2f, %3.2f)" % (i_pad, x1, y1, x2, y2))
    #print(" The button you used were: %s %s" % (eclick.button, erelease.button))
    pads[i_pad] = x1, y1, x2, y2
    i_pad += 1
    if i_pad >= 4:
        print "All four pads are selected."
        plt.close()
    

def toggle_selector(event):
    global i_pad
    print(' Key pressed.')
    if event.key in ['Q', 'q'] and toggle_selector.RS.active:
        print(' RectangleSelector deactivated.')
        toggle_selector.RS.set_active(False)
    if event.key in ['A', 'a'] and not toggle_selector.RS.active:
        print(' RectangleSelector activated.')
        toggle_selector.RS.set_active(True)
    if event.key in ['D', 'd'] and i_pad >= 1:
        i_pad -= 1


# Read data file into an array
pad = 2
run = 389
img = 1 
imFile = PIL.Image.open(
    "LA61_MgO_analysis/MgO_CSPAD/image_cspad%02d_r%04d_e%05d.tif" 
    % (pad, run, img))
nx, ny = imFile.size
imData_orig = np.array(imFile.getdata()).reshape((ny, nx)).astype(float)
# enhance contrast at low values
immin = np.min(imData_orig)
immax = np.max(imData_orig)
A = 20 # curve factor
val = ''

while val not in ['q', 'Q', 'quit', 'exit', 'Exit']:
    if val in ['s', 'S', 'save']: # save pad locations
        pads = sortpads(pads) # sort pads to meet format
        if os.path.isfile("calibration.pickle"):
            data = pickle.load(open("calibration.pickle", 'rb'))
        else:
            data = {}
        try:
            data[pad]['pads'] = pads
            data[pad]['pads_format'] = pads_format
        except KeyError:
            data[pad] = {'pads': pads, 'pads_format': pads_format}
        pickle.dump(data, open("calibration.pickle", 'wb'))
        break
    elif val != '':
        try:
            A = int(val)
        except ValueError:
            print "ValueError"
            val = raw_input("Enter contrast factor, or \'q\' to exit: ")
            continue 
    
    # Transform imData so that immin corresponds to 1 and immax to 100;
    # immin is at the center of a gaussian dip
    imData = 100. - 99.*np.exp(-A*((imData_orig-immin)/(immax-immin))**2)

    # Exponential curve
    #imData = ((immax-immin*np.exp(-A))/(1-np.exp(-A))
    #         -(immax-immin)/(1-np.exp(-A))*np.exp(-(imData_orig-immin)/(immax-immin)*A))

    fig = plt.figure(figsize=(10,10))
    ax = fig.add_subplot(1,1,1)
    ax.imshow(imData)
    toggle_selector.RS = RectangleSelector(ax, line_select_callback,
                    drawtype='box', useblit=True,
                    button=[1, 3],  # don't use middle button
                    minspanx=5, minspany=5,
                    spancoords='pixels',
                    interactive=True)
    plt.connect('key_press_event', toggle_selector)
    plt.show()

    if i_pad >= 4:
        plt.figure(figsize=(10,10))
        plt.imshow(imData)
        ax = plt.gca()
        showpads(pads, ax)

    val = raw_input("Enter contrast factor, s to save, or \'q\' to exit: ")
