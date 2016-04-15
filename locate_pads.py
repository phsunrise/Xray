import sys, os
from getopt import getopt
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import RectangleSelector
from matplotlib.patches import Rectangle
import pickle
from helpers import get_data

def sortpads(pads):
## this funciton sorts the "pads" array to meet the format
    pads1 = np.copy(pads)
    pads2 = np.copy(pads)

    # first sort y direction
    args = np.argsort(pads1[:, 1])
    pads2 = np.array([pads1[args[0]], pads1[args[1]],
                      pads1[args[2]], pads1[args[3]]])
    pads1 = np.copy(pads2)
    
    # now sort x direction
    if pads1[0,0] > pads1[1,0]:
        pads2[0], pads2[1] = pads1[1], pads1[0]
    if pads1[2,0] > pads1[3,0]:
        pads2[2], pads2[3] = pads1[3], pads1[2]
    pads1 = np.copy(pads2)

    return pads1


# Read data file into an array
pad = 2
run = 389
img = 1 
# parse command line arguments 
opts, args = getopt(sys.argv[1:], "p:r:i:")
for opt, arg in opts:
    if opt == '-p':
        pad = int(arg)
    elif opt == '-r':
        run = int(arg)
    elif opt == '-i':
        img = int(arg)
 
imData_orig = get_data(pad, run, img) 
# enhance contrast at low values
immin = np.min(imData_orig)
immax = np.max(imData_orig)
A = 20 # curve factor
val = ''

while True:
    pads = np.zeros((4,4))
    pads_format = "Format: numpy 4*4 array, each line is (x_LL, y_LL, x_UR, y_UR), and the pads are LL, LR, UL, UR"
    i_pad = 0

    # Transform imData so that immin corresponds to 1 and immax to 100;
    # immin is at the center of a gaussian dip
    imData = 100. - 99.*np.exp(-A*((imData_orig-immin)/(immax-immin))**2)

    # Exponential curve
    #imData = ((immax-immin*np.exp(-A))/(1-np.exp(-A))
    #         -(immax-immin)/(1-np.exp(-A))*np.exp(-(imData_orig-immin)/(immax-immin)*A))

    fig = plt.figure(figsize=(10,10))
    ax = fig.add_subplot(1,1,1)
    ax.imshow(imData, origin='lower')

    def onselect(eclick, erelease):
        global i_pad
        if i_pad >= 4:
            print "All four pads have been chosen! Closing plot..."
            plt.close()
        # eclick and erelease are the press and release events
        x1, y1 = eclick.xdata, eclick.ydata
        x2, y2 = erelease.xdata, erelease.ydata
        print("Pad %d selected: (%3.2f, %3.2f) --> (%3.2f, %3.2f)" % (i_pad, x1, y1, x2, y2))
        x1, x2 = min(x1, x2), max(x1, x2)
        y1, y2 = min(y1, y2), max(y1, y2)
        pads[i_pad] = x1, y1, x2, y2
        ax.add_patch(Rectangle((x1, y1), x2-x1, y2-y1, \
                               color='g', fill=True, alpha=0.5))
        i_pad += 1
        if i_pad == 4:
            print "All four pads have been chosen! Close the window to exit"

    rect = RectangleSelector(ax, onselect,
                drawtype='box', minspanx=5, minspany=5)
    plt.show()

    val = raw_input("Enter contrast factor, \'s\' to save, or \'q\' to exit: ")
    if val in ['q', 'Q', 'quit', 'exit', 'Exit']:
        break
    elif val in ['s', 'S', 'save']: # save pad locations
        pads = sortpads(pads) # sort pads to meet format
        data = pickle.load(open("parameters.pickle", 'r'))
        try:
            data[pad]['pads'] = pads
            data[pad]['pads_format'] = pads_format
        except KeyError:
            data[pad] = {'pads': pads, 'pads_format': pads_format}
        pickle.dump(data, open("parameters.pickle", 'w'))
        print "Pad positions saved!"
        break
    else:
        try:
            A = int(val)
        except ValueError:
            print "Wrong value!"
            val = raw_input("Enter contrast factor, \'s\' to save, or \'q\' to exit: ")
            continue 
