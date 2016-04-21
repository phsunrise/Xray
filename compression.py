import sys
from getopt import getopt
import numpy as np
import pickle
import matplotlib.pyplot as plt
from scipy.optimize import fsolve
from helpers import get_delaytime

def twotheta_comp(twotheta_deg, ratio):
    theta = twotheta_deg / 2. / 180. * np.pi
    theta_comp = np.arcsin(np.sin(theta) * ratio**(1./3)) 
    return 2 * theta_comp / np.pi * 180.

coldpeaks = np.array([37.3, 43.3])
coldpeaks_label = ['(111)', '(200)']
xmin, xmax = 1., 3.
ymin, ymax = 35., 65.
ratio = np.linspace(xmin, xmax, 200)
comppeaks = []
for i in xrange(len(coldpeaks)):
    comppeaks.append(twotheta_comp(coldpeaks[i], ratio))

do_debug = False
opts, args = getopt(sys.argv[1:], "d")
for opt, arg in opts:
    if opt == '-d':
        do_debug = True
     
data0 = pickle.load(open("data_0.pickle", 'rb'))
data2 = pickle.load(open("data_2.pickle", 'rb'))
for run in data2.keys():
    try:
        run = int(run)
    except ValueError:
        continue

    fig = plt.figure()
    ax1 = fig.add_subplot(1,1,1)
    intersects = [] # list to store intersects; used for xticks
    # plot measured peaks and solve for intersections
    params2 = data2[run]['params']
    for i in range(3, len(params2), 3):
        mean = params2[i+1]
        ax1.axhline(y=mean, c='k', ls='-')
        # now calculate intersection with each cold peak
        for coldpeak in coldpeaks:
            func = lambda x: twotheta_comp(coldpeak, x)-mean
            intersect = fsolve(func, 1.)
            intersects.append(intersect)
            ax1.vlines(intersect, ymin, mean, \
                       colors='k', linestyles='--')

    try:
        params0 = data0[run]['params']
        for i in range(3, len(params0), 3):
            mean = params0[i+1]
            ax1.axhline(y=mean, c='b', ls='-')
            # now calculate intersection with each cold peak
            for coldpeak in coldpeaks:
                func = lambda x: twotheta_comp(coldpeak, x)-mean
                intersect = fsolve(func, 1.)
                intersects.append(intersect)
                ax1.vlines(intersect, ymin, mean, \
                           colors='b', linestyles='--')
    except KeyError:
        pass

    ax1.set_xlabel("compression ratio")
    ax1.set_ylabel(r"$2\theta$ (deg)", color='k')

    # plot compressed cold peaks
    for i in xrange(len(coldpeaks)):
        ax1.plot(ratio, comppeaks[i], 'r-')

    # add xticks and peak labels on right 
    ax2 = ax1.twinx()
    # add xticks
    intersects.sort()
    intersects_label = ['']*len(intersects)
    for i in xrange(len(intersects)):
        if abs(intersects[i] - intersects[i-1]) < 0.02:
            pass 
        else:
            intersects_label[i] = '%.2f' % intersects[i]

    ax1.set_xticks(intersects)
    ax1.set_xticklabels(intersects_label, \
                        rotation='vertical')
    ax2.set_yticks(twotheta_comp(coldpeaks, xmax))
    ax2.set_yticklabels(coldpeaks_label, color='r')
    ax1.set_xlim(xmin, xmax)
    ax2.set_ylim(ymin, ymax)

    t = get_delaytime(run) 
    ax1.set_title("Shot %d, delay %.1f ns" % (run, t))
    plt.tight_layout()
    fig.savefig("t%.1f_r%04d_compression.pdf" % (t, run))

    if do_debug:
        plt.show()
        sys.exit(0)

    plt.close()
