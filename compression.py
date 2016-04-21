import sys
from getopt import getopt
import numpy as np
import pickle
import matplotlib.pyplot as plt
from helpers import get_delaytime

coldpeaks = np.array([37.3, 43.3])
coldpeaks_label = ['(111)', '(200)']
ratio = np.linspace(1., 3., 200)
comppeaks = []
for i in xrange(len(coldpeaks)):
    comppeaks.append(2*np.arcsin( np.sin(coldpeaks[i]/2./180.*np.pi)
                              * (ratio**(1./3)) )
                    /np.pi*180.)

#opts, args = getopt(sys.argv[1:], "p:")
#for opt, arg in opts:
#    if opt == '-p':
#        pad = int(arg)

data0 = pickle.load(open("data_0.pickle", 'rb'))
data2 = pickle.load(open("data_2.pickle", 'rb'))
for run in data2.keys():
    try:
        run = int(run)
    except ValueError:
        continue

    fig = plt.figure()
    ax1 = fig.add_subplot(1,1,1)
    # plot measured peaks
    params2 = data2[run]['params']
    for i in range(3, len(params2), 3):
        ax1.axhline(y=params2[i+1], c='k', ls='-')

    try:
        params0 = data0[run]['params']
        for i in range(3, len(params0), 3):
            ax1.axhline(y=params0[i+1], c='b', ls='-')
    except KeyError:
        pass

    ax1.set_xlabel("compression ratio")
    ax1.set_ylabel(r"$2\theta$ (deg)", color='k')

    # plot compressed cold peaks
    for i in xrange(len(coldpeaks)):
        ax1.plot(ratio, comppeaks[i], 'r-')

    # add labels on right 
    ax2 = ax1.twinx()
    xmin, xmax = ax1.get_xlim()
    ymin, ymax = ax1.get_ylim()
    ax2.set_ylim(ymin, ymax)
    ax2.set_yticks(2*np.arcsin(np.sin(coldpeaks/2./180.*np.pi)
                  * (xmax**(1./3))) / np.pi * 180.)
    ax2.set_yticklabels(coldpeaks_label, color='r')

    t = get_delaytime(run) 
    ax1.set_title("Shot %d, delay %.1f ns" % (run, t))
    fig.savefig("t%.1f_r%04d_compression.pdf" % (t, run))
    plt.close()
