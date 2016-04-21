import sys
from getopt import getopt
import numpy as np
import pickle
import matplotlib.pyplot as plt

coldpeaks = np.array([37.3, 43.3])
ratio = np.linspace(1., 3., 200)
comppeaks = []
for i in xrange(len(coldpeaks)):
    comppeaks.append(np.arctan( np.tan(coldpeaks[i]/180.*np.pi)
                              * (ratio**(1./3)) )
                    /np.pi*180.)

opts, args = getopt(sys.argv[1:], "p:r:")
for opt, arg in opts:
    if opt == '-p':
        pad = int(arg)
    #elif opt == '-r':
    #    run = int(arg)

data = pickle.load(open("data_%d.pickle" % pad, 'rb'))
for run in data.keys():
    try:
        run = int(run)
    except ValueError:
        continue

    # plot measured peaks
    params = data[run]['params']
    for i in range(3, len(params), 3):
        plt.axhline(y=params[i+1], c='k', ls='-')

    # plot compressed cold peaks
    for i in xrange(len(coldpeaks)):
        plt.plot(ratio, comppeaks[i], 'r-')

    plt.xlabel("compression ratio")
    plt.ylabel(r"$2\theta$ (deg)")
    plt.title("Shot %d" % run)
    plt.savefig("r%04d_p%d_compression.pdf" % (run, pad))
