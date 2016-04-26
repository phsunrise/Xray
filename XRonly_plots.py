import matplotlib.pyplot as plt
from matplotlib import cm
import numpy as np
import pickle

pad = 0
data = pickle.load(open("XRonly_p%d_data.pickle" % pad, 'rb'))
tt = data['tt']
tt_deg = tt / np.pi * 180.

for run in data.keys():
    if not isinstance(run, int):
        continue
    fig = plt.figure(figsize=(15, 15))
    ax1 = plt.subplot2grid((2,3), (0,0), colspan=3) # to plot data
    ax2 = plt.subplot2grid((2,3), (1,0)) # to plot mean
    ax3 = plt.subplot2grid((2,3), (1,1)) # to plot rms 
    ax4 = plt.subplot2grid((2,3), (1,2)) # to plot range

    thisrun = data[run]
    Nimg = len(thisrun.keys()) # for colormap
    ft_mean = np.zeros(Nimg)
    ft_std = np.zeros(Nimg)
    ft_range = np.zeros(Nimg)
    for img in thisrun.keys():
        ft = thisrun[img]
        ft_filtered = ft[np.isfinite(ft)]
        ax1.plot(tt_deg, ft+300.*img, c=cm.jet(img*1./Nimg))
        ft_mean[img] = np.mean(ft_filtered)
        ft_std[img] = np.std(ft_filtered)
        ft_range = sum(ft_filtered>500)

    ax2.plot(ft_mean, 'ro--')
    ax3.plot(ft_std, 'ro--')
    ax4.plot(ft_range, 'ro--')
    
    ax1.set_xlabel(r"$\phi$ (deg)")
    ax1.set_title("Shot %d" % run)

    ax2.set_xticks([])
    ax2.set_title("Mean")

    ax3.set_xticks([])
    ax3.set_title("Standard deviation")

    ax4.set_xticks([])
    ax4.set_yticks([])
    ax4.set_title("Range of >500 points")

    fig.savefig("r%04d_p%1d_XRonly.pdf" % (run, pad))
    plt.close()

## END run loop    
