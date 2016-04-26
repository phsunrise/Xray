import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np
import pickle

pad = 0
data = pickle.load(open("XRonly_p%d_data.pickle" % pad, 'rb'))
tt = data['tt']
tt_deg = tt / np.pi * 180.

with PdfPages("XRonly/all_XRonly.pdf") as pdf:
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
            if np.std(ft_filtered) < 20.:
                ax1.plot(tt_deg, ft+300.*img, c='k')
                continue # skipping potentitally no x-ray runs
            ax1.plot(tt_deg, ft+300.*img, c=cm.jet(img*1./Nimg))
            ft_mean[img] = np.mean(ft_filtered)
            ft_std[img] = np.std(ft_filtered)
            ft_range[img] = sum(ft_filtered>500)

        imgs = np.arange(Nimg)
        mask = (ft_std > 0)
        ax2.plot(imgs[mask], ft_mean[mask], 'ro--')
        ax3.plot(imgs[mask], ft_std[mask], 'ro--')
        ax4.plot(imgs[mask], ft_range[mask], 'ro--')
        mask = np.logical_not(mask)
        ax2.plot(imgs[mask], ft_mean[mask], 'bo')
        ax3.plot(imgs[mask], ft_std[mask], 'bo')
        ax4.plot(imgs[mask], ft_range[mask], 'bo')
        
        ax1.set_xlabel(r"$\phi$ (deg)")
        ax1.set_title("Shot %d" % run)

        ax2.set_xticks([])
        ax2.set_title("Mean")

        ax3.set_xticks([])
        ax3.set_title("Standard deviation")

        ax4.set_xticks([])
        ax4.set_yticks([])
        ax4.set_title(">500 range")

        fig.savefig("XRonly/r%04d_p%1d_XRonly.pdf" % (run, pad))
        pdf.savefig(fig)
        plt.close()

    ## END run loop    
## END with PdfPages
