import os, sys
sys.path.insert(1, os.path.join(sys.path[0], '..'))

import numpy as np
import matplotlib.pyplot as plt
from adjust_subpads import adjust_subpads
from helpers import get_data

pad = 2
run = 304
img = 1
run_bkgd = 401 
img_bkgd = 0
do_debug = False
fixpad = 0
figname = "test"

imData = get_data(pad, run, img)
imData_bkgd = get_data(pad, run_bkgd, img_bkgd)
(imData_polar, onPad, rr, tt, \
     twotheta_deg1, fr1, offset) = adjust_subpads(imData, pad, \
                                   fixpad=fixpad, do_debug=do_debug, \
                                   figname=figname)
(imData_polar, onPad, rr, tt, \
     twotheta_deg2, fr2, offset) = adjust_subpads(imData_bkgd, pad, \
                                   fixpad=fixpad, do_debug=do_debug, \
                                   figname=figname)
(imData_polar, onPad, rr, tt, \
     twotheta_deg3, fr3, offset) = adjust_subpads(imData-imData_bkgd, pad, \
                                   fixpad=fixpad, do_debug=do_debug, \
                                   figname=figname)

plt.plot(twotheta_deg1, fr1, 'b-', label="orig")
plt.plot(twotheta_deg2, fr2, 'g-', label="bkgd")
plt.plot(twotheta_deg3, fr3, 'r-', label="orig-bkgd")
plt.legend()
plt.xlabel(r"$2\theta$")
plt.title("Shot %d" % run)
plt.savefig("bkgdSubtract_p%d_r%04d.pdf" % (pad, run))
plt.show()
