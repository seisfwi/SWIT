###############################################################################
# SWIT v1.1: Seismic Waveform Inversion Toolbox
#
# Description: This script generates the acquisition files for the Marmousi2-marine example.
#              The acquisition files include:
#              1) true model (vp_true.npy, rho_true.npy) 
#              2) initial model (vp_init.npy, rho_init.npy)
#              3) source coordinates (src_coord.npy)
#              4) source wavelet (wavelets.npy)
#              5) receiver coordinates (rec_coord)
#    Note: The gradient mask is not needed for land acquisition.
#           The gradient damping is implemented in the code by default.
###############################################################################


## add toolbox path
import sys
sys.path.append('/homes/sep/haipeng/develop/SWIT-1.0/dev/toolbox-dev/')

import os
import numpy as np
from utils import source_wavelet
from tools import smooth2d

## create acquisition folder
if not os.path.exists('acquisition'):
    os.makedirs('acquisition', exist_ok=True)

## set model size
nx = 481
nz = 121
dx = 25
nt = 4001
dt = 0.002
x_beg = 0.
x_end = (nx-1) * dx

## save stage 1 results as stage 2 initial model
path = "/scr2/haipeng/SWIT-1.1/07-WT-tomo/stage1/fwi/model/vp_it_0010.npy"
vp_init = np.load(path)
# vp_init = smooth2d(vp_init, 5)                         # smooth tomo model a little bit
rho_init = np.power(vp_init, 0.25) * 310               # density models, (Gardner, 1974)

np.save('./acquisition/vp_tomo.npy',       vp_init)
np.save('./acquisition/rho_tomo.npy',     rho_init)