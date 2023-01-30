###############################################################################
# SWIT v1.1: Seismic Waveform Inversion Toolbox
#
# Description: This script generates the acquisition files for the Marmousi2-marine example.
#              The acquisition files include:
#              1) true model (vp_true.npy, rho_true.npy) 
#              2) initial model (vp_init.npy, rho_init.npy)
#              3) gradient mask (grad_mask.npy)
#              4) source coordinates (src_coord.npy)
#              5) source wavelet (wavelets.npy)
#              6) receiver coordinates (rec_coord)
#
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

## set the model size
nx = 481
nz = 141
dx = 25
nt = 4001
dt = 0.002
x_beg = 0.
x_end = (nx-1) * dx

## set true models for generating obs data
vp_true = np.loadtxt('./acquisition/Marmousi_481_141_25m.dat')
rho_true = np.power(vp_true, 0.25) * 310               # density models, (Gardner, 1974)

## set initial model for FWI
vp_init = np.copy(vp_true)                             # copy true model first

# vp_init[:,20:] = smooth2d(vp_init[:,20:] , span = 10)  # smooth the true model, but exclude the water columns (first 20 grids)
vp_init[:,20:] = np.loadtxt('./acquisition/Marmousi_481_121_25m_1d.dat') * 0.75  # smooth the true model, but exclude the water columns (first 20 grids)
rho_init = np.power(vp_init, 0.25) * 310               # density models, (Gardner, 1974)

## set gradient mask for FWI (1: keep gradient; 0: mute gradient)
grad_mask = np.ones((nx, nz))
grad_mask[vp_true == 1500.0] = 0.0

## set source coordinates
src_num = 21
src_coord = np.zeros((src_num, 2))
src_coord[:,0] = np.linspace( x_beg, x_end, src_num)   # linearly distributed sources
src_coord[:,1] = np.linspace(    dx,    dx, src_num)   # sources are buried at first grid depth

## set source wavelet
amp0 = 1.
f0   = 5.
wavelet = np.zeros((src_num, nt))
for i in range(src_num):
    wavelet[i,:] = source_wavelet(amp0, nt, dt, f0, 'ricker')
    # one can load their own wavelet here, e.g., wavelet[i,:] = np.loadtxt('wavelet_src1.dat')
    # the time axis should be the same as the forward modeling

## set receiver coordinates. rec_coord is a list of receiver coordinates for each source 
rec_coord = []
for i in range(src_num):
    rec_xz = np.zeros((nx, 2))
    rec_xz[:,0] = np.linspace(x_beg, x_end, nx)  # linearly distributed receivers
    rec_xz[:,1] = np.linspace(   dx,    dx, nx)  # receivers are buried at first grid depth
    rec_coord.append(rec_xz)                     # receivers can be different for different sources


## save acquisition files
np.save('./acquisition/vp_true.npy',       vp_true)
np.save('./acquisition/rho_true.npy',     rho_true)
np.save('./acquisition/vp_init.npy',       vp_init)
np.save('./acquisition/rho_init.npy',     rho_init)
np.save('./acquisition/grad_mask.npy',   grad_mask)
np.save('./acquisition/src_coord.npy',   src_coord)
np.save('./acquisition/wavelets.npy',      wavelet)
np.savez('./acquisition/rec_coord.npz', *rec_coord) # save receiver list using **np.savez**


## check if the files are saved successfully
files = ['./acquisition/vp_true.npy', './acquisition/rho_true.npy', 
         './acquisition/vp_init.npy', './acquisition/rho_init.npy', 
         './acquisition/grad_mask.npy', './acquisition/src_coord.npy', 
         './acquisition/wavelets.npy', './acquisition/rec_coord.npz']

for file in files:
    if not os.path.exists(file):
        print('Failed to save acquisition files!')
    else:
        print('Successfully save: {}'.format(file))

print('Successfully save acquisition files!\n')
print('You can now run the SWIT workflow:\n')
print('    $ python ../../toolbox-dev/SWIT.py config.yaml\n')
