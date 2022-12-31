import sys
sys.path.append('/homes/sep/haipeng/develop/SWIT-1.0/toolbox-dev/')

import os
import numpy as np
from utils import source_wavelet
from tools import smooth2d

## create acquisition folder
os.makedirs('acquisition', exist_ok=True)

## set model size
nx = 501
nz = 171
dx = 20
nt = 3001
dt = 0.0015
x_beg = 0.
x_end = (nx-1) * dx

## set true models for generating obs data
vp_true = np.loadtxt('./acquisition/Overthrust_501_171_20m.dat')
rho_true = np.power(vp_true, 0.25) * 310               # density models, (Gardner, 1974)

## set initial model for FWI
vp_init = np.copy(vp_true)                             # copy true model first
vp_init[:,20:] = smooth2d(vp_init[:,20:] , span = 10)  # smooth the true model, but exclude the water columns (first 20 grids)
rho_init = np.power(vp_init, 0.25) * 310               # density models, (Gardner, 1974)

## set gradient mask for FWI (1: keep gradient; 0: mute gradient)
grad_mask = np.ones((nx, nz))
grad_mask[vp_true == 1500.0] = 0.0


## set source coordinates
src_num = 21
src_coord = np.zeros((src_num, 2))
src_coord[:,0] = np.linspace( x_beg, x_end, src_num)   # linearly distributed sources
src_coord[:,1] = np.linspace(    dx,    dx, src_num)   # sources are buried at first grid depth

## source wavelet
amp0 = 1.
f0   = 5.
wavelet = np.zeros((src_num, nt))
for i in range(src_num):
    wavelet[i,:] = source_wavelet(amp0, nt, dt, f0, 'ricker')
    # one can load their own wavelet here, 
    # e.g., wavelet[i,:] = np.load('wavelet.npy')
    # the time axis should be the same as in the forward modeling

## set receiver coordinates. rec_coord is a list of receiver coordinates so that they can be different for different sources
rec_coord = []
for i in range(src_num):
    rec_xz = np.zeros((nx, 2))
    rec_xz[:,0] = np.linspace(x_beg, x_end, nx)  # linearly distributed receivers
    rec_xz[:,1] = np.linspace(   dx,    dx, nx)  # receivers are buried at first grid depth
    rec_coord.append(rec_xz)

## save acquisition files
np.save('./acquisition/vp_true.npy',       vp_true)
np.save('./acquisition/rho_true.npy',     rho_true)
np.save('./acquisition/vp_init.npy',       vp_init)
np.save('./acquisition/rho_init.npy',     rho_init)
np.save('./acquisition/grad_mask.npy',   grad_mask)
np.save('./acquisition/src_coord.npy',   src_coord)
np.save('./acquisition/wavelets.npy',      wavelet)
np.savez('./acquisition/rec_coord.npz', *rec_coord) # save receiver list using **np.savez**

print('Successfully save acquisition files!')
