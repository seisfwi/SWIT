import sys
sys.path.append('/homes/sep/haipeng/develop/SWIT-1.0/toolbox-dev/')

import os
import numpy as np
from utils import source_wavelet

## create acquisition folder
os.makedirs('acquisition', exist_ok=True)

## model size
nx = 481
nz = 141
dx = 25
nt = 2001
dt = 0.002
x_beg = 0.
x_end = (nx-1) * dx

## true model for generating obs data
vp_true = np.ones((nx, nz)) * 4000.
vp_true[220:240,60:80] = 4500
rho_true = np.power(vp_true, 0.25) * 310   # density models, (Gardner, 1974)

## initial model for FWI
vp_init = np.ones((nx, nz)) * 4000.
rho_init = np.power(vp_init, 0.25) * 310   # density models, (Gardner, 1974)

# save true and initial models
np.save('./acquisition/vp_true.npy', vp_true)
np.save('./acquisition/rho_true.npy', rho_true)
np.save('./acquisition/vp_init.npy', vp_init)
np.save('./acquisition/rho_init.npy', rho_init)

## set the gradient mask for FWI
# TODO: add gradient mask

## set source coordinates
srcn  = 24
src_coord = np.zeros((srcn, 2))
src_coord[:,0] = np.linspace( x_beg, x_end, srcn)
src_coord[:,1] = np.linspace(    dx,    dx, srcn)
# save source coordinates
np.save('./acquisition/src_coord.npy', src_coord)

## source wavelet
amp0 = 1.
f0   = 10.
wavelet = np.zeros((srcn, nt))
for i in range(srcn):
    # one can load their own wavelet here
    wavelet[i,:] = source_wavelet(amp0, nt, dt, f0, 'ricker')

# save source wavelet
np.save('./acquisition/wavelets.npy', wavelet)

## set receiver coordinates. rec_coord is a list of receiver coordinates for each source and they can be different
rec_coord = []
for i in range(srcn):
    rec_xz = np.zeros((nx, 2))
    rec_xz[:,0] = np.linspace(x_beg, x_end, nx)
    rec_xz[:,1] = np.linspace(   dx,    dx, nx)
    rec_coord.append(rec_xz)

## save receiver coordinates, note that we use np.savez here to save a list of arrays
np.savez('./acquisition/rec_coord.npz', *rec_coord)


print('Successfully save acquisition files!')