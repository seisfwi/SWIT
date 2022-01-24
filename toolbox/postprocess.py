###############################################################################
#
# SWIT: Seismic Waveform Inversion Toolbox
#
# by Haipeng Li at USTC, haipengl@mail.ustc.edu.cn
# 
# June, 2021  
#
# Postprocess module
#
###############################################################################

import numpy as np
import scipy.signal

from plot import plot_model2D
from tools import array2vector, smooth2d, vector2array


def grad_precond(simu, optim, grad, forw, back):
    ''' Gradient preconditioning
    '''
    nx = simu.model.nx
    nz = simu.model.nz
    vp = simu.model.vp
    vpmax = optim.vpmax
    marine_or_land = optim.marine_or_land
    grad_mute = optim.grad_mute
    grad_thred = optim.grad_thred
    grad_smooth = optim.grad_smooth
    grad = vector2array(grad, nx, nz)
    forw = vector2array(forw, nx, nz)
    back = vector2array(back, nx, nz)
    
    # plot forward and backward illumination
    plot_model2D(simu, forw.T, np.min(forw), 10 * np.median(forw), 'forw-illum-%03d' % optim.iter, 'my_seismic_cmap')
    plot_model2D(simu, back.T, np.min(back), 10 * np.median(back), 'back-illum-%03d' % optim.iter, 'my_seismic_cmap')

    # apply taper mask, land daming or water-layer masking
    if grad_mute > 0:
        grad *= grad_taper(nx, nz, tapersize = grad_mute, thred = grad_thred, marine_or_land=marine_or_land)

    #apply the inverse Hessian
    if min(nx, nz) > 40:      # set 40 grids in default
        span = 40
    else:                     # in case the grid number is less than 40
        span = int(min(nx, nz)/2)
    forw = smooth2d(forw, span)
    back = smooth2d(back, span)
    
    epsilon = 0.0001
    forw = forw / np.max(forw)
    back = back / np.max(back)
    precond = forw + back
    precond = precond / np.max(precond)
    precond[precond < epsilon] = epsilon
    grad = grad / precond

    # smooth the gradient
    if grad_smooth > 0:
        # exclude water-layer
        if marine_or_land in ['Marine', 'Offshore']: 
            grad[:,grad_mute:] = smooth2d(grad[:,grad_mute:], span=grad_smooth)
        # land gradient smooth
        else:
            grad = smooth2d(grad, span=grad_smooth)

    # scale the gradient properly
    grad *= vpmax / abs(grad).max()
  

    return array2vector(grad)


def grad_taper(nx, nz, tapersize=20, thred=0.05, marine_or_land='Marine'):
    ''' Gradient taper
    '''

    # for masking the water layer, use the zero threds
    if marine_or_land in ['Marine', 'Offshore']: 
        taper = np.ones((nx, nz))
        for ix in range(nx):
            taper[ix, :tapersize] = 0.0
            
    # for the land gradient damping, use the small threds
    else:
        H = scipy.signal.hamming(tapersize*2)  # gaussian window
        H = H[tapersize:]
        taper = np.zeros((nx, nz))
        for ix in range(nx):
            taper[ix, :tapersize] = H
        taper = smooth2d(taper, span=tapersize//2)
        taper /= taper.max()
        taper *= (1 - thred)
        taper = - taper + 1
        taper = taper * taper      # taper^2 is better than taper^1

    return taper
