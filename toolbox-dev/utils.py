###############################################################################
#
# SWIT v1.1: Seismic Waveform Inversion Toolbox
#
#   A Python package for seismic waveform inversion
#   By Haipeng Li at USTC & Stanford
#   Email: haipengl@mail.ustc.edu.cn, haipeng@stanford.edu 
#
#   Utils module specfically for waveform modeling, inversion, and migration
#
###############################################################################

import numpy as np
import scipy.signal
from tools import smooth2d


def source_wavelet(amp0, nt, dt, f0, source_type):
    ''' Source wavelet

    Parameters:
    -----------
        amp0: float
            Amplitude of the source wavelet
        nt: int
            Number of time samples
        dt: float
            Time sampling interval
        f0: float
            Dominant frequency of the source wavelet
        source_type: str
            Type of the source wavelet

    Returns:
    --------
        wavelet: 1D array
            Source wavelet
    '''

    # set up time axis and wavelet
    t = np.linspace(0, dt*nt, num=nt, endpoint=False)
    wavelet = np.zeros_like(t)

    # Ricker wavelet
    if source_type.lower() in ['ricker']:
        t0 = 1.2 / f0
        temp = (np.pi*f0) ** 2
        wavelet = amp0 * (1 - 2 * temp * (t - t0) ** 2) * np.exp(- temp * (t - t0) ** 2)
    
    # Gaussian wavelet
    elif source_type.lower() in ['gaussian']:
        # t0 = 1.2 / f0
        # wavelet = amp0 * np.exp(- (t - t0) ** 2 / (2 * t0 ** 2))
        raise NotImplementedError('Gaussian source wavelet not implemented yet.')

    # Ramp wavelet
    elif source_type.lower() in ['ramp']:
        # t0 = 1.2 / f0
        # wavelet = amp0 * 0.5 * (1. + np.tanh(t / t0))
        raise NotImplementedError('Ramp source wavelet not implemented yet.')
        
    # Unknown source type
    else:
        msg = 'Support source types: Rikcer, Guassian, Ramp. \n'
        err = 'Unknown source type: {}'.format(source_type)
        raise ValueError(msg + '\n' + err)

    return wavelet


def preconditioner(for_illum, adj_illum, epsilon = 0.0001):
    ''' Generate preconditioner

    Parameters:
    -----------
        for_illum: 2D array
            Forward illumination
        adj_illum: 2D array
            Adjoint illumination
        epsilon: float
            Small value to avoid zero division

    Returns:
    --------
        precond: 2D array
            Preconditioner
    '''

    # set proper smoothing size, 30 grids in default
    if min(for_illum.shape) > 30:
        smooth_size = 30
    # if the grid number is less than 30, use half of the grid number
    else:
        smooth_size = int(min(for_illum.shape)/2)

    # smooth the forward and adjoint illuminations
    for_illum = smooth2d(for_illum, smooth_size)
    adj_illum = smooth2d(adj_illum, smooth_size)

    # apply the approximated inverse Hessian
    precond = for_illum / np.max(for_illum) + adj_illum / np.max(adj_illum)
    precond /=  np.max(precond)
    precond[precond < epsilon] = epsilon

    return precond


def generate_mask(nx, nz, acquisition_type, threshold = 0.05, mask_size = 20):
    ''' Generate mask

    Parameters:
    -----------
        nx: int
            Number of grids in x direction
        nz: int
            Number of grids in z direction
        acquisition_type: str
            Type of the mask, 'land' or 'marine'
        threshold: float
            Threshold of the mask for tapering in the land case, 0.05 in default
        mask_size: int
            Size of the mask

    Returns:
    --------
        mask: 2D array
            Mask
    '''

    # generate a must mask for marine case, mask_size is the number of grids of water
    if acquisition_type.lower() in ['marine']:
        mask = np.ones((nx, nz))
        for ix in range(nx):
            mask[ix, :mask_size] = 0.0

    # generate a damping mask for land case to suppress the strong gradients around the source
    elif acquisition_type.lower() in ['land']:
       
        mask = np.zeros((nx, nz))
        
        # gaussian window
        H = scipy.signal.hamming(mask_size*2)[mask_size:]
        for ix in range(nx):
            mask[ix, :mask_size] = H

        # smooth the mask
        mask = smooth2d(mask, span=mask_size//2)

        # scale the mask to 0-1
        mask /= mask.max()
        mask *= (1 - threshold)
        mask = - mask + 1
        mask = mask * mask      # taper^2 is better than taper^1
    else:
        msg = 'Support mask types: land, marine. \n'
        err = 'Unknown mask type: {}'.format(acquisition_type)
        raise ValueError(msg + '\n' + err)
    
    return mask
