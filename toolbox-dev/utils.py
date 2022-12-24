###############################################################################
#
# SWIT v1.1: Seismic Waveform Inversion Toolbox
#
#   A Python package for seismic waveform inversion
#   Developed by Haipeng Li at USTC, updated on 2022-12-21 at Stanford
#   haipengl@mail.ustc.edu.cn, haipeng@stanford.edu
#
#   Solver module
#
###############################################################################


import numpy as np



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
        raise ValueError(msg + err)

    return wavelet
