###############################################################################
#
# SWIT v1.1: Seismic Waveform Inversion Toolbox
#
#   A Python package for seismic waveform inversion
#   By Haipeng Li at USTC & Stanford
#   Email: haipengl@mail.ustc.edu.cn, haipeng@stanford.edu 
#
#   Tools module, some of codes are from: https://github.com/rmodrak/seisflows
#
###############################################################################

import os

import numpy as np
import obspy
# import scipy.signal as _signal
import yaml
from obspy.io.segy.segy import _read_segy
from scipy.ndimage import gaussian_filter


def load_yaml(filename):
    ''' Load YAML file

    Parameters
    ----------
        filename : str
            file name
    '''
    with open(filename) as f:
        config = yaml.safe_load(f)
    
    return config


def save_yaml(filename, config):
    ''' Save YAML file

    Parameters
    ----------
        filename : str
            file name
        config : dict
            configuration dictionary
    '''
    with open(filename, 'w') as f:
        yaml.dump(config, f)


def load_float(filename):
    ''' Load bin file (float32)

    Parameters
    ----------
        filename : str
            file name

        Returns
        -------
        data : 1D array (float32)
            data to be read
    '''
    with open(filename, 'rb') as f:
        data = np.fromfile(f, dtype=np.float32)

    return data


def save_float(filename, data):
    ''' Save bin file (float32)

    Parameters
    ----------
        filename : str
            file name
        data : 1D array (float32)
            data to be saved
    '''
    with open(filename, 'wb') as f:
        data.astype(np.float32).tofile(f)


def load_su(filename):
    ''' Load Seismic Unix files. if nt > 32768, this function will not work. 
        
    Parameters
        ----------
        filename : str
            file name
    '''
    traces = obspy.read(filename, format='SU', byteorder='<')

    return traces


def save_su(filename, trace):
    ''' Save su file, support max_npts = 32767

    Parameters
    ----------
        filename : str
            file name
        trace : obspy trace
            trace to be saved
    '''
    
    # work around obspy data type conversion
    for tr in trace:
        tr.data = tr.data.astype(np.float32)

    max_delta = 0.065535
    dummy_delta = max_delta

    if trace[0].stats.delta > max_delta:
        for tr in trace:
            tr.stats.delta = dummy_delta
    
    # write data to file
    trace.write(filename, format='SU')


def load_segy(filename, convert_to_array = False):
    ''' Load segy file

    Parameters
    ----------
        filename : str
            file name
        convert_to_array : bool
            if True, convert segy traces to array
    '''

    traces = _read_segy(filename)

    return segy2array(traces) if convert_to_array else traces


def segy2array(segy_data):
    ''' Convert segy traces to array

    Parameters
    ----------
        segy_data : segyio.segyio.SegyFile
            segy data

    Returns
    -------
        data : 2D array (float32)
            data array
        dt : float
            time sample interval
    '''
    # retrieve data parameters
    ntr = len(segy_data.traces)
    nt = segy_data.traces[0].npts
    dt = segy_data.binary_file_header.sample_interval_in_microseconds / 1000000.  # in seconds

    # convert to array
    data = np.zeros((ntr, nt), dtype=np.float32)
    for itr in range(ntr):
        data[itr, :] = segy_data.traces[itr].data

    # return data
    return data, dt


def su2array(su_data):
    ''' Convert su traces to array

    Parameters
    ----------
        su_data : list of obspy trace
            su data (list of obspy trace)
    
    Returns
    -------
        data : 2D array (float32)
            data array
    '''
    # retrieve data parameters
    ntr = len(su_data)
    nt = len(su_data[0])
    dt = 1. / su_data[0].stats.sampling_rate  # in seconds

    # convert to array
    data = np.zeros((ntr, nt), dtype=np.float32)
    for itr in range(ntr):
        data[itr, :] = su_data[itr].data

    # return data
    return data, dt


def smooth2d(data, span = 10):
    ''' Smooths values on 1D/2D/3D rectangular grid

    Parameters
    ----------
        data : 2D array
            data to be smoothed
        span : int
            smoothing span

    Returns
    -------
        data : 2D array
            smoothed data
    '''
    
    return gaussian_filter(data, sigma=span//2)

    # import warnings
    # warnings.filterwarnings('ignore')

    # Z = np.copy(Z)

    # x = np.linspace(-2.*span, 2.*span, 2*span + 1)
    # y = np.linspace(-2.*span, 2.*span, 2*span + 1)
    # (X, Y) = np.meshgrid(x, y)
    # mu = np.array([0., 0.])
    # sigma = np.diag([span, span])**2.

    # # evaluate Gaussian over points of X,Y
    # D = sigma[0, 0]*sigma[1, 1] - sigma[0, 1]*sigma[1, 0]
    # B = np.linalg.inv(sigma)
    # X = X - mu[0]
    # Y = Y - mu[1]
    # F = B[0, 0]*X**2. + B[0, 1]*X*Y + B[1, 0]*X*Y + B[1, 1]*Y**2.
    # F = np.exp(-0.5*F)
    # F *= (2.*np.pi*np.sqrt(D))**(-1.)
    
    # F = F/np.sum(F)
    # W = np.ones(Z.shape)
    # Z = _signal.convolve2d(Z, F, 'same')
    # W = _signal.convolve2d(W, F, 'same')
    # Z = Z/W

    # return Z


def smooth1d(x, window_len = 11, window = 'hanning'):
    ''' Smooth the data using a window with requested size.
    
    Parameters
    ----------
        x : 1D array
            the input signal
        window_len : int
            the dimension of the smoothing window; should be an odd integer
        window : str
            the type of window from 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'

    Returns
    -------
        y : 1D array
            the smoothed signal
    '''

    if x.ndim != 1:
        raise ValueError("smooth only accepts 1 dimension arrays.")

    if x.size < window_len:
        raise ValueError("Input vector needs to be bigger than window size.")

    if window_len < 3:
        return x

    if not window in ['flat', 'hanning', 'hamming', 'bartlett', 'blackman']:
        raise ValueError(
            "Window is on of 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'")

    s = np.r_[x[window_len-1:0:-1], x, x[-2:-window_len-1:-1]]

    # moving average
    if window == 'flat':
        w = np.ones(window_len, 'd')
    else:
        w = eval('np.'+window+'(window_len)')

    y = np.convolve(w/w.sum(), s, mode='valid')

    return y


def load_model(file, nx, nz):
    ''' Load model from file (npy, dat, txt, bin)

    Parameters:
    -----------
        file: str
            Name of the model file
        nx: int
            Number of grid points in x direction
        nz: int
            Number of grid points in z direction

    Returns:
    --------
        model: 2D array
            Model array
    '''

    # output message
    # print('Loading model from file: {}'.format(file))

    if not os.path.isfile(file):
        raise FileNotFoundError('Model file not found: {}'.format(file))

    # npy file
    if file.endswith('.npy'):
        model = np.load(file)
    # plain text file
    elif file.endswith('.dat') or file.endswith('.txt'):
        model = np.loadtxt(file, dtype=np.float32)
    # binary file
    elif file.endswith('.bin'):
        model = load_float(file, nx, nz)
    else:
        msg = 'Support model file types: .npy, .dat, .txt, .bin. \n'
        err = 'Unknown model file type: {}'.format(file)
        raise ValueError(msg + '\n' + err)

    # reshape model
    try:
        model = model.reshape(nx, nz)
    except:
        msg = 'Model file does not match the specified grid size. \n'
        err = 'Model file: {} \n Grid size: ({}, {}) may be wrong!'.format(file, nx, nz)
        raise ValueError(msg + '\n' + err)

    return model


def load_waveform_data(path, nt):
    ''' Load waveform data from segy, su, or binary format and convert to array

    Parameters
    ----------
        path: str
            Path to the waveform data
        nt: int
            Number of time samples
    
    Returns
    -------
        data: 2D array
            Waveform data array
        dt: float
            Time sampling interval in seconds (None for binary data)
    '''

    # determine file type on the basis of file extension
    if os.path.exists(path + '.segy'):
        data = load_segy(path + '.segy')
        data, dt = segy2array(data)

    elif os.path.exists(path + '.su'):
        data = load_su(path + '.su')
        data, dt = su2array(data)

    elif os.path.exists(path + '.bin'):
        data = load_float(path + '.bin').reshape(-1, nt)
        dt = None

    else:
        msg = 'waveform data: {}.* not found.'.format(path)
        raise FileNotFoundError(msg)

    return data, dt
