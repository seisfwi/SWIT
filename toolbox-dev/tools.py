###############################################################################
#
# SWIT: Seismic Waveform Inversion Toolbox
#
# by Haipeng Li at USTC, haipengl@mail.ustc.edu.cn
#
# June, 2021
#
# Tools module, some of codes are from: https://github.com/rmodrak/seisflows
#
###############################################################################

import os
import yaml

import numpy as np
import obspy
import scipy.signal as _signal
from obspy.io.segy.segy import _read_segy



class Config_object(object):
    ''' A class to convert dictionary to object
    '''
    def __init__(self, d):
        for k, v in d.items():
            for k1, v1 in v.items():
                setattr(self, k1, v1)

            # if isinstance(k, (list, tuple)):
            #     setattr(self, k, [Config_object(x) if isinstance(x, dict) else x for x in v])
            # else:
            #     setattr(self, k, Config_object(v) if isinstance(v, dict) else v)


def load_yaml(filename):
    ''' A function to read YAML file
    '''

    try:
        with open(filename) as f:
            config = yaml.safe_load(f)
        
    except FileNotFoundError:
        raise FileNotFoundError('File not found: {}'.format(filename))

    return config


def save_yaml(filename, config):
    ''' A function to save YAML file
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
    try:
        fp = open(filename, 'rb')
        data = np.fromfile(fp, dtype=np.float32)
        fp.close()
    except FileNotFoundError:
        raise FileNotFoundError('File not found: {}'.format(filename))

    return np.array(data)


def save_float(filename, data):
    ''' Save bin file (float32)

    Parameters
    ----------
    filename : str
        file name
    data : 1D array (float32)
        data to be saved
    '''

    fp = open(filename, "wb")
    np.asarray(data, dtype=np.float32).tofile(fp)
    fp.close()


def load_su(filename):
    ''' Reads Seismic Unix files. if nt > 32768, this function will not work. 
        Find solution in Seisflows/plugins/writer.py or reader.py
    '''
    try: 
        traces = obspy.read(filename, format='SU', byteorder='<')
    except FileNotFoundError:
        raise FileNotFoundError('File not found: {}'.format(filename))

    return traces


def save_su(filename, trace):
    ''' save su file,     max_npts = 32767
    '''

    for tr in trace:
        # work around obspy data type conversion
        tr.data = tr.data.astype(np.float32)

    max_delta = 0.065535
    dummy_delta = max_delta

    if trace[0].stats.delta > max_delta:
        for tr in trace:
            tr.stats.delta = dummy_delta
    # write data to file
    trace.write(filename, format='SU')


def load_segy(filename, convert_to_array=False):
    ''' Reads segy files. 
    '''
    try:
        traces = _read_segy(filename)
    except FileNotFoundError:
        raise FileNotFoundError('File not found: {}'.format(filename))

    if convert_to_array:
        traces = segy2array(traces)

    return traces


def segy2array(segy_data):
    ''' convert segy traces to array
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
    ''' convert su traces to array
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


def gauss2(X, Y, mu, sigma, normalize=True):
    ''' Evaluates Gaussian over points of X,Y
    '''
    D = sigma[0, 0]*sigma[1, 1] - sigma[0, 1]*sigma[1, 0]
    B = np.linalg.inv(sigma)
    X = X - mu[0]
    Y = Y - mu[1]
    Z = B[0, 0]*X**2. + B[0, 1]*X*Y + B[1, 0]*X*Y + B[1, 1]*Y**2.
    Z = np.exp(-0.5*Z)

    if normalize:
        Z *= (2.*np.pi*np.sqrt(D))**(-1.)

    return Z


def smooth2d(Z, span = 10):
    ''' Smooths values on 2D rectangular grid
    '''
    import warnings
    warnings.filterwarnings('ignore')

    Z = np.copy(Z)

    x = np.linspace(-2.*span, 2.*span, 2*span + 1)
    y = np.linspace(-2.*span, 2.*span, 2*span + 1)
    (X, Y) = np.meshgrid(x, y)
    mu = np.array([0., 0.])
    sigma = np.diag([span, span])**2.
    F = gauss2(X, Y, mu, sigma)
    F = F/np.sum(F)
    W = np.ones(Z.shape)
    Z = _signal.convolve2d(Z, F, 'same')
    W = _signal.convolve2d(W, F, 'same')
    Z = Z/W

    return Z


def smooth1d(x, window_len=11, window='hanning'):
    ''' smooth the data using a window with requested size.
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

    # print(len(s))
    if window == 'flat':  # moving average
        w = np.ones(window_len, 'd')
    else:
        w = eval('np.'+window+'(window_len)')

    y = np.convolve(w/w.sum(), s, mode='valid')

    return y


def load_model(file, nx, nz):
    ''' A function to load model from file

        Parameters:
        -----------
        file: str
            Name of the model file
        nx: int
            Number of grid points in x direction
        nz: int
            Number of grid points in z direction
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
    ''' load waveform data from segy, su, or binary files, and convert the data 
        to numpy array. The file type is detected on the fly based on the file 
        extension name (.segy, .su, .bin).
    '''

    # load the waveform data in segy, su, or binary format

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
