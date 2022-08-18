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


import copy
import json
import os

import numpy as np
import obspy
import scipy.signal as _signal
from obspy.core import UTCDateTime
from scipy.signal import hilbert as analytic


class NumpyEncoder(json.JSONEncoder):
    ''' Numpy Encoder.
    '''
    def default(self, obj):
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        return json.JSONEncoder.default(self, obj)


class Struct:
    ''' Json Decoder. 
    '''

    def __init__(self, **entries):
        self.__dict__.update(entries)


def saveparjson(simu, optim):
    ''' Save all the simulation and optimaztion parameters as json files.
    '''

    path = simu.system.homepath+'parfile/'
    savejson(path + 'model.json',  simu.model)
    savejson(path + 'source.json', simu.source)
    savejson(path + 'reciver.json', simu.receiver)
    savejson(path + 'system.json', simu.system)
    savejson(path + 'optimize.json', optim)


def savejson(filename, obj):
    ''' Save object using json. 
    '''
    with open(filename, 'w') as file:
        json.dump(obj.__dict__, file, sort_keys=True, indent=4, cls=NumpyEncoder)


def loadjson(filename):
    ''' Load object using json. 
    '''
    with open(filename, 'r') as file:
        # parse the json file into diction and the convert into class
        obj = Struct(**json.load(file))
        # convert to numpy array for convience
        for att in obj.__dict__:
            setattr(obj, att, np.array(getattr(obj, att)))
        return obj


def savebinfloat32(filename, obj):
    ''' Save bin file (float32)
    '''
    fp = open(filename, "wb")
    np.asarray(obj, dtype=np.float32).tofile(fp)
    fp.close()


def loadbinfloat32(filename):
    ''' Load bin file (float32)
    '''
    fp = open(filename, 'rb')
    obj = np.fromfile(fp, dtype=np.float32)
    fp.close()

    return np.array(obj)


def loadsu(filename):
    ''' Reads Seismic Unix files. if nt > 32768, this function will not work. 
        Find solution in Seisflows/plugins/writer.py or reader.py
    '''
    traces = obspy.read(filename, format='SU', byteorder='<')

    return traces


def savesu(filename, trace):
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


def add_su_header(trace, nt, dt, isrc, channel):
    ''' add su header
    '''
    # get parameters

    irec = 0
    for tr in trace:
        # add headers
        tr.stats.network = 'FWI'
        tr.stats.station = '%d'%irec
        tr.stats.location = 'Source-%d' % isrc
        tr.stats.channel = channel
        tr.stats.starttime = UTCDateTime("2021-01-01T00:00:00")
        tr.stats.sampling_rate = 1./dt
        tr.stats.distance = tr.stats.su.trace_header.group_coordinate_x - tr.stats.su.trace_header.source_coordinate_x
        t0 = tr.stats.starttime 
        tr.trim(starttime=t0, endtime=t0 + dt*(nt-1), pad=True, nearest_sample=True, fill_value=0.)

        irec += 1
    return trace



def convert_wavelet_su(dt, wavelet, srcx):
    ''' convert wavelet to SU stream
    '''
    srcn  = np.size(wavelet, 0)

    wavelet_su = array2su(srcn, dt, wavelet)

    ishot = 0
    for iwvlt in wavelet_su:
        iwvlt.stats.distance = srcx[ishot]
        ishot+=1

    return wavelet_su


def get_offset(trace):
    ''' get offset from trace header
    '''
    offset = np.zeros(len(trace))
    irec = 0
    for tr in trace:
        offset[irec] = tr.stats.su.trace_header.group_coordinate_x - tr.stats.su.trace_header.source_coordinate_x
        irec += 1

    return offset


def array2su(recn, dt, traces_array):
    ''' convert array data to su stream
    '''
    from obspy.core.util import get_example_file

    # get a example stream and trace
    filename = get_example_file("1.su_first_trace")
    stream = obspy.read(filename, format='SU', byteorder='<')
    tr_example = stream[0]

    traces = []
    if recn > 1:
        for irec in range(recn):
            tr = copy.deepcopy(tr_example)
            tr.data = traces_array[irec,:]
            tr.stats.sampling_rate = 1./dt
            tr.stats.starttime = UTCDateTime("2021-01-01T00:00:00")
            # add trace
            traces += [tr]
    else: # signle trace
            tr = copy.deepcopy(tr_example)
            tr.data = traces_array
            tr.stats.sampling_rate = 1./dt
            tr.stats.starttime = UTCDateTime("2021-01-01T00:00:00")
            # add trace
            traces += [tr]
    # obspy stream
    return obspy.Stream(traces = traces)


def su2array(su_data):
    ''' extract shot gather data into array(srcn, recn, nt)
    '''

    recn, nt, _ = get_su_parameter(su_data)
    data = np.zeros((recn, nt), dtype=np.float32)
    for irec in range(recn):
        data[irec, :] = su_data[irec].data[:]

    return data


def get_su_parameter(trace):
    ''' get general parameters from hearder.
    '''
    recn = len(trace)
    nt = len(trace[0])
    dt = 1. / trace[0].stats.sampling_rate

    return recn, nt, dt


def hilbert(w):
    ''' hilbert transformation
    '''
    return np.imag(analytic(w))


def array2vector(array):
    ''' array to vector
    '''
    nx, nz = array.shape[0:2]

    return array.reshape(nx * nz)


def vector2array(vector, nx, nz):
    ''' vector to array
    '''
    return vector.reshape((nx, nz))


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


def smooth2d(Z, span=10):
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


def smooth1d(x,window_len=11,window='hanning'):
    ''' smooth the data using a window with requested size.
    '''

    if x.ndim != 1:
        raise ValueError("smooth only accepts 1 dimension arrays.")

    if x.size < window_len:
        raise ValueError("Input vector needs to be bigger than window size.")

    if window_len<3:
        return x

    if not window in ['flat', 'hanning', 'hamming', 'bartlett', 'blackman']:
        raise ValueError("Window is on of 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'")

    s=np.r_[x[window_len-1:0:-1],x,x[-2:-window_len-1:-1]]

    #print(len(s))
    if window == 'flat': #moving average
        w=np.ones(window_len,'d')
    else:
        w=eval('np.'+window+'(window_len)')

    y=np.convolve(w/w.sum(),s,mode='valid')

    return y


def savetxt(filename, data, datatype):
    ''' save txt file
    '''
    if os.path.exists(filename):
        data = np.append(np.loadtxt(filename), data)

    data = [data]
    if datatype == 'int':
        np.savetxt(filename, np.array(data), fmt='%d', delimiter='\n')
    elif datatype == 'float':
        np.savetxt(filename, np.array(data), fmt='%.4e', delimiter='\n')
    else:
        raise ValueError('unsupport data type')


def cleandata(datapath):
    ''' clean data
    '''
    if os.listdir(path = datapath) != []:
        os.system('rm -r %s' % datapath)
        os.system('mkdir %s' % datapath)



def model_misfit(true_model, inv_model):
    ''' Estimate the model misfit
    '''
    num = np.linalg.norm(true_model-inv_model, ord=2)
    den = np.linalg.norm(true_model, ord=2)

    return num/den


def save_inv_scheme(simu, optim, inv_scheme):
    ''' save current inversion results
    '''
    nx = simu.model.nx
    nz = simu.model.nz
    it = optim.iter
    homepath = simu.system.homepath
    outputpath = homepath + 'outputs'

    savebinfloat32(outputpath+'/gradient/grad-%d.bin' % it, inv_scheme['g_now'])
    savebinfloat32(outputpath+'/direction/dirc-%d.bin' % it, inv_scheme['d_now'])
    savebinfloat32(outputpath+'/velocity/vp-%d.bin' % it, inv_scheme['v_now'])
    if optim.iter == 1:
        savetxt(outputpath+'/misfit_data.dat', inv_scheme['f_old'], 'float')
    savetxt(outputpath+'/misfit_data.dat', inv_scheme['f_now'], 'float')
    savetxt(outputpath+'/line_search_step.dat', inv_scheme['ls_step'], 'float')
    savetxt(outputpath+'/line_search_iteration.dat', inv_scheme['ls_iter'], 'int')

    if  'waveform_misfit' in inv_scheme.keys():
        savetxt(outputpath+'/waveform_misfit.dat', inv_scheme['waveform_misfit'], 'float')


def loadfile_gui(filename, nx, nz):
    ''' load file for GUI
    '''

    data = np.zeros((nx, nz))
    
    if filename.endswith('bin'):
        data = loadbinfloat32(filename).reshape((nx,nz))
    elif filename.endswith('dat'):
        data = np.loadtxt(filename)
    else:
        raise ValueError('not supported model file, use: bin, txt, dat files instead')

    return data