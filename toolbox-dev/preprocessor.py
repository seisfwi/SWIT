###############################################################################
#
# SWIT v1.1: Seismic Waveform Inversion Toolbox
#
#   A Python package for seismic waveform inversion
#   Developed by Haipeng Li at USTC, updated on 2022-12-21 at Stanford
#   haipengl@mail.ustc.edu.cn, haipeng@stanford.edu
#
#   Data preprocessing module
#
###############################################################################


import time
from multiprocessing import Pool

import numpy as np
import os
from scipy.signal import butter, filtfilt

from tools import load_waveform_data, save_float


class Preprocessor(object):
    '''  preprocessor class describes the data preprocessing

    Parameters
    ----------
    filter : str
        Type of the data filter, 'lowpass', 'bandpass', 'highpass' or 'none'
    filter_low : float
        Low frequency of the data filter in Hz
    filter_high : float
        High frequency of the data filter in Hz
    mute_late_arrival : bool
        Whether to mute the late arrivals after the first break
    mute_late_size : float
        Time window of the late arrivals to be muted in seconds
    normalize_data : bool
        Whether to normalize the trace by its maximum amplitude (trace by trace)
    mute_near_offset : bool
        Whether to mute the near offset traces
    mute_near_distance : float  
        Distance of the near offset traces to be muted in meters
    mute_far_offset : bool 
        Whether to mute the far offset traces
    mute_far_distance : float
        Distance of the far offset traces to be muted in meters
    '''

    def __init__(self, filter_data = 'bandpass', filter_low = 5.0, filter_high = 10.0, 
                 mute_late_arrival = False, mute_late_size = 0.5, 
                 normalize_data = False,
                 mute_near_offset = False, mute_near_distance = 500, 
                 mute_far_offset = False, mute_far_distance = 8000):
        ''' Initialize preprocessor class
        '''

        # data filter
        self.filter_data = filter_data
        self.filter_low = filter_low
        self.filter_high = filter_high

        # mute later arrivals after the first break
        self.mute_late_arrival = mute_late_arrival
        self.mute_late_size = mute_late_size

        # data normalization
        self.normalize_data = normalize_data

        # data offset mute
        self.mute_near_offset = mute_near_offset
        self.mute_far_offset = mute_far_offset
        self.mute_near_distance = mute_near_distance
        self.mute_far_distance = mute_far_distance

        # set one thread for scipy to perform the filtering
        os.environ["OMP_NUM_THREADS"] = "1"    # $ export OMP_NUM_THREADS=1
    
        # check the parameters
        self.__check__()

        # print the information of the preprocessor
        self.__info__()


    def __info__(self):
        ''' Print the information of the preprocessor
        '''
        print('\nPreprocessor information:')
        print('    Filter data       : {}'.format(self.filter_data))
        print('    Lowcut frequency  : {} Hz'.format(self.filter_low))
        print('    Highcut frequency : {} Hz'.format(self.filter_high))
        print('    Mute late arrivals: {}'.format(self.mute_late_arrival))
        print('    Mute time window  : {} s'.format(self.mute_late_size))
        print('    Mute near offset  : {}'.format(self.mute_near_offset))
        print('    Mute near dist    : {} m'.format(self.mute_near_distance))
        print('    Mute far offset   : {}'.format(self.mute_far_offset))
        print('    Mute far dist     : {} m'.format(self.mute_far_distance))
        print('    Data normalization: {}'.format(self.normalize_data))


    def __check__(self):

        # check the data filter
        if self.filter_data.lower() not in ['lowpass', 'bandpass', 'highpass', 'none']:
            msg = 'The filter type is Lowpass, Bandpass, Highpass or None'
            err = 'Not supported filter type: {}'.format(self.filter)
            raise ValueError(msg + '\n' + err)

        # check the low and high frequency
        if self.filter_low > self.filter_high:
            raise ValueError('The low frequency is larger than the high frequency')

        # check the late arrival mute
        if self.mute_late_arrival and self.mute_late_size <= 0:
            raise ValueError('The late arrival mute window should be larger than 0')


        # check the data offset mute
        if self.mute_near_offset > self.mute_far_distance:
            raise ValueError('The short offset distance is larger than the long offset distance')


    def run(self, data_path = None, src_num = None, mpi_num = 1, nt = None, dt = None, src_coord = None, rec_coord = None):
        ''' run the preprocessor to process the data in the provided directory

        Parameters
        ----------
        data_path: str
            the path of the data directory
        src_num: int
            the number of sources
        mpi_num: int
            the number of MPI processes
        nt: int
            the number of time samples
        dt: float
            the time sampling interval in seconds
        src_coord: 2D array of float
            the source coordinates in meters
        rec_coord: lisft of 2D array of float
            the receiver coordinates in meters
        '''

        # print('Start preprocessing the data in: {} ...'.format(data_path))
        # create a pool of processes
        pool = Pool(mpi_num)

        # loop over the sources
        for isrc in range(src_num):
            # get the directory of the data
            load_path = os.path.join(data_path, 'src{}/sg'.format(isrc+1))
            
            # get the offset
            offset = src_coord[isrc, 0] - rec_coord[isrc][:, 0]
            
            # process the data
            pool.apply_async(self.process_workflow_it, (load_path, nt, dt, offset,) )
            
            # wait for a while
            time.sleep(0.01)

        # close the pool and wait for all processes to be done
        pool.close()

        # block at this line until all processes are done
        pool.join()


    def process_workflow_it(self, load_path, nt, dt, offset):
        ''' process the data in the provided directory
        '''

        # load data
        trace, _ = load_waveform_data(load_path, nt)

        # get parameters
        ntrace, nt = trace.shape

        if len(offset) != ntrace:
            raise RuntimeError('offset and trace are not consistant\n')

        # dafault parameters
        length = 200
    
        # mute short traces
        if self.mute_near_offset:
            trace[abs(offset) < self.mute_near_distance,:] = 0.
        
        # mute far traces
        if  self.mute_far_offset:
            trace[abs(offset) > self.mute_far_distance,:] = 0.
        
        # pick the first arrival and mute late arrival
        if self.mute_late_arrival:

            # calculate the first break in grid size
            pick = brutal_picker(trace) + np.ceil(self.mute_late_size/dt)
            
            # safeguard window
            itmin = (pick  - length/2).astype(int)
            itmax = (itmin + length  ).astype(int)

            # construct the taper
            taper = np.ones(nt)
            win = np.sin(np.linspace(0, np.pi, 2*length)) [0:length]
            
            # mute late arrivals
            for i in range(ntrace):
                trace[i,:] *= (1.0 - custom_mask(itmin[i], itmax[i], nt, length, taper, win))

        # scipy filtering
        if self.filter_data.lower() in ['bandpass']:
            for i in range(ntrace):
                trace[i] = bandpass_filter(trace[i], self.filter_low, self.filter_high, dt, order=4)
        
        elif self.filter_data.lower() in ['lowpass']:
            for i in range(ntrace):
                trace[i] = lowpass_filter(trace[i], self.filter_low, dt, order=4)
        
        elif self.filter_data.lower() in ['highpass']:
            for i in range(ntrace):
                trace[i] = highpass_filter(trace[i], self.filter_high, dt, order=4)
        
        elif self.filter_data.lower() in ['none']:
            pass
        
        else:
            raise RuntimeError('Not supported filter type: {}\n'.format(self.filter_data))
        
        # normalization
        if self.normalize_data:
            for i in range(ntrace):
                w = 0.
                w = np.max(abs(trace[i]))
                if w > 0:    
                    trace[i] /= w
        
        # always save the processed data as binary file, note that the data is flattened
        save_float(load_path + '_processed.bin', trace.flatten())


def brutal_picker(trace, threshold=0.001):
    ''' pick the first arrival based on the amplitude of the trace 

    Parameters
    ----------
    trace: 2D array of float
        the seismic traces
    
    Returns
    -------
    pick: 1D array of int
        the index of the first arrival
    '''

    pick = [(abs(trace[i]) > threshold * np.max(abs(trace[i]))).argmax(axis=-1) 
        for i in range(len(trace))]

    return np.array(pick)


def custom_mask(itmin, itmax, nt, length, mask, win):
    ''' Constructs tapered mask that can be applied to trace to
        mute early or late arrivals.
    
    Parameters
    ----------
    itmin: int
        the start index of the mask
    itmax: int
        the end index of the mask
    nt: int 
        the number of time samples
    length: int
        the length of the taper
    mask: 1D array of float
        the mask to be applied
    win: 1D array of float
        the taper

    Returns
    -------
    mask: 1D array of float
        the mask to be applied
    '''

    if 1 < itmin < itmax < nt:
        mask[0:itmin] = 0.
        mask[itmin:itmax] = win*mask[itmin:itmax]
    elif itmin < 1 <= itmax:
        mask[0:itmax] = win[length-itmax:length]*mask[0:itmax]
    elif itmin < nt < itmax:
        mask[0:itmin] = 0.
        mask[itmin:nt] = win[0:nt-itmin]*mask[itmin:nt]
    elif itmin > nt:
        mask[:] = 0.

    return mask


def bandpass_filter(data, lowcut, highcut, dt, order=4):
    '''Butterworth bandpass filter from scipy Cookbook
    
    Parameters
    ----------
    data: 1D array of float
        the data to be filtered
    lowcut: float
        the lowcut frequency
    highcut: float
        the highcut frequency
    dt: float
        the sampling interval
    order: int
        the order of the filter (default: 4)
    '''

    fs = 1.0 / dt
    nyq = 0.5 * fs
    low = lowcut / nyq
    high = highcut / nyq
    b, a = butter(order, [low, high], btype='band')
    data_bp = filtfilt(b, a, data)

    return np.asarray(data_bp, dtype=np.float32)


def lowpass_filter(data, highcut, dt, order=4):
    '''Butterworth lowpass filter from scipy Cookbook

    Parameters
    ----------
    data: 1D array of float
        the data to be filtered
    highcut: float
        the highcut frequency
    dt: float
        the sampling interval
    order: int
        the order of the filter (default: 4)
    '''

    fs = 1.0 / dt
    nyq = 0.5 * fs
    high = highcut / nyq
    b, a = butter(order, high, btype='low')
    data_lp = filtfilt(b, a, data)

    return np.asarray(data_lp, dtype=np.float32)



def highpass_filter(data, lowcut, dt, order=4):
    '''Butterworth highpass filter from scipy Cookbook
    
    Parameters
    ----------
    data: 1D array of float
        the data to be filtered
    lowcut: float
        the lowcut frequency
    dt: float
        the sampling interval
    order: int
        the order of the filter (default: 4)
    '''

    fs = 1.0 / dt
    nyq = 0.5 * fs
    low = lowcut / nyq
    b, a = butter(order, low, btype='hp')
    data_hp = filtfilt(b, a, data)

    return np.asarray(data_hp, dtype=np.float32)


