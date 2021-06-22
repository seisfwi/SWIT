###############################################################################
#
# SWIT: Seismic Waveform Inversion Toolbox
#
# by Haipeng Li at USTC, haipengl@mail.ustc.edu.cn
# 
# June, 2021  
#
# Preprocess module, some of codes are from: https://github.com/rmodrak/seisflows 
#
###############################################################################

import time
from multiprocessing import Pool

import numpy as np
import obspy
from scipy.signal import butter, filtfilt

from tools import add_su_header, get_su_parameter, loadsu, savesu, su2array


def process_workflow(simu, optim, simu_type = 'syn', use_first_break_in_su_header = False):
    ''' process all shot gathers
    '''
    
    ## get prapred 
    nt = simu.model.nt
    dt = simu.model.dt 
    nproc = simu.system.mpiproc
    srcn = simu.source.n
    datapath = simu.system.homepath + 'data/' + simu_type + '/'

    # scipy would open too mnewspaper threads, reduce the Pool number
    pool = Pool(nproc)
    for isrc in range(srcn):
        loadpath = datapath + 'src%d_sg.su'      % ((isrc+1))
        savepath = datapath + 'src%d_sg_proc.su' % ((isrc+1))
        pool.apply_async(process_workflow_serial, (optim, loadpath, savepath, isrc, nt, dt, use_first_break_in_su_header,) )
        time.sleep(0.01)
    pool.close()  # close pool
    pool.join()   # block at this line until all processes are done


def process_workflow_serial(optim, loadpath, savepath, isrc, nt, dt, use_first_break_in_su_header):
    ''' process shot gather in serial
    '''

    ## load the data
    trace = loadsu(loadpath)
    trace = add_su_header(trace, nt, dt, isrc, 'pressure')

    ## process workflow
    apply_normalize(optim, trace)
    apply_mute(optim, trace, use_first_break_in_su_header)
    apply_filter(optim, trace)
    apply_normalize(optim, trace)

    ## save the data
    savesu(savepath, trace)


def apply_filter(optim, trace):
    
    # get parameters
    recn, nt, dt = get_su_parameter(trace)
    fre_filter = optim.fre_filter
    fre_low = optim.fre_low
    fre_high = optim.fre_high

    # scipy filtering
    if fre_filter.lower() in ['none']:
        pass
    else:
        # convert into numpy array
        data = su2array(trace)

        if fre_filter.lower() in ['bandpass']:
            data_proc = butter_bandpass_filter(data, fre_low, fre_high, dt, order=4)
        elif fre_filter.lower() in ['lowpass']:
            data_proc = butter_lowpass_filter(data, fre_low, dt, order=4)
        elif fre_filter.lower() in ['highpass']:
            data_proc = butter_highpass_filter(data, fre_high, dt, order=4)

        # replace the original waveform in the SU stream
        for irec in range(recn):
            trace[irec].data[:] = data_proc[irec]

        # detrend and taper
        for tr in trace:
            tr.detrend('demean')
            tr.detrend('linear')
            tr.taper(0.05, type='hann')


def apply_mute(optim, trace, use_first_break_in_su_header):
    ''' apply time window and offset window mute.
    '''
    # set parameters
    length = 200
    recn, nt, dt = get_su_parameter(trace)

    mute_late_arrival = optim.mute_late_arrival
    mute_late_window  = optim.mute_late_window
    mute_offset_short = optim.mute_offset_short
    mute_offset_long  = optim.mute_offset_long
    mute_offset_short_dis = optim.mute_offset_short_dis         # (units: m)
    mute_offset_long_dis  = optim.mute_offset_long_dis          # (units: m)

    for tr in trace:        
        # mute short offset traces
        if mute_offset_short:
            mute_offset(tr, mute_offset_short_dis, 'short')

        # mute long offset traces
        if mute_offset_long:
            mute_offset(tr, mute_offset_long_dis, 'long')

    # pick the first arrival and mute the late arrival
    if mute_late_arrival:
        # use the firsy arrival in the header
        if use_first_break_in_su_header:
            itrace = 0
            pick = np.zeros(len(trace))
            for tr in trace:
                pick[itrace] = int((tr.stats.su.trace_header.mute_time_start_time_in_ms/1000 + mute_late_window)/dt)
                itrace +=1
        # pick the first arrival in a brutal manner 
        else:
            pick = brutal_picker(trace) + np.ceil(mute_late_window/dt)

        itrace = 0
        for tr in trace:
            # set mute window
            itmin = int(pick[itrace] - length/2)
            itmax = int(itmin + length)
            # apply mute
            mute_arrival(tr, itmin, itmax, 'late', nt, length)
            itrace +=1


def apply_normalize(optim, trace):
    ''' apply normalize using L1 or L2 norms on whole event or single trace
    '''

    # parameter
    normalize = optim.normalize

    if 'L1-Trace' in normalize:
        # normalize each trace by its L1 norm
        for tr in trace:
            w = 0.
            w = np.linalg.norm(tr.data, ord=1)
            if w > 0:
                tr.data /= w
    elif 'L2-Trace' in normalize:
        # normalize each trace by its L2 norm
        for tr in trace:
            w = 0.
            w = np.linalg.norm(tr.data, ord=2)
            if w > 0:
                tr.data /= w
    elif 'Max-Trace' in normalize:
        # normalize each trace by the max value
        for tr in trace:
            w = 0.
            w = np.max(abs(tr.data))
            if w > 0:
                tr.data /= w
    else:
        pass
    
    if 'L1-Event' in normalize:
        # normalize event by L1 norm of all data
        w = 0.
        for tr in trace:
            w += np.linalg.norm(tr.data, ord=1)
        for tr in trace:
            tr.data /= w
    elif 'L2-Event' in normalize:
        # normalize event by L2 norm of all data
        w = 0.
        for tr in trace:
            w += np.linalg.norm(tr.data, ord=2)
        for tr in trace:
            tr.data /= w
    else:
        pass


def brutal_picker(trace):
    ''' pick the first arrival 
    '''
    # set parameter
    # threds = 0.000001
    # trace = su2array(trace)
    # threds *= np.median(abs(trace))

    threds = 0.01
    trace = su2array(trace)

    return (abs(trace) > threds).argmax(axis=-1)


def mute_offset(trace, mute_dist, mutetype):
    ''' Mutes traces having

         mutetype = short:   || s - r || < mute_dist  or 
         mutetype = long :   || s - r || > mute_dist
         
        where || s - r || is the offset between source and receiver and
        DIST is a user-supplied cutoff
    '''

    # offset
    offset = abs(trace.stats.distance)

    if mutetype in ['short']:
        if offset < mute_dist:
            trace.data[:] = 0.
    elif mutetype in ['long']:
        if offset > mute_dist:
            trace.data[:] = 0.
    else:
        raise ValueError('Wrong mutetype')


def mute_arrival(trace, itmin, itmax, mutetype, nt, length):
    ''' Applies tapered mask to record section, muting early or late arrivals
    '''

    # apply tapered mask
    if mutetype in ['early']:
        trace.data *= mask(itmin, itmax, nt, length)
    elif mutetype in ['late']:
        trace.data *= (1. - mask(itmin, itmax, nt, length))
    else:
        raise ValueError('Wrong mutetype')


# functions acting on individual traces
def mask(itmin, itmax, nt, length):
    ''' Constructs tapered mask that can be applied to trace to
        mute early or late arrivals.
    '''
    mask = np.ones(nt)
    # construct taper
    win = np.sin(np.linspace(0, np.pi, 2*length))
    win = win[0:length]

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

# scipy bandpass filter
def butter_bandpass_filter(data, lowcut, highcut, dt, order=4):
    '''Butterworth bandpass filter from scipy Cookbook
    '''
    fs = 1.0 / dt
    nyq = 0.5 * fs
    low = lowcut / nyq
    high = highcut / nyq
    b, a = butter(order, [low, high], btype='band')
    data_bp = filtfilt(b, a, data)

    return data_bp

# scipy lowpass filter
def butter_lowpass_filter(data, highcut, dt, order=4):
    '''Butterworth lowpass filter from scipy Cookbook
    '''

    fs = 1.0 / dt
    nyq = 0.5 * fs
    high = highcut / nyq
    b, a = butter(order, high, btype='low')
    data_lp = filtfilt(b, a, data)

    return data_lp

# scipy highpass filter
def butter_highpass_filter(data, lowcut, dt, order=4):
    '''Butterworth highpass filter from scipy Cookbook
    '''

    fs = 1.0 / dt
    nyq = 0.5 * fs
    low = lowcut / nyq
    b, a = butter(order, low, btype='hp')
    data_hp = filtfilt(b, a, data)

    return data_hp

