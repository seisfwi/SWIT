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


import numpy as np
import os

class Preprocessor(object):
    ''' 
        preprocessor class describes the data preprocessing
    '''

    def __init__(self, fre_filter, fre_low, fre_high, 
                 mute_late_arrival, mute_late_window, normalize,
                 mute_offset_short, mute_offset_long, 
                 mute_offset_short_dis, mute_offset_long_dis):
        '''
            initialize source class

            input:
                fre_filter: filter the data or not
                fre_low: low frequency of the bandpass filter
                fre_high: high frequency of the bandpass filter
                mute_late_arrival: mute late arrivals or not
                mute_late_window: window length for late arrival mute
                normalize: normalize the data or not
                mute_offset_short: mute short offset or not
                mute_offset_long: mute long offset or not
                mute_offset_short_dis: short offset distance    (units: m)
                mute_offset_long_dis: long offset distance      (units: m)
        '''

        # data filter
        self.fre_filter = fre_filter
        self.fre_low = fre_low
        self.fre_high = fre_high

        # pick first break and mute later arrivals
        self.mute_late_arrival = mute_late_arrival
        self.mute_late_window = mute_late_window

        # data normalization
        self.normalize = normalize

        # data offset mute
        self.mute_offset_short = mute_offset_short
        self.mute_offset_long = mute_offset_long
        self.mute_offset_short_dis = mute_offset_short_dis           
        self.mute_offset_long_dis = mute_offset_long_dis             

        # set one thread for scipy to perform the filtering, which is equivalent to 
        #   $ export OMP_NUM_THREADS=1
        os.environ["OMP_NUM_THREADS"] = "1" 

        # print the information of the preprocessor
        self.__info__()

    def __info__(self):
        '''
            print the information of the preprocessor
        '''

        print('Preprocessor information:')
        print('    Data filter: {}'.format(self.fre_filter))
        print('    Low frequency: {} Hz'.format(self.fre_low))
        print('    High frequency: {} Hz'.format(self.fre_high))
        print('    Mute late arrivals: {}'.format(self.mute_late_arrival))
        print('    Late arrival mute window: {} s'.format(self.mute_late_window))
        print('    Data normalization: {}'.format(self.normalize))
        print('    Mute short offset: {}'.format(self.mute_offset_short))
        print('    Mute long offset: {}'.format(self.mute_offset_long))
        print('    Short offset distance: {} m'.format(self.mute_offset_short_dis))
        print('    Long offset distance: {} m'.format(self.mute_offset_long_dis))



    def __check__(self):

        # check the data filter
        if self.fre_filter not in ['Lowpass', 'Bandpass', 'Highpass', 'None']:
            raise ValueError(' The filter type is Lowpass, Bandpass, Highpass or None\n  \
                Not supported frequency filter: {}'.format(self.fre_filter))

        if self.fre_low > self.fre_high:
            raise ValueError('The low frequency is larger than the high frequency')

        # check the late arrival mute
        if self.mute_late_arrival and self.mute_late_window <= 0:
            raise ValueError('The late arrival mute window should be larger than 0')

        # check the data normalization
        if self.normalize not in ['Max-Trace', 'L1-Event', 'L2-Event', 'L1-Trace', 'L2-Trace', 'None']:
            raise ValueError('The normalization type is Max-Trace, L1-Event, L2-Event, L1-Trace, L2-Trace or None \n \
                Not supported normalization type: {}'.format(self.normalize))

        # check the data offset mute
        if self.mute_offset_short_dis > self.mute_offset_long_dis:
            raise ValueError('The short offset distance is larger than the long offset distance')


    def apply_filter(optim, trace):
    ''' apply filter
    '''
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
        length = 100
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

        # pick the first arrival and mute the late or early arrival
        if mute_late_arrival:

            # decide mute late or early arrivals W.R.T first break
            if mute_late_window > 0.0:
                mutetype = 'late'
            else:
                mutetype = 'early'
            
            # mute window in second
            mute_window = abs(mute_late_window)
            
            # use the firsy arrival in the header
            if use_first_break_in_su_header:
                itrace = 0
                pick = np.zeros(len(trace))
                for tr in trace:
                    pick[itrace] = int((tr.stats.su.trace_header.mute_time_start_time_in_ms/1000 + mute_window)/dt)
                    itrace +=1
            # pick the first arrival in a brutal manner 
            else:
                pick = brutal_picker(trace) + np.ceil(mute_window/dt)

            itrace = 0
            for tr in trace:
                # set mute window
                itmin = int(pick[itrace] - length/2)
                itmax = int(itmin + length)
                # apply mute
                mute_arrival(tr, itmin, itmax, mutetype, nt, length)
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
        
        trace  = su2array(trace)
        threds = 0.001 * np.max(abs(trace), axis=-1)

        pick = [(abs(trace[i,:]) > threds[i]).argmax(axis=-1) for i in range(trace.shape[0])]

        return np.array(pick)



    def mute_offset(trace, mute_dist, mutetype):
        ''' mutes traces according to offset
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
        ''' applies tapered mask to record section, muting early or late arrivals
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
        ''' constructs tapered mask that can be applied to trace to
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


    def source_wavelet_process(stf, stf_dt=0.001, use_dt = 0.001, use_nt=4001, shift = 0.0, lowpass = 0, highpass = 50, taper_beg=0.005):
        ''' source wavelet process
        '''
        # convert to su
        stf = array2su(1, stf_dt, stf)
        t0 = stf[0].stats.starttime 
        stf.resample(1.0/use_dt)
        stf.detrend('demean')
        stf.detrend('linear')
        # filter
        if lowpass == 0:
            stf.filter('lowpass', freq=highpass)#, corners=4, zerophase=True)
        elif 0 < lowpass < highpass:
            stf.filter('bandpass', freqmin=lowpass, freqmax=highpass)#, corners=4, zerophase=True)
        else:
            pass
        stf.trim(starttime=t0, endtime=t0 + use_dt*(use_nt-1), 
                pad=True, nearest_sample=True, fill_value=0.)
        #stf[0].data[:] = smooth1d(stf[0].data[:], window_len = 5)
        # apply shift
        if shift > 0:
            shift = int(shift//use_dt)
            stf[0].data[0:-1-shift] =  stf[0].data[shift:-1]
        
        stf.taper(taper_beg, type='hann', side='left')

        stf = stf[0].data[:]
        if stf.size != use_nt:
            ValueError('wrong size of the source wavelet')

        return stf