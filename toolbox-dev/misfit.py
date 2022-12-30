###############################################################################
#
# SWIT v1.1: Seismic Waveform Inversion Toolbox
#
#   A Python package for seismic waveform inversion
#   By Haipeng Li at USTC & Stanford
#   Email: haipengl@mail.ustc.edu.cn, haipeng@stanford.edu 
#
#   Misfit module for calculating misfit functions and adjoint sources
#   Some of codes are from:  https://github.com/rmodrak/seisflows 
#                          & https://github.com/pysit/pysit
#
###############################################################################

import os

import numpy as np
from scipy import fftpack
from scipy.signal import hilbert
from tools import load_waveform_data, save_float


def calculate_adjoint_misfit_is(isrc, data_path, dt, nt, misfit_type):
    ''' Caculate the adjoint source for a single source

    Parameters
    ----------
        isrc : int
            source index
        data_path : str
            working directory
        dt : float
            sampling rate
        nt : int
            number of time samples
        misfit_type : str
            misfit type, 'waveform', 'envelope', 'crosscorrelation' or 'globalcorrelation'
    '''

    # load observed and synthetic data
    obs_path = os.path.join(data_path, 'data/obs/src{}/sg_processed'.format(isrc+1))
    syn_path = os.path.join(data_path, 'data/syn/src{}/sg_processed'.format(isrc+1))

    # load observed data
    obs, obs_dt = load_waveform_data(obs_path, nt)
    if obs_dt is None:
        obs_dt = dt

    # load synthetic data
    syn, syn_dt = load_waveform_data(syn_path, nt)
    if syn_dt is None:
        syn_dt = dt

    # check the consistency of the obs and syn data
    if obs_dt != syn_dt or obs_dt != dt:
        raise ValueError('The sampling rate of the obs and syn are not consistent.')

    if obs.shape != syn.shape:
        raise ValueError('The shape of the obs and syn are not consistent.')

    if misfit_type.lower() in ['waveform']:
        adj, rsd = misfit_waveform(obs, syn, dt)

    elif misfit_type.lower() in ['envelope']:
        adj, rsd = misfit_envelope(obs, syn, dt)

    elif misfit_type.lower() in ['crosscorrelation']:
        adj, rsd = misfit_crosscorrelation(obs, syn, dt)

    elif misfit_type.lower() in ['globalcorrelation']:
        adj, rsd = misfit_globalcorrelation(obs, syn, dt)
    
    else:
        msg = 'Only support the following misfit types: waveform, envelope, traveltime, globalcorrelation'
        err = 'The misfit type {} is not supported.'.format(misfit_type)
        raise ValueError(msg + '\n' + err)
    
    # write the adjoint source to file
    adj_path = os.path.join(data_path, 'config/wavelet/src{}_adj.bin'.format(isrc+1))
    save_float(adj_path, adj)

    # return summed misfit function over all traces for a single source
    return np.sum(rsd)


def misfit_waveform(obs, syn, dt):
    ''' Waveform difference L2-norm (Tarantola, 1984)
    '''

    # parameters
    rec_num, nt = obs.shape
    
    # initialize the adjoint source and misfit
    adj = np.zeros((rec_num, nt))
    rsd = np.zeros(rec_num)

    # loop over all traces
    for irec in range(rec_num):
        obs_trace = obs[irec,:]
        syn_trace = syn[irec,:]
        obs_norm = np.linalg.norm(obs_trace, ord=2)
        syn_norm = np.linalg.norm(syn_trace, ord=2)
        if obs_norm > 0. and syn_norm > 0.:
            # calculate the adjoint source
            adj[irec,:] = syn_trace - obs_trace

            # calculate the misfit
            rsd[irec] = np.sum(adj[irec,:]*adj[irec,:]*dt)

    # return adjoint source and misfit
    return adj, rsd


def misfit_envelope(obs, syn, dt):
    ''' Envelope difference (Wu et al., 2014; Yuan et al., 2015)
    '''
    # envelope_power
    p = 2.0
    
    # parameters
    rec_num, nt = obs.shape

    # initialize the adjoint source and misfit
    adj = np.zeros((rec_num, nt))
    rsd = np.zeros(rec_num)

    # loop over all traces
    for irec in range(rec_num):
        obs_trace = obs[irec,:]
        syn_trace = syn[irec,:]
        obs_norm = np.linalg.norm(obs_trace, ord=2)
        syn_norm = np.linalg.norm(syn_trace, ord=2)

        adj_trace = np.zeros(nt)
        rsd_envelope = np.zeros(nt)
        
        if obs_norm > 0. and syn_norm > 0.:
            syn_Hilbert = hilbert(syn_trace, axis=0).imag
            obs_Hilbert = hilbert(obs_trace, axis=0).imag

            syn_envelope = syn_trace**2.0 + syn_Hilbert**2.0
            obs_envelope = obs_trace**2.0 + obs_Hilbert**2.0
            rsd_envelope = syn_envelope**(p/2.0) - obs_envelope**(p/2.0)

            denvelope_ddata = p * syn_envelope**(p/2.0 - 1.0) * syn_trace
            adj_trace = denvelope_ddata * rsd_envelope

            denvelope_ddataH = p * syn_envelope**(p/2.0 - 1.0) * syn_Hilbert 
            
            # calculate the adjoint source
            adj_trace += (-hilbert(denvelope_ddataH * rsd_envelope, axis=0)).imag
            adj[irec,:] = adj_trace

            # calculate misfit
            rsd[irec] = np.linalg.norm(rsd_envelope)**2  * dt

    # return adjoint source and misfit
    return adj, rsd



def misfit_crosscorrelation(obs, syn, dt):
    ''' Cross correlation traveltime (Luo & Schuster, 1991; Tromp et al., 2005)
    '''
    # parameters
    rec_num, nt = obs.shape

    # initialize the adjoint source and misfit
    adj = np.zeros((rec_num, nt))
    rsd = np.zeros(rec_num)

    # compute the cross-correlation
    ccmax = cross_correlate_max(obs, syn, nt)

    # loop over all traces
    for irec in range(rec_num):
        obs_trace = obs[irec, :]
        syn_trace = syn[irec, :]
        obs_norm = np.linalg.norm(obs_trace, ord=2)
        syn_norm = np.linalg.norm(syn_trace, ord=2)

        adj_trace = np.zeros(nt)
        adj_trace[1:-1] = (syn_trace[2:] - syn_trace[0:-2])/(2.*dt)

        if obs_norm > 0. and syn_norm > 0. and np.sum(abs(adj_trace)) > 0.:
            # calculate adjoint source
            adj_trace *= 1./(np.sum(adj_trace*adj_trace)*dt)
            adj_trace *= (ccmax[irec]-nt+1)*dt
            
            # calculate misfit
            rsd[irec] = 0.5 * np.power(ccmax[irec]-nt+1, 2)*dt
        else:
            adj_trace *= 0.

        # fill the adjoint source
        adj[irec, :] = adj_trace

    # return adjoint source and misfit
    return adj, rsd



def misfit_globalcorrelation(obs, syn, dt):
    ''' Normalized global-correlation coefficient (Choi & Alkhalifah, 2012)
    '''
    # parameters
    rec_num, nt = obs.shape

    # initialize the adjoint source and misfit
    adj = np.zeros((rec_num, nt))
    rsd = np.zeros(rec_num)

    # loop over all traces
    for irec in range(rec_num):
        obs_trace = obs[irec, :]
        syn_trace = syn[irec, :]
        obs_norm = np.linalg.norm(obs_trace, ord=2)
        syn_norm = np.linalg.norm(syn_trace, ord=2)

        adj_trace = np.zeros(nt)
        if obs_norm > 0. and syn_norm > 0.:
            obs_trace_norm = obs_trace / obs_norm
            syn_trace_norm = syn_trace / syn_norm

            # calculate adjoint source
            adj_trace = 1.0/syn_norm * (syn_trace_norm * np.corrcoef(syn_trace_norm, obs_trace_norm)[0,1] - obs_trace_norm)

            # calculate misfit
            rsd[irec] = - np.corrcoef(syn_trace_norm, obs_trace_norm)[0,1] * dt
        else:
            adj_trace *= 0.

        # fill the adjoint source
        adj[irec, :] = adj_trace

    # return adjoint source and misfit
    return adj, rsd


def cross_correlate_max(obs, syn, nt):
    ''' Calculate the cross-correlation lag between two traces
    '''

    a =   fftpack.fft(obs)
    b = - fftpack.fft(syn).conjugate()

    cc = np.argmax(np.abs(fftpack.ifft(a*b)), -1) - 1
    cc[np.where(cc<nt//2)] = cc[np.where(cc<nt//2)] + nt

    return cc
