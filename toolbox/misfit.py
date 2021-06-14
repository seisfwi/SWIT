###############################################################################
#
# SWIT: Seismic Waveform Inversion Toolbox
#
# by Haipeng Li at USTC, haipengl@mail.ustc.edu.cn
# 
# June, 2021  
#
# Misfit module, some of codes are from: https://github.com/rmodrak/seisflows 
#                                      & https://github.com/pysit/pysit
#
###############################################################################


from multiprocessing import Pool
import numpy as np
from scipy import fftpack
from scipy.signal import hilbert

from tools import get_su_parameter, loadsu, su2array


def misfit(simu, misfit_type):
    ''' Calculate the misfit function
    '''
    homepath = simu.system.homepath    
    nproc = simu.system.mpiproc
    srcn = simu.source.n

    pool = Pool(nproc)
    misfit = [pool.apply_async(misfit_serial, (homepath, isrc, misfit_type, ))  for isrc in range(srcn)]
    pool.close()
    misfits = [p.get() for p in misfit]
    pool.join()

    return np.sum(misfits)


def misfit_serial(homepath, isrc, misfit_type):
    ''' Calculate the misfit function for a single shot
    '''

    obs = loadsu(homepath + 'data/obs/src%d_sg_proc.su'%(isrc+1))
    syn = loadsu(homepath + 'data/syn/src%d_sg_proc.su'%(isrc+1))

    if get_su_parameter(obs) != get_su_parameter(syn):
        raise ValueError('obs and syn are not consistent.')

    # Waveform difference L2-norm (Tarantola, 1984)
    if misfit_type.lower() in ['waveform']:
        rsd = misfit_waveform(obs, syn)
    
    # Envelope difference (Wu et al., 2014; Yuan et al., 2015)
    elif misfit_type.lower() in ['envelope']:
        rsd = misfit_envelope(obs, syn)

    # Cross correlation traveltime (Luo & Schuster, 1991; Tromp et al., 2005)
    elif misfit_type.lower() in ['traveltime']:
        rsd = misfit_traveltime(obs, syn)

    # Normalized global-correlation coefficient (Choi & Alkhalifah, 2012)
    elif misfit_type.lower() in ['globalcorrelation']:
        rsd = misfit_global_correlation(obs, syn)

    return rsd



def adjoint_source(simu, misfit_type):
    ''' Caculate the adjoint source
    '''

    homepath = simu.system.homepath    
    nproc = simu.system.mpiproc
    srcn = simu.source.n

    pool = Pool(nproc)
    adj = [pool.apply_async(adjoint_source_serial, (homepath, isrc, misfit_type, ))  for isrc in range(srcn)]
    pool.close()
    adjs = [p.get() for p in adj]
    pool.join()

    return np.array(adjs)



def adjoint_source_serial(homepath, isrc, misfit_type):
    ''' Caculate the adjoint source for a single source
    '''

    obs = loadsu(homepath + 'data/obs/src%d_sg_proc.su'%(isrc+1))
    syn = loadsu(homepath + 'data/syn/src%d_sg_proc.su'%(isrc+1))

    if get_su_parameter(obs) != get_su_parameter(syn):
        raise ValueError('obs and syn are not consistent.')

    # Waveform difference L2-norm (Tarantola, 1984)
    if misfit_type.lower() in ['waveform']:
        adj = adjoint_source_waveform(obs, syn)

    # Envelope difference (Wu et al., 2014; Yuan et al., 2015)
    elif misfit_type.lower() in ['envelope']:
        adj = adjoint_source_envelope(obs, syn)

    # Cross correlation traveltime (Luo & Schuster, 1991; Tromp et al., 2005)
    elif misfit_type.lower() in ['traveltime']:
        adj = adjoint_source_traveltime(obs, syn)

    # Normalized global-correlation coefficient (Choi & Alkhalifah, 2012)
    elif misfit_type.lower() in ['globalcorrelation']:
        adj = adjoint_source_global_correlation(obs, syn)

    return adj


def adjoint_source_waveform(obs, syn):
    ''' Waveform difference L2-norm (Tarantola, 1984)
    '''
    # parameters
    recn, nt, _ = get_su_parameter(obs)
    adj = np.zeros((recn, nt))

    for irec in range(recn):
        obs_trace = obs[irec].data
        syn_trace = syn[irec].data
        obs_norm = np.linalg.norm(obs_trace, ord=2)
        syn_norm = np.linalg.norm(syn_trace, ord=2)

        adj_trace = np.zeros(nt)

        if obs_norm > 0. and syn_norm > 0.:
            adj_trace = syn_trace - obs_trace

        adj[irec,:] = adj_trace

    return adj


def adjoint_source_envelope(obs, syn):
    ''' Envelope difference (Wu et al., 2014; Yuan et al., 2015)
    '''
    # parameters
    p = 2.0            # envelope_power
    recn, nt, _ = get_su_parameter(obs)
    adj = np.zeros((recn, nt))

    for irec in range(recn):
        obs_trace = obs[irec].data
        syn_trace = syn[irec].data
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
            adj_trace += (-hilbert(denvelope_ddataH * rsd_envelope, axis=0)).imag

        adj[irec,:] = adj_trace

    return adj



def adjoint_source_traveltime(obs, syn):
    ''' Cross correlation traveltime (Luo & Schuster, 1991; Tromp et al., 2005)
    '''
    # parameters
    recn, nt, dt = get_su_parameter(obs)
    adj = np.zeros((recn, nt))

    # compute the cross-correlation
    ccmax = cross_correlate_max(obs, syn)

    for irec in range(recn):
        obs_trace = obs[irec].data
        syn_trace = syn[irec].data
        obs_norm = np.linalg.norm(obs_trace, ord=2)
        syn_norm = np.linalg.norm(syn_trace, ord=2)

        adj_trace = np.zeros(nt)
        adj_trace[1:-1] = (syn_trace[2:] - syn_trace[0:-2])/(2.*dt)

        if obs_norm > 0. and syn_norm > 0. and np.sum(abs(adj_trace)) > 0.:
            adj_trace *= 1./(np.sum(adj_trace*adj_trace)*dt)
            adj_trace *= (ccmax[irec]-nt+1)*dt
        else:
            adj_trace *= 0.

        adj[irec, :] = adj_trace

    return adj



def adjoint_source_global_correlation(obs, syn):
    ''' Normalized global-correlation coefficient (Choi & Alkhalifah, 2012)
    '''
    # parameters
    recn, nt, dt = get_su_parameter(obs)
    adj = np.zeros((recn, nt))

    for irec in range(recn):
        obs_trace = obs[irec].data
        syn_trace = syn[irec].data
        obs_norm = np.linalg.norm(obs_trace, ord=2)
        syn_norm = np.linalg.norm(syn_trace, ord=2)

        adj_trace = np.zeros(nt)
        if obs_norm > 0. and syn_norm > 0.:
            obs_trace_norm = obs_trace / obs_norm
            syn_trace_norm = syn_trace / syn_norm
            adj_trace = 1.0/syn_norm * (syn_trace_norm * np.corrcoef(syn_trace_norm, obs_trace_norm)[0,1] - obs_trace_norm)
        else:
            adj_trace *= 0.

        adj[irec, :] = adj_trace

    return adj


def misfit_waveform(obs, syn):
    ''' Waveform difference L2-norm (Tarantola, 1984)
    '''
    # parameters
    recn, nt, dt = get_su_parameter(obs)
    rsd = np.zeros(1)

    for irec in range(recn):
        obs_trace = obs[irec].data
        syn_trace = syn[irec].data
        obs_norm = np.linalg.norm(obs_trace, ord=2)
        syn_norm = np.linalg.norm(syn_trace, ord=2)

        rsd_trace = np.zeros(nt)

        if obs_norm > 0. and syn_norm > 0.:
            rsd_trace = syn_trace - obs_trace
            rsd += np.sqrt(np.sum(rsd_trace*rsd_trace*dt))

    return rsd


def misfit_envelope(obs, syn):
    ''' Envelope difference (Wu et al., 2014; Yuan et al., 2015)
    '''
    # parameters
    p = 2.0            # envelope_power
    recn, nt, dt = get_su_parameter(obs)
    rsd = np.zeros(1)

    for irec in range(recn):
        obs_trace = obs[irec].data
        syn_trace = syn[irec].data
        obs_norm = np.linalg.norm(obs_trace, ord=2)
        syn_norm = np.linalg.norm(syn_trace, ord=2)

        if obs_norm > 0. and syn_norm > 0.:
            syn_Hilbert = hilbert(syn_trace, axis=0).imag
            obs_Hilbert = hilbert(obs_trace, axis=0).imag

            syn_envelope = syn_trace**2.0 + syn_Hilbert**2.0
            obs_envelope = obs_trace**2.0 + obs_Hilbert**2.0
            rsd_envelope = syn_envelope**(p/2.0) - obs_envelope**(p/2.0)

            rsd += np.linalg.norm(rsd_envelope)**2

    return rsd * dt


def misfit_traveltime(obs, syn):
    ''' Cross correlation traveltime (Luo & Schuster, 1991; Tromp et al., 2005)
    '''
    # parameters
    recn, nt, dt = get_su_parameter(obs)
    rsd = np.zeros(1)

    # compute the cross-correlation in an awkward paralle scheme
    ccmax =  cross_correlate_max(obs, syn)

    for irec in range(recn):
        obs_trace = obs[irec].data
        syn_trace = syn[irec].data
        obs_norm = np.linalg.norm(obs_trace, ord=2)
        syn_norm = np.linalg.norm(syn_trace, ord=2)

        adj_trace = np.zeros(nt)
        adj_trace[1:-1] = (syn_trace[2:] - syn_trace[0:-2])/(2.*dt)

        if obs_norm > 0. and syn_norm > 0. and np.sum(abs(adj_trace)) > 0.:
            rsd += 0.5 * np.power(ccmax[irec]-nt+1, 2)*dt

    return rsd



def misfit_global_correlation(obs, syn):
    ''' Normalized global-correlation coefficient (Choi & Alkhalifah, 2012)
    '''
    # get parameters
    recn, nt, dt = get_su_parameter(obs)
    rsd = np.zeros(1)

    for irec in range(recn):
        obs_trace = obs[irec].data
        syn_trace = syn[irec].data
        obs_norm = np.linalg.norm(obs_trace, ord=2)
        syn_norm = np.linalg.norm(syn_trace, ord=2)
        if obs_norm > 0. and syn_norm > 0.:
            obs_trace_norm = obs_trace/obs_norm
            syn_trace_norm = syn_trace/syn_norm
            rsd += - np.corrcoef(syn_trace_norm, obs_trace_norm)[0,1]

    return rsd



def cross_correlate_max(obs, syn):
    ''' calculate the cross-correlation lag between two traces
    '''

    obs_data = su2array(obs)
    syn_data = su2array(syn)
    nt = np.size(obs_data, -1)

    a =   fftpack.fft(obs_data)
    b = - fftpack.fft(syn_data).conjugate()

    cc = np.argmax(np.abs(fftpack.ifft(a*b)), -1) - 1
    cc[np.where(cc<nt//2)] = cc[np.where(cc<nt//2)] + nt

    return cc
