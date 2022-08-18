###############################################################################
#
# SWIT: Seismic Waveform Inversion Toolbox
#
# by Haipeng Li at USTC, haipengl@mail.ustc.edu.cn
# 
# June, 2021  
#
# Field data module
#
###############################################################################


from multiprocessing import Pool, cpu_count
import numpy as np
import obspy

from plot import plot_wavelet
from tools import array2su, get_offset, get_su_parameter, loadsu, savesu, smooth1d, su2array


def load_field_data(simu, data_list, x_beg, x_end):
    ''' Load field data (pre-processed), vertical geophone
    '''
    
    nproc    = np.min((simu.system.mpiproc, cpu_count()//2))
    homepath = simu.system.homepath

    # get prepared
    srcn = len(data_list)
    dt = simu.model.dt
    nt = simu.model.nt
    
    print('Field data: loading data and adjusting to current acquisition')

    pool = Pool(nproc)
    receiver_range = [pool.apply_async(load_field_data_serial, (homepath, data_list[isrc], isrc, x_beg, x_end, dt, nt, )) for isrc in range(srcn)]
    pool.close()
    receiver_range = np.array([p.get() for p in receiver_range])

    pool.join()

    simu.receiver.rec_beg = receiver_range[:,0]
    simu.receiver.rec_end = receiver_range[:,1]

    print('Field data: %d shotgathers in total' %(srcn))
    print('Field data: %d time steps, %.2f ms samping interval\n'%(nt, dt*1000))
    print('*****************************************************\n')



def load_field_data_serial(homepath, loadfile, isrc, x_beg, x_end, dt, nt):

        # set path
        savefile = homepath + 'data/obs/src%d_sg.su' %(isrc + 1)

        # load data
        trace = loadsu(loadfile)
        _, _, dt_data = get_su_parameter(trace)

        # select traces in calculate domain
        trace = obspy.Stream(traces = [tr for tr in trace
                if  tr.stats.su.trace_header.group_coordinate_x - x_beg >= 0.0
                and tr.stats.su.trace_header.group_coordinate_x - x_end <= 0.0])
        
        # modify header to adjust to current acquisition geometry and for obspy plotting
        for tr in trace:
            tr.stats.su.trace_header.group_coordinate_x  = int(tr.stats.su.trace_header.group_coordinate_x  - x_beg)
            tr.stats.su.trace_header.source_coordinate_x = int(tr.stats.su.trace_header.source_coordinate_x - x_beg)
            tr.stats.distance = tr.stats.su.trace_header.group_coordinate_x - tr.stats.su.trace_header.source_coordinate_x
        
        # interpolate if necessary
        if dt_data != dt:
            for tr in trace:
                tr.resample(1.0/dt)

        # ensure the same time steps
        for tr in trace:
            t0 = tr.stats.starttime 
            tr.trim(starttime=t0, endtime=t0 + dt*(nt-1), pad=True, nearest_sample=True, fill_value=0.0)
        
        # save file to the working folder
        savesu(savefile, trace)

        # record the receiver range
        x0 = trace[ 0].stats.su.trace_header.group_coordinate_x
        x1 = trace[-1].stats.su.trace_header.group_coordinate_x

        return np.array([x0, x1])


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
    stf[0].data[:] = smooth1d(stf[0].data[:], span = 5)
    # apply shift
    if shift > 0:
        shift = int(shift//use_dt)
        stf[0].data[0:-1-shift] =  stf[0].data[shift:-1]
    
    stf.taper(taper_beg, type='hann', side='left')

    stf = stf[0].data[:]
    if stf.size != use_nt:
        ValueError('wrong size of the source wavelet')

    return stf


def source_inversion(simu, inv_offset=10000):
    ''' source signature inversion by (Parrat, 1999)
    '''
    # get prepared
    homepath = simu.system.homepath
    stf_now = simu.source.wavelet
    srcx = simu.source.xyz[:,0]
    recx = simu.receiver.xyz[:,:,0]
    srcn = simu.source.n
    dt = simu.model.dt

    stf_inv = np.zeros_like(stf_now)

    for isrc in range(srcn):

        # load data and set offset
        obs = loadsu(homepath + 'data/obs/src%d_sg_proc.su'%(isrc+1))
        syn = loadsu(homepath + 'data/syn/src%d_sg_proc.su'%(isrc+1))
        if not np.array_equal(get_offset(obs), get_offset(syn)):
            raise ValueError("offset not consistant.")
        offset = get_offset(syn)
        obs = su2array(obs)
        syn = su2array(syn)

        # select data for source inversion
        rec_used = np.argwhere(abs(offset) < inv_offset)
        n1 = int(rec_used[0])
        n2 = int(rec_used[-1])
        data_obs = np.squeeze(obs[n1:n2, :]).T
        data_syn = np.squeeze(syn[n1:n2, :]).T
        Do = np.fft.fft(data_obs, axis=0)  # frequency domain
        Dm = np.fft.fft(data_syn, axis=0)  # frequency domain
        
        # current source wavelet
        src = np.squeeze(stf_now[isrc, :])
        S = np.fft.fft(np.squeeze(src), axis=0) # frequency domain

        # check
        if abs(np.sum(Dm * np.conj(Dm))) == 0:
            raise ValueError("No trace for source inversion, check for the reason.")
        
        # source inversion
        A = np.sum(np.conj(Dm)*Do, axis=1) / np.sum(Dm * np.conj(Dm), axis=1)
        temp = np.real(np.fft.ifft(A*S[:]))
        temp = temp / np.max(abs(temp))
        stf_inv[isrc, :] = temp

    # process the source wavelet
    stf_now = convert_wavelet_su(dt, stf_now, srcx)
    stf_inv = convert_wavelet_su(dt, stf_inv, srcx)

    # save source wavelet plots
    plot_wavelet(simu, stf_now, 'stf_now', scale=1.0, color='k', plot_dx=5000, t_end = 1.0)
    plot_wavelet(simu, stf_inv, 'stf_inv', scale=1.0, color='r', plot_dx=5000, t_end = 1.0)
   
    # save source wavelet data
    np.savetxt(homepath+'outputs/stf_now.dat', su2array(stf_now))
    np.savetxt(homepath+'outputs/stf_inv.dat', su2array(stf_inv))

    print('Source inversion finished\n')


def set_frequency_band(zmax, hmax, f_begin, band_num):
    ''' Multi-scale implemation based on Boonyasiriwat et al, 2009.
        f(n+1)    = f(n) / alpha_min,       where 
        alpha_min = z/sqrt( h**2 + z**2),   and 
        z is maximum depth to be imaged,    and
        h is maximum half-offset.
    '''

    alpha_min = zmax / np.sqrt(zmax**2 + hmax**2)
    fre_band = np.zeros(band_num)
    fre_band[0] = f_begin
    for iband in range(band_num):
        if iband > 0:
            fre_band[iband] = fre_band[iband-1] / alpha_min

    # set two decimal
    fre_band = np.around(fre_band, decimals=2)

    return fre_band


def convert_wavelet_su(dt, wavelet, srcx):
    ''' Convert wavelet array to the SU streame
    '''
    srcn  = np.size(wavelet, 0)

    wavelet_su = array2su(srcn, dt, wavelet)

    ishot = 0
    for iwvlt in wavelet_su:
        iwvlt.stats.distance = srcx[ishot]
        ishot+=1


    return wavelet_su
