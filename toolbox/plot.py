###############################################################################
#
# SWIT: Seismic Waveform Inversion Toolbox
#
# by Haipeng Li at USTC, haipengl@mail.ustc.edu.cn
# 
# June, 2021  
#
# Plot module
#
###############################################################################

import time
import matplotlib
matplotlib.use('Agg') # Must be before importing matplotlib.pyplot or pylab!

import numpy as np
import matplotlib.pyplot as plt
from multiprocessing import Pool

from tools import loadsu, add_su_header, convert_wavelet_su

### Plot acquisition geometry
def plot_geometry(simu):
    ''' Plot source and receiver acquisition geometry
    '''
    srcn  = simu.source.n
    srcxz = simu.source.xz
    recxz = simu.receiver.xz
    figpath = simu.system.homepath + 'figures/'

    # Acquisition geometry
    fig = plt.figure()
    ax  = fig.add_subplot(111)
    for isrc in range(srcn):
        if np.mod(isrc,1) == 0:
            plt.scatter(recxz[isrc,:,0], 0 * recxz[isrc,:, 1] + 1 * isrc + 1, c = 'green', marker='o', s = 2)
            plt.scatter(srcxz[isrc, 0],  0 * srcxz[isrc,   1] + 1 * isrc + 1, c = 'red',   marker='^', s = 6)

    ax.set_ylim(0,  srcn+1)
    plt.xlabel('Distance-x (m)', fontsize=12)
    plt.ylabel('Shot #', fontsize=12)
    plt.title('Acquisition Geometry, %d sources\n' % (srcn), fontsize=14)
    plt.savefig(figpath + 'Acquisition-Geometry.png', dpi=300)
    plt.close()


### Plot source time function (STF).
def plot_stf(simu, isrc=1, stf_type='obs', t_end=1.0):
    ''' Plot source time function in the time and Frequency domanin
    '''

    if isrc < 1 or isrc > simu.source.n: 
        raise ValueError('isrc exceeds source range: 1~%d.\n'%simu.source.n)
    else:
        ISRC = isrc  
        isrc = isrc - 1
        
    # set data for plot
    stf_time = simu.source.wavelet[isrc, :]
    stf_spectrum = np.abs(np.fft.fft(stf_time))
    freqs = np.fft.fftfreq(len(stf_time), simu.model.dt)
    idx = np.argsort(freqs)
    idx = idx[int(len(idx) / 2):]

    nt = np.array(t_end/simu.model.dt)
    nt = nt.astype(int)

    # t0, t1 in ms; relative amps
    WAVELET_EXTENT = (0, simu.model.t[nt+1], -1.2, 1.2)
    SPECTRUM_EXTENT = (0, 50, 0, 1.2)  # f0, f1, p0, p1

    # Figure
    fig = plt.figure()
    # subplot1: time domain
    ax1 = fig.add_subplot(2, 1, 1)
    ax1.xaxis.set_label_text('Time (s)', fontsize=12)
    ax1.yaxis.set_label_text('Amplitude', fontsize=12)
    ax1.set_title('Source wavelet (normalized)- source %d - %s' % (ISRC, stf_type), fontsize=16)
    ax1.axis(WAVELET_EXTENT)
    ax1.plot(simu.model.t[0:nt+1], stf_time[0:nt+1] / abs(stf_time[0:nt+1]).max(), 'g-')
    # subplot2: frequency domain
    ax2 = fig.add_subplot(2, 1, 2)
    ax2.xaxis.set_label_text('Frequency (Hz)', fontsize=12)
    ax2.yaxis.set_label_text('Norm. Amplitude', fontsize=12)
    ax2.set_title('Amplitude spectrum', fontsize=16)
    ax2.axis(SPECTRUM_EXTENT)
    ax2.fill(freqs[idx], stf_spectrum[idx] / abs(stf_spectrum).max(), 'c')
    ax2.plot(freqs[idx], stf_spectrum[idx] / abs(stf_spectrum).max(), 'b')
    fig.tight_layout()
    plt.savefig(simu.system.homepath + 'figures/STF-%s-src%d.png' % (stf_type, ISRC), dpi=300)
    plt.close()


### Plot material (2D): velocity, gradient
def plot_model2D(simu, data, vmin, vmax, filename, colormap = 'jet'):
    ''' Plot model material (2D), e.g. vp, vs, and rho.
    '''
    xx = simu.model.xx
    zz = simu.model.zz
    figpath = simu.system.homepath + 'figures/model/'
    figaspect = simu.system.figaspect

    ## Figure: model_2D.
    fig = plt.figure(figsize=(10, 5))
    ax = fig.add_subplot(111)

    if colormap in ['my_seismic_cmap']:
        plotopts = {
            'vmin': vmin,
            'vmax': vmax,
            'cmap': my_seismic_cmap(),  # my colormap
            'extent': (xx[0], xx[-1], zz[-1], 0)
        }
    elif colormap in ['seismic']: 
        plotopts = {
            'vmin': vmin,
            'vmax': vmax,
            'cmap': plt.cm.bwr,  # my colormap
            'extent': (xx[0], xx[-1], zz[-1], 0)
        }
    else:
        plotopts = {
            'vmin': vmin,
            'vmax': vmax,
            'cmap': plt.cm.jet,  # my colormap
            'extent': (xx[0], xx[-1], zz[-1], 0)
        }


    im = ax.imshow(data, **plotopts)
    fig.colorbar(im, shrink=0.5, extend='both')
    ax.xaxis.set_label_text('Distance (m)')
    ax.yaxis.set_label_text('Depth (m)')
    ax.set_title(filename)
    ax.set_aspect(figaspect)
    plt.savefig(figpath + filename + '.png', dpi=300)
    plt.close()


def plot_trace(simu, filename, simu_type='syn', suffix='', src_space=5, trace_space=5, scale=1.0, color='r', plot_dx=2000):
    ''' Plot trace for SU stream data.
    '''

    # get prepared
    nt = simu.model.nt
    dt = simu.model.dt 
    srcn = simu.source.n
    srcn = simu.source.n
    nproc = simu.system.mpiproc
    homepath = simu.system.homepath
    figpath = homepath + 'figures/waveform/'

    src = list(range(0, srcn))
    if src_space > 1:
        src = src[0:-1:src_space]

    # submit plot
    pool = Pool(nproc)
    for isrc in range(srcn):
        datapath = homepath + 'data/' + simu_type + '/src%d_sg%s.su'%(isrc+1, suffix)
        figname = figpath + 'shot%03d-' % isrc + filename + '.png'
        pool.apply_async(plot_trace_serial, args=(datapath, figname, trace_space, scale, color, 
                         plot_dx,nt, dt, isrc, 'pressure', ))
        time.sleep(0.001)
    pool.close()
    pool.join()



def plot_trace_serial(datapath, figname, trace_space, scale, color, plot_dx, nt, dt, isrc, comp):
    
    trace = loadsu(datapath)

    trace = add_su_header(trace, nt, dt, isrc, comp)

    offset_min = trace[0].stats.distance
    offset_max = trace[-1].stats.distance

    fig = plt.figure(figsize=(10, 6))
    trace[0:-1:trace_space].plot(
        # recordstart = t1, recordlength = t2,
        type='section', scale=scale, time_down=True,
        offset_min=offset_min-200, offset_max=offset_max+200, plot_dx=plot_dx,
        fillcolors=(color, ''),  morm_method='trace', method='full',
        linewidth=0.25, grid_width=0.5, alpha=1.0, dpi=200, fig=fig)
    plt.savefig(figname, dpi=200)
    plt.close()



def plot_wavelet(simu, wavelet, filename, scale=1.0, color='r', plot_dx=1000, t_end = 1.0):
    ''' plot wavelet
    '''
    dt = simu.model.dt
    srcx = simu.source.xz[:,0]
    wavelet = convert_wavelet_su(dt, simu.source.wavelet, srcx)


    figpath = simu.system.homepath + 'figures/'

    fig = plt.figure(figsize=(10, 6))
    wavelet.plot(type='section', scale=scale, time_down=True, plot_dx=plot_dx,
            recordstart = 0., recordlength = t_end,
            fillcolors=(color, ''),  morm_method='trace', method='full',
            linewidth=0.25, grid_width=0.5, alpha=1.0, dpi=200, fig=fig)

    plt.savefig(figpath + filename + '.png', dpi=300)
    plt.close()



def plot_misfit(simu, misfit, mistype):
    ''' plot misfit
    '''
    figpath = simu.system.homepath + 'figures/'

    fig = plt.figure(figsize=(10, 6))
    plt.plot(misfit, marker='o')
    plt.xlabel('iteration', fontsize=14)
    plt.ylabel('misfit %s' % mistype, fontsize=14)
    plt.savefig(figpath + 'misfit_%s.png' % mistype, dpi=300)
    plt.close()


def plot_inv_scheme(simu, optim, inv_scheme):
    ''' plot many outputs
    '''
    nx = simu.model.nx
    nz = simu.model.nz
    it = optim.iter
    vpmin = optim.vpmin
    vpmax = optim.vpmax

    vp = simu.model.vp
    grad = inv_scheme['g_now']
    dire = inv_scheme['d_now']

    grad_caxis = np.max(abs(grad)) * 0.9
    dirc_caxis = np.max(abs(dire)) * 0.9

    plot_model2D(simu, vp.T, vpmin, vpmax, 'vp-%03d' % it)
    plot_model2D(simu, grad.reshape(nx, nz).T, -grad_caxis, grad_caxis, 'grad-%03d' % it, colormap = 'seismic')
    plot_model2D(simu, dire.reshape(nx, nz).T, -dirc_caxis, dirc_caxis, 'dire-%03d' % it, colormap = 'seismic')

    if optim.iter == 1 :
        plot_trace(simu, 'syn-initial-model-proc', simu_type = 'syn', suffix='_proc', src_space=1, trace_space=5, scale=0.8, color='k')
    elif optim.iter == optim.maxiter:
        data_misfit = np.loadtxt('./outputs/misfit_data.dat')
        data_misfit = data_misfit / data_misfit[0]
        plot_misfit(simu, data_misfit, 'data')
        plot_trace(simu, 'syn-final-model-proc', simu_type = 'syn', suffix='_proc', src_space=1, trace_space=5, scale=0.8, color='b')
    else:
        pass


def my_seismic_cmap():
    ''' my seismic cmap
    '''
    
    cdict = {'red': ((0.0, 0.0, 0.0),
                    (0.1, 0.5, 0.5),
                    (0.2, 0.0, 0.0),
                    (0.4, 0.2, 0.2),
                    (0.6, 0.0, 0.0),
                    (0.8, 1.0, 1.0),
                    (1.0, 1.0, 1.0)),
            'green':((0.0, 0.0, 0.0),
                    (0.1, 0.0, 0.0),
                    (0.2, 0.0, 0.0),
                    (0.4, 1.0, 1.0),
                    (0.6, 1.0, 1.0),
                    (0.8, 1.0, 1.0),
                    (1.0, 0.0, 0.0)),
            'blue': ((0.0, 0.0, 0.0),
                    (0.1, 0.5, 0.5),
                    (0.2, 1.0, 1.0),
                    (0.4, 1.0, 1.0),
                    (0.6, 0.0, 0.0),
                    (0.8, 0.0, 0.0),
                    (1.0, 0.0, 0.0))}
    my_cmap = matplotlib.colors.LinearSegmentedColormap('my_colormap',cdict,256)

    return my_cmap