###############################################################################
#
# SWIT: Seismic Waveform Inversion Toolbox
#
# by Haipeng Li at USTC, haipengl@mail.ustc.edu.cn
# 
# June, 2021  
#
# Workflow
#
########################################################################################

# import modules
import numpy as np
import base
from inversion import inversion, source_inversion
from plot import plot_geometry, plot_model2D, plot_stf, plot_trace
from preprocess import process_workflow
from solver import forward, source_wavelet
from tools import saveparjson, smooth2d, loadbinfloat32

# set frequency band
freq   = [ 2, 4, 6, 8, 10]
smooth = [10, 6, 4, 2,  2]

for istage in range(5):

    ### system setup
    homepath  = '/data1/SWIT/case6-overthrust-land-multiscale/stage-%d'%(istage+1)   # working path
    mpiproc   = 41                                     # mpi process for fd2dmpi
    figaspect = 1.0                                    # Figure aspect

    ### model setup
    nx,  nz      = [401,  101]                         # Grid number along x and z directions
    pml, fs      = [50,  True]                         # Grid number for PML layers (use a large one) and free surface
    dx,  dt,  nt = [25, 0.002, 4001]                   # Grid size, time interval, and time step

    # velocity and density
    vp_true = np.zeros((nx, nz))
    vp_init = np.zeros((nx, nz))
    rho_true = np.zeros((nx, nz))
    rho_init = np.zeros((nx, nz))

    # true model
    model_folder = '/home/haipeng/Nutstore Files/Nutstore/MyCodes/SWIT-1.0/examples/case6-overthrust-land-multiscale/'
    vp_true = np.loadtxt(model_folder + 'model/Overthrust_401_101_25m_True.dat')

    # initial model
    if istage == 0:
        vp_init = np.copy(vp_true)
        vp_init = smooth2d(vp_true, span=30)
        #vp_init[:,:] = np.linspace(1500,5500,nz)
        #np.savetxt(model_folder + 'model/Overthrust_401_101_25m_1D.dat', vp_init)
    else:
        vp_init = loadbinfloat32('/data1/SWIT/case6-overthrust-land-multiscale/stage-%d/outputs/velocity/vp-20.bin'%(istage)).reshape(nx, nz)    
        vp_init = smooth2d(vp_init, span = 5)


    # density models, (Gardner, 1974)
    rho_true = np.power(vp_true, 0.25) * 310 
    rho_init = np.power(vp_init, 0.25) * 310 

    ### sources setup 
    f0    = 5.0                                        # Dominant frequency in Hz
    srcxz = np.loadtxt(model_folder + 'model/source_coordinate.dat')
    srcn  = srcxz.shape[0]                             # Source number along x axis
    wavelet  = np.zeros((srcn, nt))                    # source wavelet
    for isrc in range(srcn):
        wavelet[isrc,:] = source_wavelet(nt, dt, f0, 'ricker')

    ### receivers setup
    temp = np.loadtxt(model_folder + 'model/receiver_coordinate.dat')
    recn  = temp.shape[0]                              # receiver number
    recxz = np.zeros((srcn, recn, 2))                  # receiver positions
    for isrc in range(srcn):
        recxz[isrc,:,0] = temp[:,0]                    # receiver x position (m)
        recxz[isrc,:,1] = temp[:,1]                    # receiver z position (m)


    ### inversion parameter
    misfit_type = 'Waveform'                           # 'Traveltime', 'Waveform', 'Globalcorrelation'
    scheme      = 'NLCG'                               # 'NLCG', 'LBFGS'
    maxiter     = 20                                   # maximum iteration number,  i.e., 20
    step_length = 0.01                                 # maximum update percentage, i.e., 0.05
    vpmax       = 5500                                 # maximum allowed velocity
    vpmin       = 1500                                 # minimum allowed velocity
    marine_or_land = 'Land'                          # 'Land' or 'Marine'

    # gradient postprocess
    grad_mute = 10                                     # mute source energy for 'Land', or water mask for 'Marine'
    grad_smooth = smooth[istage]                       # gradient smooth radius 

    # data filter
    fre_filter = 'Lowpass'                             # 'Bandpass', 'Lowpass', 'Highpass', 'None'
    fre_low  = freq[istage]                            # low  frequency corner (units: Hz)
    fre_high = 20                                      # high frequency corner (units: Hz)

    # mute later arrival
    mute_late_arrival = False                          # pick first break and mute later arrival
    mute_late_window = 1.0                             # mute time window (units: time)

    # mute offset
    mute_offset_short  = False                         # mute short-offset traces 
    mute_offset_long   = False                         # mute long-offset traces 
    mute_offset_short_dis = 500                        # mute short-offset distance (units: m)
    mute_offset_long_dis  = 5000                       # mute long-offset distance (units: m)

    # data normalize
    normalize = ['None']                               # 'Max-Trace', 'L1-Trace', 'L2-Trace', 'L1-Event', 'L2-Event', 'None'


    ### simulate setup 
    sys  = base.system(homepath, mpiproc, figaspect=figaspect)
    mod  = base.model(nx, nz, dx, dt, nt, fs, pml, vp_true, rho_true)
    src  = base.source(f0, srcn, srcxz, wavelet)
    rec  = base.receiver(recn, recxz)
    simu = base.simulate(mod, src, rec, sys)

    ### optimize setup 
    optim = base.optimize(misfit_type, scheme, maxiter, step_length, vpmax, vpmin, marine_or_land,
                    grad_mute, grad_smooth,
                    fre_filter, fre_low, fre_high, 
                    mute_late_arrival, mute_late_window, normalize,
                    mute_offset_short, mute_offset_long, 
                    mute_offset_short_dis, mute_offset_long_dis)

    ### Save parameter as json
    saveparjson(simu, optim)

    ### Plots
    plot_geometry(simu)
    plot_stf(simu, isrc=1,  stf_type='in-use', t_end = 2.0)
    plot_model2D(simu, vp_true.T, vpmin, vpmax, 'vp-obs', colormap = 'jet')
    plot_model2D(simu, vp_init.T, vpmin, vpmax, 'vp-ini', colormap = 'jet')

    ### forward data
    forward(simu, simu_type='obs', savesnap=0)
    process_workflow(simu, optim, simu_type='obs')
    plot_trace(simu, 'obs',      simu_type='obs', suffix='',      src_space=1, trace_space=5, scale = 0.8, color='k')
    plot_trace(simu, 'obs-proc', simu_type='obs', suffix='_proc', src_space=1, trace_space=5, scale = 0.8, color='k')

    ### begin inversion
    inversion(simu, optim, {'vp':vp_init,'rho':rho_init})

    ### source inversion
    source_inversion(simu, inv_offset=20000)
