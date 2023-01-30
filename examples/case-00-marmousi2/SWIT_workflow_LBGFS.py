###############################################################################
#
# SWIT: Seismic Waveform Inversion Toolbox
#
# by Haipeng Li at USTC, haipengl@mail.ustc.edu.cn
# 
# May, 2022  
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
from tools import saveparjson, smooth2d


modelpath = '/data/haipeng/SWIT-1.0/examples/case-01-marmousi2/model/'   # Model path, i.e., velocity, source & receiver coordinates
homepath  = '/data/haipeng/SWIT-1.0/examples/case-01-marmousi2/LBFGS/'   # Working path
mpiproc   = 49                                     # MPI process for fd2dmpi
figaspect = 1.0                                    # Figure aspect

### model setup
nx,  nz      = [481,  141]                         # Grid number along x and z directions
pml, fs      = [50,  True]                         # Grid number for PML layers (use a large one) and free surface
dx,  dt,  nt = [25, 0.002, 4001]                   # Grid size, time interval, and time step

# velocity and density
vp_true = np.zeros((nx, nz))
vp_init = np.zeros((nx, nz))
rho_true = np.zeros((nx, nz))
rho_init = np.zeros((nx, nz))

# true model
vp_true = np.loadtxt(modelpath + 'Marmousi_481_141_25m.dat')

# initial model
vp_init = np.copy(vp_true)
vp_init[:,20:] = smooth2d(vp_true[:,20:], span = 10)
np.savetxt(modelpath + 'Marmousi_481_141_25m_smooth10.dat', vp_init)

# density models, (Gardner, 1974)
rho_true = np.power(vp_true, 0.25) * 310
rho_init = np.power(vp_true, 0.25) * 310

### sources setup 
f0    = 5.0                                              # Dominant frequency in Hz
srcxz = np.loadtxt(modelpath + 'source_coordinate.dat')  # Source coordinates
srcn  = srcxz.shape[0]                                   # Source number along x axis
wavelet  = np.zeros((srcn, nt))                          # Source wavelet
for isrc in range(srcn):
    wavelet[isrc,:] = source_wavelet(nt, dt, f0, 'ricker')

### receivers setup
temp = np.loadtxt(modelpath + 'receiver_coordinate.dat') # Receiver coordinates
recn  = temp.shape[0]                                    # Receiver number
recxz = np.zeros((srcn, recn, 2))                        # Receiver positions
for isrc in range(srcn):
    recxz[isrc,:,0] = temp[:,0]                          # Receiver x position (m)
    recxz[isrc,:,1] = temp[:,1]                          # Receiver z position (m)

### inversion parameter
misfit_type = 'Waveform'                           # 'Traveltime', 'Waveform', 'Globalcorrelation'
scheme      = 'LBFGS'                              # 'NLCG', 'LBFGS'
maxiter     = 50                                   # maximum iteration number,  i.e., 20
step_length = 0.02                                 # maximum update percentage, i.e., 0.05
vpmax       = 5000                                 # maximum allowed velocity
vpmin       = 1000                                 # minimum allowed velocity
marine_or_land = 'Marine'                          # 'Land' or 'Marine'

# gradient postprocess
grad_mute = 20                                     # mute source energy for 'Land', or water mask for 'Marine'
grad_smooth = 0                                    # gradient smooth radius 

# data filter
fre_filter = 'None'                                # 'Bandpass', 'Lowpass', 'Highpass', 'None'
fre_low  = 10                                      # low  frequency corner (units: Hz)
fre_high = 20                                      # high frequency corner (units: Hz)

# mute later arrival
mute_late_arrival = False                          # pick first break and mute later arrival
mute_late_window = 1.0                             # mute time window (units: time)

# mute offset
mute_offset_short  = False                         # mute short-offset traces 
mute_offset_long   = False                         # mute long-offset traces 
mute_offset_short_dis = 500                        # mute short-offset distance (units: m)
mute_offset_long_dis  = 10000                      # mute long-offset distance (units: m)

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

### posterior source inversion
source_inversion(simu, inv_offset=20000)
