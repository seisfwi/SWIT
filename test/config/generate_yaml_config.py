
import yaml


def read_yaml(filename):
    ''' A function to read YAML file
    '''
    with open(filename) as f:
        system = yaml.safe_load(f)
 
    return system


def write_yaml(system, filename):
    ''' A function to write YAML file
    '''
    with open(filename, 'w') as f:
        yaml.dump(system, f)


# create an empty dictionary
system = {}

# system parameters
system['Config'] = {}
system['Config']['job_type'] = 'forward'                               # 'forward', 'FWI', 'RTM'
system['Config']['path'] = '/scr2/haipeng/SWIT-1.1/01_fwi/'            # Working path
system['Config']['mpi_num'] = 40                                       # Number of MPI processes
system['Config']['max_cpu_num'] = 48                                   # Maximum number of CPUs to use
system['Config']['fig_aspect'] = 1.0                                   # Aspect ratio of the figure

# model parameters
system['Model'] = {}
system['Model']['nx'] = 481
system['Model']['nz'] = 141
system['Model']['pml'] = 40
system['Model']['dx'] = 25
system['Model']['dt'] = 0.002
system['Model']['nt'] = 2001
system['Model']['acquisition_type'] = 'land'                       # 'land', 'marine'
system['Model']['vp_true'] = '/homes/sep/haipeng/develop/SWIT-1.0/examples-dev/case-01-Marmousi2/acquisition/vp_true.npy'
system['Model']['rho_true'] = '/homes/sep/haipeng/develop/SWIT-1.0/examples-dev/case-01-Marmousi2/acquisition/rho_true.npy'

# source parameters
system['Source'] = {}
system['Source']['f0'] = 5.0
system['Source']['amp0'] = 1.0
system['Source']['type'] = 'ricker'   # 'ricker', 'gaussian', 'ramp', 'file'
system['Source']['coord_file'] = '/homes/sep/haipeng/develop/SWIT-1.0/examples-dev/case-01-Marmousi2/acquisition/src_coord.npy'
system['Source']['wavelet_file'] = '/homes/sep/haipeng/develop/SWIT-1.0/examples-dev/case-01-Marmousi2/acquisition/wavelets.npy'

# receiver parameters
system['Receiver'] = {}
system['Receiver']['coord_file'] = '/homes/sep/haipeng/develop/SWIT-1.0/examples-dev/case-01-Marmousi2/acquisition/rec_coord.npz'
system['Receiver']['comp'] = 'p'   # 'vx', 'vz', 'p'

# preprocessor parameters
system['Preprocessor'] = {}
system['Preprocessor']['filter_data'] = 'none'
system['Preprocessor']['filter_low'] = 5.0
system['Preprocessor']['filter_high'] = 10.0
system['Preprocessor']['mute_late_arrival'] = True
system['Preprocessor']['mute_late_size'] = 0.5
system['Preprocessor']['normalize_data'] = False
system['Preprocessor']['mute_near_offset'] = True
system['Preprocessor']['mute_near_distance'] = 500
system['Preprocessor']['mute_far_offset'] = True
system['Preprocessor']['mute_far_distance'] = 8000

# optimizer parameters
system['Optimizer'] = {}
system['Optimizer']['vp_init']  = '/homes/sep/haipeng/develop/SWIT-1.0/examples-dev/case-01-Marmousi2/acquisition/vp_init.npy'
system['Optimizer']['rho_init'] = '/homes/sep/haipeng/develop/SWIT-1.0/examples-dev/case-01-Marmousi2/acquisition/rho_init.npy'
system['Optimizer']['misfit_type'] = 'waveform'   # 'waveform', 'amplitude'
system['Optimizer']['method'] = 'LBFGS'   # 'LBFGS', 'BFGS', 'CG', 'Newton-CG', 'Powell', 'TNC', 'COBYLA', 'SLSQP', 'trust-constr', 'dogleg', 'trust-ncg', 'trust-exact', 'trust-krylov'
system['Optimizer']['niter_max'] = 20
system['Optimizer']['bound'] = False
system['Optimizer']['vp_max'] = 5000
system['Optimizer']['vp_min'] = 1000
system['Optimizer']['update_vpmax'] = 50
system['Optimizer']['grad_smooth_size'] = 0
system['Optimizer']['grad_mask'] = None
system['Optimizer']['debug'] = False


write_yaml(system, 'system.yaml')


# work path

# with open('system.yaml') as file:
#     try:
#         databaseConfig = yaml.safe_load(file)   
#         print(databaseConfig)
#     except yaml.YAMLError as exc:
#         print(exc)


# # read_categories.py file

# import yaml

# work_path = '/scr2/haipeng/SWIT-1.1/00_forward/'   # Working path
# mpi_num   = 49                                     # MPI process for fd2dmpi

# ### model setup
# nx,  nz, pml = [481,  141, 40]                     # Grid number along x and z directions, Grid number for PML layers (use a large one)
# dx,  dt,  nt = [25, 0.002, 2001]                   # Grid size, time interval, and time step

# # true model
# modelpath = '/homes/sep/haipeng/develop/SWIT-1.0/examples/case-01-marmousi2/model/'
# vp_true = modelpath + 'Marmousi_481_141_25m.dat'
# rho_true = modelpath + 'Marmousi_481_141_25m_rho.dat'

# ### sources setup 
# f0    = 5.0                                              # Dominant frequency in Hz
# amp0  = 1.0                                              # Amplitude of the source wavelet
# srcxz = modelpath + 'source_coordinate.dat'  # Source coordinates

# ### receivers setup
# recxz = modelpath + 'receiver_coordinate.dat' # Receiver coordinates


# # with open(r'system.yaml') as file:
# #     documents = yaml.full_load(file)

# #     for item, doc in documents.items():
# #         print(item, ":", doc)


