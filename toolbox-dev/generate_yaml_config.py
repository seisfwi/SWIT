
import yaml


def read_yaml(filename):
    ''' A function to read YAML file
    '''
    with open(filename) as f:
        config = yaml.safe_load(f)
 
    return config


def write_yaml(config, filename):
    ''' A function to write YAML file
    '''
    with open(filename, 'w') as f:
        yaml.dump(config, f)


# create an empty dictionary
config = {}

# config parameters
config['Config'] = {}
config['Config']['job_type'] = 'forward'                               # 'forward', 'FWI', 'RTM'
config['Config']['path'] = '/scr2/haipeng/SWIT-1.1/01_fwi/'            # Working path
config['Config']['mpi_num'] = 40                                       # Number of MPI processes
config['Config']['cpu_max_num'] = 48                                   # Maximum number of CPUs to use
config['Config']['fig_aspect'] = 1.0                                   # Aspect ratio of the figure

# model parameters
config['Model'] = {}
config['Model']['nx'] = 481
config['Model']['nz'] = 141
config['Model']['pml'] = 40
config['Model']['dx'] = 25
config['Model']['dt'] = 0.002
config['Model']['nt'] = 2001
config['Model']['acquisition_type'] = 'land'                       # 'land', 'marine'
config['Model']['vp_true'] = '/homes/sep/haipeng/develop/SWIT-1.0/examples-dev/case-01-Marmousi2/acquisition/vp_true.npy'
config['Model']['rho_true'] = '/homes/sep/haipeng/develop/SWIT-1.0/examples-dev/case-01-Marmousi2/acquisition/rho_true.npy'

# source parameters
config['Source'] = {}
config['Source']['f0'] = 5.0
config['Source']['amp0'] = 1.0
config['Source']['type'] = 'ricker'   # 'ricker', 'gaussian', 'ramp', 'file'
config['Source']['coord_file'] = '/homes/sep/haipeng/develop/SWIT-1.0/examples-dev/case-01-Marmousi2/acquisition/src_coord.npy'
config['Source']['wavelet_file'] = '/homes/sep/haipeng/develop/SWIT-1.0/examples-dev/case-01-Marmousi2/acquisition/wavelets.npy'

# receiver parameters
config['Receiver'] = {}
config['Receiver']['coord_file'] = '/homes/sep/haipeng/develop/SWIT-1.0/examples-dev/case-01-Marmousi2/acquisition/rec_coord.npz'
config['Receiver']['comp'] = 'p'   # 'vx', 'vz', 'p'

# preprocessor parameters
config['Preprocessor'] = {}
config['Preprocessor']['filter_data'] = 'none'
config['Preprocessor']['filter_low'] = 5.0
config['Preprocessor']['filter_high'] = 10.0
config['Preprocessor']['mute_late_arrival'] = True
config['Preprocessor']['mute_late_size'] = 0.5
config['Preprocessor']['normalize_data'] = False
config['Preprocessor']['mute_near_offset'] = True
config['Preprocessor']['mute_near_distance'] = 500
config['Preprocessor']['mute_far_offset'] = True
config['Preprocessor']['mute_far_distance'] = 8000

# optimizer parameters
config['Optimizer'] = {}
config['Optimizer']['vp_init']  = '/homes/sep/haipeng/develop/SWIT-1.0/examples-dev/case-01-Marmousi2/acquisition/vp_init.npy'
config['Optimizer']['rho_init'] = '/homes/sep/haipeng/develop/SWIT-1.0/examples-dev/case-01-Marmousi2/acquisition/rho_init.npy'
config['Optimizer']['misfit_type'] = 'waveform'   # 'waveform', 'amplitude'
config['Optimizer']['method'] = 'LBFGS'   # 'LBFGS', 'BFGS', 'CG', 'Newton-CG', 'Powell', 'TNC', 'COBYLA', 'SLSQP', 'trust-constr', 'dogleg', 'trust-ncg', 'trust-exact', 'trust-krylov'
config['Optimizer']['niter_max'] = 20
config['Optimizer']['bound'] = False
config['Optimizer']['vp_max'] = 5000
config['Optimizer']['vp_min'] = 1000
config['Optimizer']['update_vpmax'] = 50
config['Optimizer']['grad_smooth_size'] = 0
config['Optimizer']['grad_mask'] = None
config['Optimizer']['debug'] = False


write_yaml(config, 'config.yaml')


# work path

# with open('config.yaml') as file:
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


# # with open(r'config.yaml') as file:
# #     documents = yaml.full_load(file)

# #     for item, doc in documents.items():
# #         print(item, ":", doc)


