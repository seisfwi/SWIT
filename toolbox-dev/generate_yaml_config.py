
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
config['Config']['job_type'] = 'forward'   # 'forward', 'FWI', 'RTM'
config['Config']['work_path'] = '/scr2/haipeng/SWIT-1.1/00_forward/'   # Working path
config['Config']['mpi_num'] = 40   # Number of MPI processes

# model parameters
config['Model'] = {}
config['Model']['nx'] = 481
config['Model']['nz'] = 141
config['Model']['pml'] = 40
config['Model']['dx'] = 25
config['Model']['dt'] = 0.002
config['Model']['nt'] = 2001
config['Model']['vp_file'] = '/homes/sep/haipeng/develop/SWIT-1.0/examples/case-01-marmousi2/model/Marmousi_481_141_25m.dat'
config['Model']['rho_file'] = '/homes/sep/haipeng/develop/SWIT-1.0/examples/case-01-marmousi2/model/Marmousi_481_141_25m.dat'

# source parameters
config['Source'] = {}
config['Source']['f0'] = 5.0
config['Source']['amp0'] = 1.0
config['Source']['type'] = 'ricker'   # 'ricker', 'gaussian', 'ramp', 'file'
config['Source']['coord_file'] = '/homes/sep/haipeng/develop/SWIT-1.0/examples/case-01-marmousi2/source/source_coordinate.npy'

# receiver parameters
config['Receiver'] = {}
config['Receiver']['recxz'] = '/homes/sep/haipeng/develop/SWIT-1.0/examples/case-01-marmousi2/source/receiver_coordinate.npy'
config['Receiver']['rec_comp'] = 'vz'   # 'vx', 'vz', 'p'

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


