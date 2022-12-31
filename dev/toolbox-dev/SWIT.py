###############################################################################
#
# SWIT v1.1: Seismic Waveform Inversion Toolbox
#
#   A Python package for seismic waveform inversion
#   By Haipeng Li at USTC & Stanford
#   Email: haipengl@mail.ustc.edu.cn, haipeng@stanford.edu 
#
#   SWIT class defines the job flow control
#
###############################################################################

import sys
sys.path.append('/homes/sep/haipeng/develop/SWIT-1.0/toolbox-dev/')

import numpy as np

try:
    from survey import System, Model, Receiver, Source
    from optimizer import Optimizer
    from preprocessor import Preprocessor
    from solver import Solver
    from tools import load_model, load_yaml
    from utils import source_wavelet
    from workflow import FWI, RTM, Configuration
except:
    raise ImportError('SWIT modules not found.')


class SWIT(Configuration):
    '''
        SWIT class defines the job flow control

    parameters:
    ----------
        config_file: str
            path to the config file (yaml format)
    '''

    def __init__(self, config_file):
        
        # parse the config file by inheriting Configuration class
        Configuration.__init__(self, load_yaml(config_file))

        # check the job type
        self.check_config()

    
    def check_config(self):
        ''' check the config file
        '''
        # parameter list
        par_list  = ['path', 'job_type', 'max_cpu_num', 'mpi_num', 'fig_aspect']
        par_list += ['dt', 'dx', 'nt', 'nx', 'nz', 'pml', 'vp_file', 'rho_file']
        par_list += ['rec_comp', 'rec_coord_file']
        par_list += ['amp0', 'f0', 'src_type', 'src_coord_file', 'wavelet_file']
        par_list += ['vp_init_file', 'rho_init_file', 'misfit_type', 'method', 
                     'niter_max', 'bound', 'vp_min', 'vp_max', 'grad_smooth_size', 
                     'update_vpmax', 'grad_mask', 'debug']
        par_list += ['filter_data', 'filter_high', 'filter_low', 'mute_near_offset', 
                     'mute_near_distance', 'mute_far_offset']
        
        # check the parameters
        for par in par_list:
            if par not in self.__dict__:
                raise ValueError('Parameter [{}] is not found in config file.'.format(par))
        
        # check the job type
        job_str = ''
        self.job_type = [job.upper() for job in self.job_type]
        for job in self.job_type:
            if job not in ['FORWARD', 'FWI', 'RTM']:
                raise ValueError('Job type {} not supported.'.format(job))
            job_str += job + ' '

        # print the job workflow
        print('****************************************************************')
        print('          SEISMIC WAVEFORM INVERSION TOOLBOX (V-1.1)          \n')
        print('  Job workflow: {}                                              '.format(job_str))
        print('  Job workpath: {}                                              '.format(self.path))
        print('****************************************************************\n')


    def init_solver(self):
        ''' initialize the solver class
        '''

        # load model files
        vp  = load_model(self.vp_file,  self.nx, self.nz)  # load vp  model from file
        rho = load_model(self.rho_file, self.nx, self.nz)  # load rho model from file
        
        # load source coordinates
        src_coord = np.load(self.src_coord_file)

        # load source wavelet from file or generate it from the system file
        if self.src_type.lower() in ['file']:
            wavelet = np.load(self.wavelet_path)
        else:
            src_num = src_coord.shape[0]            
            wavelet  = np.zeros((src_num, self.nt))                          # Source wavelet
            for isrc in range(src_num):
                wavelet[isrc,:] = source_wavelet(self.amp0, self.nt, self.dt, self.f0, self.src_type)

        # load receiver coordinates in the format of npz file
        rec_file = np.load(self.rec_coord_file)
        rec_coord = []
        for key in rec_file.keys():
            rec_coord.append(rec_file[key])

        # configuration, model, source, receiver are required for all jobs
        system = System(self.path, self.mpi_num, max_cpu_num = self.max_cpu_num, fig_aspect = self.fig_aspect)
        model = Model(self.nx, self.nz, self.dx, self.dt, self.nt, self.pml, vp, rho)
        source = Source(src_coord, wavelet, self.f0)
        receiver = Receiver(rec_coord, self.rec_comp)
        
        # initialize solver
        self.solver = Solver(system, model, source, receiver)


    def init_optimizer(self):
        ''' initialize the optimizer class
        '''

        # load initial model files for FWI (vp and rho)
        vp_init = load_model(self.vp_init_file, self.nx, self.nz)
        rho_init = load_model(self.rho_init_file, self.nx, self.nz)
     
        # load gradient mask if a path to mask is specified, otherwise set to None
        if self.grad_mask is None or len(self.grad_mask) == 0:
            grad_mask = None
        else:
            grad_mask = load_model(self.grad_mask, self.nx, self.nz)

        # initialize optimizer
        self.optimizer = Optimizer(vp_init = vp_init, rho_init = rho_init, 
                misfit_type = self.misfit_type, method = self.method, niter_max = self.niter_max,
                bound = self.bound, vp_max = self.vp_max, vp_min = self.vp_min, update_vpmax = self.update_vpmax,
                grad_smooth_size = self.grad_smooth_size, grad_mask = grad_mask, debug = self.debug)


    def init_preprocessor(self):
        ''' initialize the processor class
        '''

        ## optimizer
        self.preprocessor = Preprocessor(filter_data = self.filter_data, 
                                         filter_low = self.filter_low, 
                                         filter_high = self.filter_high, 
                                         mute_late_arrival = self.mute_late_arrival, 
                                         mute_late_size = self.mute_late_size, 
                                         normalize_data = self.normalize_data,
                                         mute_near_offset = self.mute_near_offset, 
                                         mute_near_distance = self.mute_near_distance, 
                                         mute_far_offset = self.mute_far_offset, 
                                         mute_far_distance = self.mute_far_distance)


    def run(self):
        ''' run the job workflow
        '''

        # solver is required for all jobs
        self.init_solver()

        # Forward Modeling
        if 'FORWARD' in self.job_type:
            self.solver.run()

        # Full Waveform Inversion
        if 'FWI' in self.job_type:
            self.init_optimizer()
            self.init_preprocessor()
            fwi = FWI(self.solver, self.optimizer, self.preprocessor)
            fwi.run()
        
        # Reverse Time Migration
        if 'RTM' in self.job_type:
            self.init_preprocessor()

            # load velocity and density model for RTM
            vp_rtm  = load_model(self.vp_init_file, self.nx, self.nz)
            rho_rtm = load_model(self.rho_init_file, self.nx, self.nz)
    
            rtm = RTM(self.solver, self.preprocessor)
            rtm.run(vp = vp_rtm, rho = rho_rtm)


if __name__ == '__main__':
    
    # check if all required parameters are given
    if len(sys.argv) != 2:
        print("Usage: SWIT.py config.yaml")
        sys.exit(1)

    # initilize SWIT with config file
    swit = SWIT(sys.argv[1])
    
    # run SWIT workflow
    swit.run()
