###############################################################################
#
# SWIT v1.1: Seismic Waveform Inversion Toolbox
#
#   A Python package for seismic waveform inversion
#   Developed by Haipeng Li at USTC, updated on 2022-12-21 at Stanford
#   haipengl@mail.ustc.edu.cn, haipeng@stanford.edu
#
#   Main module
#
###############################################################################

import sys
sys.path.append('/homes/sep/haipeng/develop/SWIT-1.0/toolbox-dev/')

# import modules
import sys
import numpy as np

try:
    from acquisition import Config, Model, Source, Receiver
    from solver import Solver
    from optimizer import Optimizer
    from preprocessor import Preprocessor
    from utils import source_wavelet
    from tools import load_yaml, load_model, Config_object
    from workflow import FWI, RTM
except:
    raise ImportError('SWIT modules not found.')


class SWIT(Config_object):
    '''
        SWIT class defines the main class for seismic waveform inversion

    parameters:
    ----------
        config_file: str
            path to the config file (yaml format)
    '''

    def __init__(self, config_file):
        
        # inherit Config_object containing the config parameters
        Config_object.__init__(self, load_yaml(config_file))

        # check the job type
        self.check_job()

        # initialize solver
        self.init_solver()

        if 'FWI' in self.job_type:

            # initialize optimizer
            self.init_optimizer()

            # initialize preprocessor
            self.init_preprocessor()

        elif 'RTM' in self.job_type:
            raise NotImplementedError('RTM is not implemented yet.')
 

    def check_job(self):
        ''' check the job type and print the message
        '''

        # check the job type and print the message
        job_str = ''
        self.job_type = [job.upper() for job in self.job_type]
        for job in self.job_type:
            if job not in ['FORWARD', 'FWI', 'RTM']:
                raise ValueError('Job type {} not supported.'.format(job))
            job_str += job + '  '

        print('*****************************************************')
        print('\n        Seismic Waveform Inversion Toolbox        \n')
        print('        Job Flow: {}                                  '.format(job_str))
        print('*****************************************************\n')


    def init_solver(self):
        ''' initialize the solver class
        '''

        # load model files from config file
        vp  = load_model(self.vp_file,  self.nx, self.nz)  # load vp  model from file
        rho = load_model(self.rho_file, self.nx, self.nz)  # load rho model from file
        
        # load source coordinates
        src_coord = np.load(self.src_coord_file)

        # load source wavelet from file or generate it from the config file
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
        config = Config(self.path, self.mpi_num, cpu_max_num = self.cpu_max_num, fig_aspect = self.fig_aspect)
        model = Model(self.nx, self.nz, self.dx, self.dt, self.nt, self.pml, vp, rho, acquisition_type = self.acquisition_type)
        source = Source(src_coord, wavelet, self.f0)
        receiver = Receiver(rec_coord, self.rec_comp)
        
        # initialize solver
        self.solver = Solver(config, model, source, receiver)


    def init_optimizer(self):
        ''' initialize the optimizer class
        '''

        # load initial model files
        vp_init = load_model(self.vp_init_file, self.nx, self.nz)            # load initial vp model from file
        rho_init = load_model(self.rho_init_file, self.nx, self.nz)          # load initial rho model from file
     
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
    
        # Forward simulation
        if 'FORWARD' in self.job_type:
            self.solver.run()
        
        # Full Waveform Inversion
        if 'FWI' in self.job_type:
            fwi = FWI(self.solver, self.optimizer, self.preprocessor)
            fwi.run()
        
        # Reverse Time Migration
        if 'RTM' in self.job_type:
            rtm = RTM(self.solver, self.preprocessor)
            rtm.run()


if __name__ == '__main__':
    
    # check if all required parameters are given
    if len(sys.argv) != 2:
        print("Usage: SWIT.py config.yaml")
        sys.exit(1)

    # parse config file
    config_file = sys.argv[1]

    # initilize SWIT with config
    swit = SWIT(config_file)
    
    # run SWIT workflow
    swit.run()
