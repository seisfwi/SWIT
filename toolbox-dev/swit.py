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


# import modules
import datetime, time
import sys
import numpy as np
from utils import read_yaml, load_model, source_wavelet
from acquisition import Config, Model, Source, Receiver
from solver import Solver

FWI = None
RTM = None


class SWIT(object):

    def __init__(self, config_file):
        
        print('*****************************************************')
        print('\n    Seismic Waveform Inversion Toolbox           \n')
        print('*****************************************************\n')

        # read config file
        try:
            par = read_yaml(config_file)
        except:
            raise FileNotFoundError('Config file not found: {}'.format(config_file))

        # job type
        self.job_type = par['Config']['job_type'].lower()

        # retrieve parameters
        try:
            # config parameters
            work_path = par['Config']['work_path']
            mpi_num = par['Config']['mpi_num']

            # model parameters
            nx = par['Model']['nx']
            nz = par['Model']['nz']
            dx = par['Model']['dx']
            dt = par['Model']['dt']
            nt = par['Model']['nt']
            pml = par['Model']['pml']
            vp = load_model(par['Model']['vp_path'], nx, nz)    # load model from file
            rho = load_model(par['Model']['rho_path'])          # load model from file 
            
            # source parameters
            f0 = par['Source']['f0']
            amp0 = par['Source']['amp0']
            src_type = par['Source']['type']

            # load source coordinates
            src_coord_file = par['Source']['coord_file']
            try:
                src_coord = np.loadtxt(src_coord_file)
            except:
                msg = 'Source coordinate file not found: {}'.format(src_coord_file)
                raise FileNotFoundError(msg)

            # set up source wavelet
            src_num = src_coord.shape[0]            
            wavelet  = np.zeros((src_num, nt))                          # Source wavelet
            if src_type in ['file']:
                try:
                    wavelet = np.load(par['Source']['wavelet_path'])
                except:
                    msg = 'Source wavelet file not found: {}'.format(par['Source']['wavelet_path'])
                    raise FileNotFoundError(msg)
            else:
                for isrc in range(src_num):
                    wavelet[isrc,:] = source_wavelet(amp0, nt, dt, f0, src_type)

            # receiver parameters
            rec_type = par['Receiver']['type']

            # load receiver coordinates
            rec_coord_file = par['Receiver']['coord_file']
            try:
                rec_coord = np.load(rec_coord_file)
            except:
                msg = 'Receiver coordinate file not found: {}'.format(rec_coord_file)
                raise FileNotFoundError(msg)

        except:
            raise KeyError('Some parameters are missing or wrong in the config file.')


        # configuration, model, source, receiver are required for all jobs
        self.config = Config(work_path, mpi_num)
        self.model = Model(nx, nz, dx, dt, nt, pml, vp, rho)
        self.source = Source(src_coord, wavelet, f0)
        self.receiver = Receiver(rec_coord, rec_type)
        self.solver = Solver(self.config, self.model, self.source, self.receiver)

        # processor, optimizer are optionally required for FWI or RTM
        if self.job_type in ['fwi', 'rtm']:

            # processor parameters
            self.processor = None

            # optimizer parameters
            if self.job_type in ['fwi']:
                self.optimizer = None


    def run(self, simu_tag = 'obs', save_snap = False):
    
        # Forward simulation
        if self.job_type in ['forward']:
            self.solver.forward(simu_tag = simu_tag, save_snap = save_snap)
        
        # Full Waveform Inversion
        elif self.job_type in ['fwi']:
            fwi = FWI(self.solver, self.processor, self.optimizer)
            fwi.run()
        
        # Reverse Time Migration
        elif self.job_type in ['rtm']:
            rtm = RTM(self.solver, self.processor)
            rtm.run()
            
        # Unknown job type
        else:
            msg = 'Support job types: Forward, FWI, RTM. \n'
            err = 'Unknown job type: {}'.format(self.config['Config']['job_type'])
            raise ValueError(msg + '\n' + err)


if __name__ == '__main__':
    
    # check if all required parameters are given
    if len(sys.argv) != 2:
        print("Usage: SWIT.py config.yaml")
        sys.exit(1)

    # read config file
    config_file = sys.argv[1]

    # initilize SWIT
    swit = SWIT(config_file)
    
    # set timer
    start_time = time.time()

    # run SWIT
    swit.run()

    # print running time
    consumed_time = datetime.timedelta(seconds = time.time() - start_time)
    print('SWIT running time: {}'.format(consumed_time))