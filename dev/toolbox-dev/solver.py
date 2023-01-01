###############################################################################
#
# SWIT v1.1: Seismic Waveform Inversion Toolbox
#
#   A Python package for seismic waveform inversion
#   By Haipeng Li at USTC & Stanford
#   Email: haipengl@mail.ustc.edu.cn, haipeng@stanford.edu 
#
#   Solver class for calling fd2mpi wavefield solver (2D acoustic wave equation)
#
###############################################################################

import os
import subprocess

import numpy as np
from scipy import integrate
from survey import Survey
from tools import save_float, load_float


class Solver(Survey):
    ''' Solver class for calling fd2mpi wavefield solver (2D acoustic wave equation)
        One can replace this class with their own solvers

    parameters:
    ----------
        system: System class
            system configuration
        model: Model class
            model parameters
        source: Source class
            source parameters
        receiver: Receiver class
            receiver parameters
    '''

    def __init__(self, system, model, source, receiver):
        ''' Initialize the solver class
        '''

        # inherit Survey class
        Survey.__init__(self, system, model, source, receiver)

        # print solver information
        self.__info__()


    def run(self, simu_type = 'forward', simu_tag = 'obs', data_format = 'bin', save_snap = False, save_boundary = False):
        ''' Forward solver

        Parameters
        ----------  
            simu_type : str
                simulation type, 'forward', 'adjoint' or 'gradient'
            simu_tag : str
                simulation tag that makes the simulation unique, 'obs', 'syn' or 'adj'
            data_format : str
                data format, 'bin' or 'su'
            save_snap : boolean
                save snapshot data or not
            save_boundary : boolean
                save boundary data or not, for wavefield reconstruction

        Returns (only for 'gradient' simulation)
        -------
            gradient : 2D array
                gradient of misfit function
            for_illum : 2D array
                forward illumination
            adj_illum : 2D array
                adjoint illumination
        '''

        # check the specification of the simulation type
        simu_tapes = ['forward', 'adjoint', 'gradient']
        if simu_type not in simu_tapes:
            msg = 'Solver: simu_type must be one of {}'.format(simu_tapes)
            err = 'Solver: unknown simulation type: {}'.format(simu_type)
            raise ValueError(msg + '\n' + err)

        # check the specification of the data format
        formats = ['bin', 'su']
        if data_format not in formats:
            msg = 'Solver: data_format must be one of {}'.format(formats)
            err = 'Solver: unknown data format: {}'.format(data_format)
            raise ValueError(msg + '\n' + err)

        # check mpi and wavefield solver fd2dmpi
        for cmd in ['mpirun', 'fd2dmpi']:
            if subprocess.getstatusoutput('which ' + cmd)[0] != 0:
                raise ValueError('Solver: cannot find {}, check its installation.'.format(cmd))

        # create working directory
        for isrc in range(self.source.num):
            folder = self.system.path + 'data/{}/src{}'.format(simu_tag, isrc+1)
            if not os.path.exists(folder):
                os.system('mkdir -p {}'.format(folder))

        # prepare the system file for the solver
        self.write_config(simu_type, simu_tag, data_format, save_snap, save_boundary)

        # run the solver with mpi
        cmd = 'mpirun -np {}  fd2dmpi config={}'.format(self.system.mpi_num, self.system.path + 'config/solver.config')
        status = subprocess.getstatusoutput(cmd)

        # check the final status of the solver
        if status[0]:
            print(status[1])
            raise ValueError('Solver: failed, check the solver configuration.')

        # read gradient and wavefield illumination for gradient simulation
        if simu_type == 'gradient':
            nx = self.model.nx
            nz = self.model.nz
            gradient  = np.zeros((nx, nz))
            for_illum = np.zeros((nx, nz))
            adj_illum = np.zeros((nx, nz))

            # read gradient from binary files and sum over all sources
            for isrc in range(self.source.num):
                data_path  = os.path.join(self.system.path, 'data/syn/src{}/'.format(isrc+1))
                gradient  += load_float(data_path + 'vp_gradient.bin').reshape(nx, nz)
                for_illum += load_float(data_path + 'forward_illumination.bin').reshape(nx, nz)
                adj_illum += load_float(data_path + 'adjoint_illumination.bin').reshape(nx, nz)

            # return
            return gradient, for_illum, adj_illum


    def set_model(self, vp = None, rho = None):
        ''' Reset model parameters

        Parameters
        ----------
            vp : 1D vector or 2D array
                P-wave velocity for reset model
            rho : 1D vector or 2D array
                density for reset model
        '''

        if vp is not None:
            try:
                self.model.vp = vp.reshape(self.model.vp.shape)
            except:
               raise ValueError('Solver: the shape of vp does not match the shape of the model')

        if rho is not None:
            try:
                self.model.rho = rho.reshape(self.model.rho.shape)
            except:
               raise ValueError('Solver: the shape of rho does not match the shape of the model')
        

    def write_config(self, simu_type, simu_tag, data_format, save_snap, save_boundary):
        ''' Write configuration files for forward solver

                Note: the adjoint source wavelet should be saved as: config/wavelet/src1_adj.bin, ... 
                before running adjoint or gradient job
        '''

        # save source wavelet files (src1.bin, ...)
        # the source wavelet is integrated to count for the derivative effect of 1st order FD
        for isrc in range(self.source.num):
            src = integrate.cumtrapz(self.source.wavelet[isrc,:], axis=-1, initial=0)
            src_path = os.path.join(self.system.path, 'config/wavelet/src{}.bin'.format(isrc+1))
            save_float(src_path, src)

        # save P-wave velocity and density files (vp.bin, rho.bin)
        save_float(os.path.join(self.system.path, 'config/vp.bin'),  self.model.vp)
        save_float(os.path.join(self.system.path, 'config/rho.bin'), self.model.rho)

        # save geometry file (geometry.config)
        fp = open(os.path.join(self.system.path + 'config/geometry.config'), "w")
        for isrc in range(self.source.num):
            for irec in range(len(self.receiver.coord[isrc])):
                # column:       1        2      3     4     5     6        7
                # parameter: S_index  R_index  Sx    Sz    Rx    Rz  S_is_alive(0/1))
                fp.write('%6i %6i %10.1f %10.1f %10.1f %10.1f %10.1f\n'%(
                    isrc+1, irec+1, self.source.coord[isrc, 0], self.source.coord[isrc, 1], 
                    self.receiver.coord[isrc][irec, 0], self.receiver.coord[isrc][irec, 1], 1))
        fp.close()

        # save solver system file (solver.config)
        parpath = os.path.join(self.system.path, 'config/solver.config')
        fp = open(parpath, "w")
        fp.write('######################################### \n')
        fp.write('#                                         \n')
        fp.write('#     fd2dmpi input parameter file        \n')
        fp.write('#                                         \n')
        fp.write('######################################### \n')
        fp.write('                                          \n')
        fp.write('jobtype=%s\n'         % simu_type)
        fp.write('COORD_FILE=%s\n'      % (self.system.path + 'config/geometry.config'))
        fp.write('SOURCE_FILE=%s\n'     % (self.system.path + 'config/wavelet/src'))
        fp.write('DATA_OUT=%s\n'        % (self.system.path + 'data/{}/src'.format(simu_tag)))
        fp.write('DATA_COMP=%s\n'       % (self.receiver.comp.lower()))
        fp.write('VEL_IN=%s\n'          % (self.system.path + 'config/vp.bin'))
        fp.write('DENSITYFILE=%s\n'     % (self.system.path + 'config/rho.bin'))
        fp.write('FILEFORMAT=%s\n'      % data_format)
        fp.write('NX=%d\n'              % self.model.nx)
        fp.write('NZ=%d\n'              % self.model.nz)
        fp.write('DX=%f\n'              % self.model.dx)
        fp.write('NPML=%d\n'            % self.model.pml)
        fp.write('NT_WORK=%d\n'         % self.model.nt)
        fp.write('DT_WORK=%f\n'         % self.model.dt)
        fp.write('FREESURFACE=1\n')                              # free surface is always on
        fp.write('STORE_SNAP=%d\n'      % save_snap)             # save snapshot or not
        fp.write('STORE_STEP=%d\n'      % self.model.save_step)  # save snapshot every 10 time steps (default)
        fp.write('STORE_BOUNDAARY=%d\n' % save_boundary)         # save boundary wavefield for reconstruction
        fp.close()


    def __info__(self):
        ''' Print solver information
        '''

        print('\nSolver information:')
        print('    Model dimension   : nx = {}, nz = {}, dx = {} m, x = 0 ~ {:.2f} km,  z = 0 ~ {:.2f} km'.format(self.model.nx, self.model.nz, self.model.dx, self.model.x[-1]/1000, self.model.z[-1]/1000))
        print('    Data  dimension   : nt = {}, dt = {:.2f} ms, t = 0 ~ {:.2f} s'.format(self.model.nt, self.model.dt*1000, self.model.t[-1]))
        print('    Data acquisition  : {} sources, {}-component receivers'.format(self.source.num, self.receiver.comp))
        print('    P-wave velocity   : {:.2f} ~ {:.2f} m/s'.format(self.model.vp.min(), self.model.vp.max()))
        print('    MPI information   : {} tasks run in parallel'.format(self.system.mpi_num))
