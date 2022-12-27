###############################################################################
#
# SWIT v1.1: Seismic Waveform Inversion Toolbox
#
#   A Python package for seismic waveform inversion
#   Developed by Haipeng Li at USTC, updated on 2022-12-21 at Stanford
#   haipengl@mail.ustc.edu.cn, haipeng@stanford.edu
#
#   Solver module
#
###############################################################################


import os
import subprocess
import numpy as np
from scipy import integrate

from tools import save_float

class Solver(object):
    ''' Acoustic Wavefiled Solver class
    '''

    def __init__(self, config, model, source, receiver):
        ''' Initialize solver
        '''
        # basic class variables for solver
        self.config = config
        self.model = model
        self.source = source
        self.receiver = receiver

        # check the solver
        self.__check__()

        # print solver information
        self.__info__()


    def forward(self, simu_tag = 'obs', save_snap = 0):
        ''' Forward solver
        '''

        # create working directory
        for isrc in range(self.source.num):
            ifolder = self.config.path + 'data/{}/src{}'.format(simu_tag, isrc+1)
            if not os.path.exists(ifolder) and save_snap == 1:
                os.system('mkdir -p %s' % ifolder)

        # # prepare the forward source
        # src = integrate.cumtrapz(self.source.wavelet, axis=-1, initial=0)

        # prepare the config file for forward solver
        self.prepare_configfile('forward', simu_tag, save_snap)

        # check mpi and wavefield solver
        if os.system('which mpirun') != 0:
            raise ValueError('Cannot find mpirun, please check your mpi installation.')
        if os.system('which fd2dmpi') != 0:
            raise ValueError('Cannot find fd2dmpi, please check your wavefield solver installation.')

        # run the forward solver with mpi
        cmd = 'mpirun -np {}  fd2dmpi config={}'.format(self.config.mpi_num, self.config.path + 'config/solver.config')
        status = subprocess.getstatusoutput(cmd)

        # check the status of forward solver
        if status[0]:
            print(status[1])
            raise ValueError('Forward solver crash')


    def prepare_configfile(self, simu_type, simu_tag, save_snap):
        ''' Prepare configuration files for forward solver, including:
                1. source wavelet files: src1.bin, src2.bin, ...
                2. velocity and density files: vp.bin, rho.bin
                3. geometry config file: geometry.config
                4. solver config file: solver.config
        '''

        # create configfile directory
        if not os.path.exists(self.config.path + 'config'):
            os.system('mkdir -p %s' % (self.config.path + 'config'))
            os.system('mkdir -p %s' % (self.config.path + 'config/wavelet'))

        # save source wavelet files: src1.bin, src2.bin, ...
        for isrc in range(self.source.num):
            srcpath = os.path.join(self.config.path, 'config/wavelet/src{}.bin'.format(isrc+1))
            save_float(srcpath, self.source.wavelet[isrc, :])
                
        # save P-wave velocity and density files: vp.bin, rho.bin
        save_float(os.path.join(self.config.path, 'config/vp.bin'), self.model.vp)
        save_float(os.path.join(self.config.path, 'config/rho.bin'), self.model.rho)

        # write geometry config file: geometry.config
        fp = open(os.path.join(self.config.path + 'config/geometry.config'), "w")
        for isrc in range(self.source.num):
            for irec in range(len(self.receiver.coord[isrc])):
                # column:       1        2      3     4     5     6        7
                # parameter: S_index  R_index  Sx    Sz    Rx    Rz  S_is_alive(0/1))
                fp.write('%6i %6i %10.1f %10.1f %10.1f %10.1f %10.1f\n'%(
                    isrc+1, irec+1, self.source.coord[isrc, 0], self.source.coord[isrc, 1], 
                    self.receiver.coord[isrc][irec, 0], self.receiver.coord[isrc][irec, 1], 1))
        fp.close()

        # write solver config file: solver.config
        parpath = os.path.join(self.config.path, 'config/solver.config')
        fp = open(parpath, "w")
        fp.write('######################################### \n')
        fp.write('#                                         \n')
        fp.write('#     fd2dmpi input parameter file        \n')
        fp.write('#                                         \n')
        fp.write('######################################### \n')
        fp.write('                                          \n')

        # specify the simulation type
        if simu_type in ['forward']:
            fp.write('jobtype=forward\n')
        elif simu_type in ['adjoint']:
            fp.write('jobtype=adjoint\n')
        elif simu_type in ['gradient']:
            fp.write('jobtype=gradient\n')
        else:
            msg = "simu_type must be 'forward', 'adjoint' or 'gradient' \n"
            err = 'Unknown simulation type: {}'.format(simu_type)
            raise ValueError(msg + err)

        fp.write('COORD_FILE=%s\n'    % (self.config.path + 'config/geometry.config'))
        fp.write('SOURCE_FILE=%s\n'   % (self.config.path + 'config/wavelet/src'))
        fp.write('DATA_OUT=%s\n'      % (self.config.path + 'data/{}/src'.format(simu_tag)))
        fp.write('VEL_IN=%s\n'        % (self.config.path + 'config/vp.bin'))
        fp.write('DENSITYFILE=%s\n'   % (self.config.path + 'config/rho.bin'))
        fp.write('FILEFORMAT=bin\n')
        fp.write('NX=%d\n'            % self.model.nx)
        fp.write('NZ=%d\n'            % self.model.nz)
        fp.write('DX=%f\n'            % self.model.dx)
        fp.write('NPML=%d\n'          % self.model.pml)
        fp.write('NT_WORK=%d\n'       % self.model.nt)
        fp.write('DT_WORK=%f\n'       % self.model.dt)
        fp.write('FREESURFACE=1\n')             # Free surface is always on
        fp.write('STORE_SNAP=%d\n'    % save_snap)
        fp.write('STORE_STEP=%d\n'    % (10)) # save snapshot every 10 time steps
        fp.close()



    # def adjoint(self, savesnap = 0):
    #     ''' Adjoint solver
    #     '''
    
    #     simu_type = 'adj'
        
    #     # change path and clean previous data (always, even empty)
    #     cleandata(self.config.homepath + 'data/adj/')

    #     # create working directory
    #     for isrc in range(self.source.srcn):
    #         ifolder = self.config.homepath + 'data/%s/src%d_snapshot'%(simu_type, isrc+1)
    #         if not os.path.exists(ifolder) and savesnap == 1:
    #             os.system('mkdir %s' % ifolder)

    #     # prepare the forward source
    #     src = integrate.cumtrapz(self.source.wavelet, axis=-1, initial=0)

    #     # write parameters and model files
    #     self.write_parfile(simu_type, src, savesnap=savesnap)

    #     # check mpi and wavefield solver
    #     if os.system('which mpirun') != 0:
    #         raise ValueError('Cannot find mpirun, please check your mpi installation.')
    #     if os.system('which fd2dmpi') != 0:
    #         raise ValueError('Cannot find fd2dmpi, please check your wavefield solver installation.')

    #     # submit job
    #     os.chdir(self.config.homepath)
    #     solver_cmd = 'mpirun -np %d  fd2dmpi par=%s' % (self.config.mpiproc, 
    #                 self.config.homepath + 'parfile/forward_parfile/parfile')
    #     status = subprocess.getstatusoutput(solver_cmd)
    #     if status[0]:
    #         print(status[1])
    #         raise ValueError('Forward solver crash')

  
    # def dot_product(self, simu_type = 'obs'):
    #     ''' Calculate dot product of two wavefields
    #     '''
    #     pass


    def __check__(self):
        ''' Check the solver
        '''
        
        # check the stability condition (4-th order FD): dt <= sqrt(3/8) * dx / vmax
        dt0 = np.sqrt(3.0/8.0) * self.model.dx / np.max(self.model.vp)
        if dt0 <= self.model.dt:
            raise ValueError('the stability condition of 4-th order FD method \
            is not satisfied: dt = %.4f s > dt_required = %.4f s'.format(self.model.dt, dt0))

        # check the numerical dispersion condition: dx <= vmin/(10*f0)
        dx0 = np.min(self.model.vp) / self.source.f0 / 10.
        f00 = np.min(self.model.vp) / self.model.dx / 10.

        if dx0 < self.model.dx:
            print('Warning: modeling dispersion, dx = %.2f m > dx_required =  %.2f m' %(self.model.dx, dx0))
            print('Warning: modeling dispersion, f0 = %.2f Hz > f0_required = %.2f Hz' %(self.source.f0, f00))
  
        # check the source location
        if (self.source.coord[:,0].min() < self.model.x.min() or 
            self.source.coord[:,0].max() > self.model.x.max() or
            self.source.coord[:,1].min() < self.model.z.min() or 
            self.source.coord[:,1].max() > self.model.z.max()):
            raise ValueError('Source location out of model range.')
        
        # check the receiver
        if len(self.receiver.coord) != self.source.num:
            raise ValueError('Receiver coord list is not consistent with the source number.')

        # check the receiver location
        for rec_coord in self.receiver.coord:

            if (rec_coord[:, 0].min() < self.model.x.min() or 
                rec_coord[:, 0].max() > self.model.x.max() or
                rec_coord[:, 1].min() < self.model.z.min() or 
                rec_coord[:, 1].max() > self.model.z.max()):
                raise ValueError('Receiver location out of model range.')

        # check the mpi number and the source number
        if self.config.mpi_num > self.source.num:
            print('Warning: mpi number is larger than source number, reset to source number.')
            self.config.mpi_num = self.source.num


    def __info__(self):
        ''' Print solver information
        '''
        print('*****************************************************')
        print('\n    Seismic Waveform Inversion Toolbox: Forward    \n')
        print('*****************************************************\n')
        print('Solver: nt = {}, dt = {} ms, time = {} s'.format(self.model.nt, self.model.dt * 1000, self.model.t[-1]))
        print('Solver: nx = {}, nz = {}, dx = {} m'.format(self.model.nx, self.model.nz, self.model.dx))
        print('Solver: x  = 0 ~ {} km'.format(self.model.x[-1]/1000))
        print('Solver: z  = 0 ~ {} km'.format(self.model.z[-1]/1000))        
        print('Solver: vp = {} ~ {} m/s'.format(self.model.vp.min(), self.model.vp.max()))
        print('Solver: Source number = {}'.format(self.source.num))
        print('Solver: Receiver component = {}'.format(self.receiver.comp))
        print('Solver: {} task in parallel\n'.format(self.config.mpi_num))
