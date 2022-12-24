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


    def forward(self, simu_type = 'obs', save_snap = 0):
        ''' Forward solver
        '''

        if simu_type not in ['obs', 'syn']:
            raise ValueError('Forward: unsupport simulation type: {}.'.format(simu_type))

        # change path and clean previous data (always, even empty)
        if simu_type in ['syn']: 
            cleandata(self.config.path + 'data/syn/')

        # create working directory
        for isrc in range(self.source.srcn):
            ifolder = self.config.path + 'data/%s/src%d_snapshot'%(simu_type, isrc+1)
            if not os.path.exists(ifolder) and save_snap == 1:
                os.system('mkdir %s' % ifolder)

        # prepare the forward source
        src = integrate.cumtrapz(self.source.wavelet, axis=-1, initial=0)

        # write parameters and model files
        self.write_parfile(simu_type, src, save_snap=save_snap)

        # check mpi and wavefield solver
        if os.system('which mpirun') != 0:
            raise ValueError('Cannot find mpirun, please check your mpi installation.')
        if os.system('which fd2dmpi') != 0:
            raise ValueError('Cannot find fd2dmpi, please check your wavefield solver installation.')

        # submit job
        os.chdir(self.config.homepath)
        solver_cmd = 'mpirun -np %d  fd2dmpi par=%s' % (self.config.mpiproc, 
                    self.config.homepath + 'parfile/forward_parfile/parfile')
        status = subprocess.getstatusoutput(solver_cmd)
        if status[0]:
            print(status[1])
            raise ValueError('Forward solver crash')


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
        
        # check the receiver location
        if (self.receiver.coord[:,0].min() < self.model.x.min() or
            self.receiver.coord[:,0].max() > self.model.x.max() or
            self.receiver.coord[:,1].min() < self.model.z.min() or
            self.receiver.coord[:,1].max() > self.model.z.max()):
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
