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
from scipy import integrate

from tools import cleandata
from multiprocessing import cpu_count


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



    def __check__(self):
        ''' Check the solver
        '''
    
        # check the source location
        if (self.source.coord[:,0].min() < self.model.x.min() or 
            self.source.coord[:,0].max() > self.model.x.max() or
            self.source.coord[:,1].min() < self.model.z.min() or 
            self.source.coord[:,1].max() > self.model.z.max()):
            raise ValueError('Source location out of model range.')
        
        # check the receiver location
        if (self.receiver.xz[:,0].min() < self.model.x.min() or
            self.receiver.xz[:,0].max() > self.model.x.max() or
            self.receiver.xz[:,1].min() < self.model.z.min() or
            self.receiver.xz[:,1].max() > self.model.z.max()):
            raise ValueError('Receiver location out of model range.')

        # check the mpiproc and set the number of CPUs
        self.config.mpiproc = min([self.config.mpiproc, self.source.n, cpu_count() // 2])
                

    def __info__(self):
        ''' Print solver information
        '''
        print('*****************************************************')
        print('\n        Seismic Waveform Inversion Toolbox         \n')
        print('*****************************************************\n')
        print('Forward modeling : nx, nz = {}, {}'.format(self.model.nx, self.model.nz) )
        print('Forward modeling : dx = {.1f} m'.format(self.model.dx))
        print('Forward modeling : dt = {.2f} ms, {} steps'.format(self.model.dt * 1000, self.model.nt))
        print('Forward modeling : %d shots run in mpi, %d CPU available'.format(self.system.mpiproc, cpu_count() // 2))



    def forward(self, simu_type = 'obs', savesnap = 0):
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
            if not os.path.exists(ifolder) and savesnap == 1:
                os.system('mkdir %s' % ifolder)

        # prepare the forward source
        src = integrate.cumtrapz(self.source.wavelet, axis=-1, initial=0)

        # write parameters and model files
        self.write_parfile(simu_type, src, savesnap=savesnap)

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


    def adjoint(self, savesnap = 0):
        ''' Adjoint solver
        '''
    
        simu_type = 'adj'
        
        # change path and clean previous data (always, even empty)
        cleandata(self.config.homepath + 'data/adj/')

        # create working directory
        for isrc in range(self.source.srcn):
            ifolder = self.config.homepath + 'data/%s/src%d_snapshot'%(simu_type, isrc+1)
            if not os.path.exists(ifolder) and savesnap == 1:
                os.system('mkdir %s' % ifolder)

        # prepare the forward source
        src = integrate.cumtrapz(self.source.wavelet, axis=-1, initial=0)

        # write parameters and model files
        self.write_parfile(simu_type, src, savesnap=savesnap)

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

  
    def dot_product(self, simu_type = 'obs'):
        ''' Calculate dot product of two wavefields
        '''
        pass
