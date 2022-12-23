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


class Solver(object):
    ''' Acoustic Wavefiled Solver class
    '''

    def __init__(self, config, model, source, receiver):
        ''' Initialize solver
        '''
        self.config = config
        self.model = model
        self.source = source
        self.receiver = receiver


    def forward(self, simu_type = 'obs', savesnap = 0):
        ''' Forward solver
        '''

        if simu_type not in ['obs', 'syn']:
            raise ValueError('Forward: unsupport simulation type: {}.'.format(simu_type))

        # change path and clean previous data (always, even empty)
        if simu_type in ['syn']: 
            cleandata(self.config.homepath + 'data/syn/')

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
