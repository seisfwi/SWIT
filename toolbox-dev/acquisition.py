###############################################################################
#
# SWIT v1.1: Seismic Waveform Inversion Toolbox
#
#   A Python package for seismic waveform inversion
#   Developed by Haipeng Li at USTC, updated on 2022-12-21 at Stanford
#   haipengl@mail.ustc.edu.cn, haipeng@stanford.edu
#
#   Acquisition module defines the basic classes for seismic acquisition,
#   including Config, Model, Source, Receiver, and Solver
#
###############################################################################

import os
from multiprocessing import cpu_count
import numpy as np


class Config(object):
    '''
        config class describes the configuration
    '''

    def __init__(self, path, mpi_num, fig_aspect = 1.0):
        '''
            initialize config class

            input:
                path: path to perform the modeling or inversion
                mpi_num: number of mpi processes
                fig_aspect: aspect ratio of the figure
        '''

        # read config file
        self.path = path
        self.mpi_num = mpi_num
        self.fig_aspect = fig_aspect

        # check the config parameters
        self.__check__()


    def __check__(self):
        '''
            check the config parameters
        '''
        # check path format
        if self.path[-1] != '/':
            self.path += '/'

        # check the path
        if not os.path.exists(self.path):
            print('Warning: working path {} does not exist, creating it now.'.format(self.path))
            os.makedirs(self.path)

        # check the number of mpi processes
        if self.mpi_num > cpu_count()// 2:
            print('Warning: number of mpi processes {} is larger than the number of CPUs, setting it maximum available CPUs to {}.'.format(self.mpi_num, cpu_count()//2))
            self.mpi_num = cpu_count() // 2
    

class Model(object):
    '''
        model class describes the model for wavefield simulation
    '''

    def __init__(self, nx, nz, dx, dt, nt, pml, vp, rho):
        '''
            initialize model class

            input:
                nx: number of grid points in x direction
                nz: number of grid points in z direction
                dx: grid spacing in x and z directions
                dt: time step
                nt: number of time steps
                pml: width of pml boundary
                vp: p-wave velocity model, in kg/m^3
                rho: density model, in m/s
        '''

        # basic parameters
        self.nx = nx
        self.nz = nz
        self.dx = dx
        self.dt = dt
        self.nt = nt
        self.pml = pml
        self.rho = rho.copy()
        self.vp = vp.copy()

        # pml parameters, free surface at the top
        self.nx_pml = self.nx + self.pml * 2
        self.nz_pml = self.nz + self.pml

        # coordinates of grid points
        self.x = np.arange(0, self.nx * self.dx, self.dx)
        self.z = np.arange(0, self.nz * self.dx, self.dx)

        # time axis
        self.t = np.linspace(0, self.dt*self.nt, num=self.nt, endpoint=False)

        # default setting for output wavefield
        self.save_snap = 0
        self.save_step = 10

        # check the model parameters
        self.__check__()


    def __check__(self):
        '''
            check the model parameters
        '''

        # check the dimensions of models
        if np.shape(self.vp) != (self.nx, self.nz):
            raise ValueError('The dimensions of vp are not consistant with nx = {}, nz = {}'.format(self.nx, self.nz))
       
        if np.shape(self.rho) != (self.nx, self.nz):
            raise ValueError('The dimensions of rho are not consistant with nx = {}, nz = {}'.format(self.nx, self.nz))


class Source(object):
    '''
        source class describes the source geometry
    '''

    def __init__(self, coord, wavelet, f0):
        '''
            initialize source class

            input:
                coord: coordinates of sources, 2D array
                wavelet: source wavelet
                f0: dominant frequency of the source wavelet
        '''

        # basic parameters
        self.coord = coord
        self.wavelet = wavelet
        self.f0 = f0

        # number of sources
        self.num = len(self.coord)


    def __check__(self):
        '''
            check the source parameters
        '''

        # check coordinates of sources
        if self.coord.shape[-1] != 2:
            raise ValueError('The coordinates of sources must be 2D array')

        # check the wavelet
        if self.wavelet.shape[0] != self.num:
            raise ValueError('The number of sources is not consistant with \
                the number of source wavelets')


class Receiver(object):
    '''
        receiver class describes the receiver geometry
    '''

    def __init__(self, coord, comp):
        '''
            initialize receiver class

            input:
                coord: coordinates of receivers, 2D array
                comp: components of receivers, ['vx', 'vz', 'p']
        '''

        # basic parameters
        self.coord = coord
        self.comp = comp

        # check the receiver parameters
        self.__check__()


    def __check__(self):
        '''
            check the receiver parameters
        '''
        # check coordinates of receivers
        if not isinstance(self.coord, list):
            raise ValueError('coord must be a list of the length of the source number, with a 2D array for each element')

        # components of receivers must be a list
        if not isinstance(self.comp, list):
            raise ValueError('comp must be a list')

        # check the components of receivers
        for comp in self.comp:
            if comp not in ['vx', 'vz', 'p']:
                raise ValueError('comp can only be vx, vz or p, or a combination of them')





