###############################################################################
#
# SWIT v1.1: Seismic Waveform Inversion Toolbox
#
#   A Python package for seismic waveform inversion
#   Developed by Haipeng Li at USTC, updated on 2022-12-21 at Stanford
#   haipengl@mail.ustc.edu.cn, haipeng@stanford.edu
#
#   Base module, defining the basic classes and functions
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





class optimizer(object):
    '''
        optimizer class describes the optimization algorithm
    '''

    def __init__(self, misfit_type, optim_scheme, max_iter, step_length, 
                vp_max, vp_min, acqusition_type, grad_mute_size, 
                grad_smooth_radius, grad_mask = None):
        '''
            initialize optimizer class

            input:
                misfit_type: type of misfit function
                optim_scheme: optimization scheme
                max_iter: maximum number of iterations
                step_length: step length for gradient descent
                vp_max: maximum value of p-wave velocity
                vp_min: minimum value of p-wave velocity
                acqusition_type:  marine or land
                grad_mute_size: mute the gradient at the top of the model
                grad_smooth_radius: smooth the gradient with Gaussian filter
                grad_mask: mask the gradient, 2D array with 0 and 1
        '''

        self.misfit_type = misfit_type
        self.optim_scheme = optim_scheme
        self.max_iter = max_iter
        self.step_length = step_length
        self.vp_max = vp_max
        self.vp_min = vp_min
        self.acqusition_type = acqusition_type
        self.grad_mute_size = grad_mute_size
        self.grad_smooth_radius = grad_smooth_radius
        self.grad_mask = grad_mask


    def __check__(self):
        '''
            check the optimizer parameters
        '''

        # check the type of misfit function
        if self.misfit_type.lower() not in ['waveform', 'envelope', 'traveltime', 'clobalcorrelation']:
            raise ValueError('The misfit function can be waveform, envelope, traveltime or clobalcorrelation \n \
            Not supported misfit function: {}'.format(self.misfit_type))

        # check the optimization scheme
        if self.optim_scheme.lower() not in ['descent', 'nlcg', 'lbfgs']:
            raise ValueError('The optimization scheme can be descent, nlcg or lbfgs \n \
            Not supported optimization scheme: {}'.format(self.optim_scheme))

        #  check the maximum number of iterations
        if self.max_iter <= 0:
            raise ValueError('The maximum number of iterations should be larger than 0')
        
        # check the step length
        if self.step_length <= 0 and self.step_length > 0.1:
            raise ValueError('The step length should be larger than 0, and is recomanded to be smaller than 0.1')

        # check the maximum and minimum value of p-wave velocity
        if self.vp_max <= self.vp_min:
            raise ValueError('The maximum value of p-wave velocity should be larger than the minimum value')

        # check the acquisition type
        if self.acqusition_type.lower() not in ['marine', 'land']:
            raise ValueError('The acquisition type can be marine or land \n \
            Not supported acquisition type: {}'.format(self.acqusition_type))

        # check the mute size of the gradient
        if self.grad_mute_size < 0:
            raise ValueError('The mute size of the gradient should be no less than 0')
        
        # check the smooth radius of the gradient
        if self.grad_smooth_radius < 0:
            raise ValueError('The smooth radius of the gradient should be no less than 0')

        # check the mask of the gradient
        if self.grad_mask is not None:
            raise NotImplementedError('The mask of the gradient is not implemented yet')