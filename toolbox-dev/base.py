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


class config(object):
    '''
        config class describes the configuration
    '''

    def __init__(self, path, mpi_num, fig_aspect):
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


class model(object):
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
        self.rho = rho
        self.vp = vp

        # pml parameters
        self.nx_pml = self.nx + self.pml * 2
        self.nz_pml = self.nz + self.pml  # free surface at the top

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
            raise ValueError('The dimensions of vp are not consistant with \
                nx = {}, nz = {}'.format(self.nx, self.nz))
       
        if np.shape(self.rho) != (self.nx, self.nz):
            raise ValueError('The dimensions of rho are not consistant with \
                nx = {}, nz = {}'.format(self.nx, self.nz))

        # check path format
        if self.path[-1] != '/':
            self.path += '/'

        # check the stability condition of 4-th order FD method: dt <= sqrt(3/8) * dx / vmax
        dt_req = np.sqrt(3.0/8.0) * self.dx / np.max(self.vp)

        if dt_req <= self.dt:
            raise ValueError('the stability condition of 4-th order FD method \
            is not satisfied: dt = %.4f ms > dt_required = %.4f ms'.format(dt*1000, dt_req*1000))

        # check the numerical dispersion condition: dx <= vmin/(10*f0)
        dx_req = np.min(self.vp) / self.f0 / 10.
        f0_req = np.min(self.vp) / self.dx / 10.

        if dx_req < self.dx:
            print('Warning: modeling dispersion, dx = %.2f m > dx_required =  %.2f m' %(self.dx, dx_req))
            print('Warning: modeling dispersion, f0 = %.2f Hz > f0_required = %.2f Hz' %(self.f0, f0_req))
  


class receiver(object):
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

        # number of receivers
        self.num = len(self.coord)

        # check the receiver parameters
        self.__check__()


    def __check__(self):
        '''
            check the receiver parameters
        '''

        # check the dimensions of models
        if self.rec_comp not in ['vx', 'vz', 'p']:
            raise ValueError('rec_comp can only be vx, vz or p, or a combination of them')

        
class source(object):
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

        # check the dimensions of models
        if self.wavelet.shape[0] != self.num:
            raise ValueError('The number of sources is not consistant with \
                the number of source wavelets')


class preprocessor(object):
    ''' 
        preprocessor class describes the data preprocessing
    '''

    def __init__(self, fre_filter, fre_low, fre_high, 
                 mute_late_arrival, mute_late_window, normalize,
                 mute_offset_short, mute_offset_long, 
                 mute_offset_short_dis, mute_offset_long_dis):
        '''
            initialize source class

            input:
                fre_filter: filter the data or not
                fre_low: low frequency of the bandpass filter
                fre_high: high frequency of the bandpass filter
                mute_late_arrival: mute late arrivals or not
                mute_late_window: window length for late arrival mute
                normalize: normalize the data or not
                mute_offset_short: mute short offset or not
                mute_offset_long: mute long offset or not
                mute_offset_short_dis: short offset distance    (units: m)
                mute_offset_long_dis: long offset distance      (units: m)
        '''

        # data filter
        self.fre_filter = fre_filter
        self.fre_low = fre_low
        self.fre_high = fre_high

        # pick first break and mute later arrivals
        self.mute_late_arrival = mute_late_arrival
        self.mute_late_window = mute_late_window

        # data normalization
        self.normalize = normalize

        # data offset mute
        self.mute_offset_short = mute_offset_short
        self.mute_offset_long = mute_offset_long
        self.mute_offset_short_dis = mute_offset_short_dis           
        self.mute_offset_long_dis = mute_offset_long_dis             


    def __check__(self):

        # check the data filter
        if self.fre_filter not in ['Lowpass', 'Bandpass', 'Highpass', 'None']:
            raise ValueError(' The filter type is Lowpass, Bandpass, Highpass or None\n  \
                Not supported frequency filter: {}'.format(self.fre_filter))

        if self.fre_low > self.fre_high:
            raise ValueError('The low frequency is larger than the high frequency')

        # check the late arrival mute
        if self.mute_late_arrival and self.mute_late_window <= 0:
            raise ValueError('The late arrival mute window should be larger than 0')

        # check the data normalization
        if self.normalize not in ['Max-Trace', 'L1-Event', 'L2-Event', 'L1-Trace', 'L2-Trace', 'None']:
            raise ValueError('The normalization type is Max-Trace, L1-Event, L2-Event, L1-Trace, L2-Trace or None \n \
                Not supported normalization type: {}'.format(self.normalize))

        # check the data offset mute
        if self.mute_offset_short_dis > self.mute_offset_long_dis:
            raise ValueError('The short offset distance is larger than the long offset distance')



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