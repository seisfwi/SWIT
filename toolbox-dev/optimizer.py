###############################################################################
#
# SWIT v1.1: Seismic Waveform Inversion Toolbox
#
#   A Python package for seismic waveform inversion
#   Developed by Haipeng Li at USTC, updated on 2022-12-21 at Stanford
#   haipengl@mail.ustc.edu.cn, haipeng@stanford.edu
#
#   Optimizer module defines the nonlinear optimization algorithm
#
###############################################################################


import numpy as np

from optimization.optimization import Optimization


class Optimizer(Optimization):
    ''' Optimizer class describes the optimization algorithm

    Parameters
    ----------
    vp_init : 2D array (float32)
        initial velocity model
    rho_init : 2D array (float32)
        initial density model (optional)
    misfit_type : str
        misfit type, 'waveform', 'envelope', 'crosscorrelation', 'globalcorrelation'
    method : str
        optimization method, 'SD', 'CG', 'LBFGS', 'PLBFGS', 'TRN', 'PTRN'
    niter_max : int
        maximum number of iterations
    bound : bool
        whether to use the bound constraint
    vp_max : float
        maximum velocity in the bound constraint
    vp_min : float
        minimum velocity in the bound constraint
    update_vpmax : float 
        allowed maximum velocity update in each iteration
    grad_smooth_size : int
        size of the Gaussian filter to smooth the gradient, 0 for no smoothing
    grad_mask : 2D array (float32)
        mask the gradient to avoid the update in certain regions denoted by 0, 
        1 for no mask (optional)
    debug : bool
        whether to output the debug information
    '''

    def __init__(self, vp_init = None, rho_init = None, 
                misfit_type = 'waveform', method = 'SD', niter_max = 20,
                bound = False, vp_max = 5000, vp_min = 1000, update_vpmax = 400, 
                grad_smooth_size = 0, grad_mask = None, debug = False):
        ''' initialize the optimizer class
        '''
        # set up the lower and upper bounds
        if bound:
            lb = np.ones_like(vp_init) * vp_min
            ub = np.ones_like(vp_init) * vp_max
        else:
            lb = None
            ub = None

        # initialize the SEISCOPE Optimization class
        Optimization.__init__(self, niter_max = niter_max, conv = 1e-6, method = method,
        bound = bound, lb = lb, ub = ub, lbfgs_memory = 5, debug = debug)

        # set other parameters
        self.vp_init = vp_init
        self.rho_init = rho_init
        self.misfit_type = misfit_type
        self.vp_max = vp_max
        self.vp_min = vp_min
        self.update_vpmax = update_vpmax
        self.grad_smooth_size = grad_smooth_size
        self.grad_mask = grad_mask

        # check the optimizer parameters
        self.__check__()

        # output the optimizer parameters
        self.__info__()


    def __check__(self):
        ''' check the optimizer parameters
        '''

        # check the initial model
        if self.vp_init is None:
            raise ValueError('The initial vp model should be provided')
        if self.rho_init is None:
            raise ValueError('The initial rho model should be provided')
        if self.vp_init.shape != self.rho_init.shape:
            raise ValueError('The shape of the initial vp and rho should be the same')

        # check the type of misfit function
        supported_misfit = ['waveform', 'envelope', 'crosscorrelation', 'globalcorrelation']
        if self.misfit_type.lower() not in supported_misfit:
            msg = 'Only support the following misfit types: {}'.format(supported_misfit)
            err = 'The misfit type {} is not supported.'.format(self.misfit_type)
            raise ValueError(msg + '\n' + err)
    
        # check the optimization scheme
        supported_scheme = ['sd', 'cg', 'lbfgs', 'plbfgs', 'trn', 'ptrn']

        if self.method.lower() not in supported_scheme:
            msg = 'Only support the following optimization schemes: {}'.format(supported_scheme)
            err = 'The optimization scheme {} is not supported.'.format(self.method)
            raise ValueError(msg + '\n' + err)

        # check the maximum update value
        if self.update_vpmax <= 0 or self.update_vpmax > 500:
            raise ValueError('The maximum update value should be larger than 0 and is no larger than 500 m/s (empirically)')

        #  check the maximum number of iterations
        if self.niter_max <= 0:
            raise ValueError('The maximum number of iterations should be larger than 0')

        # check the smooth radius of the gradient
        if self.grad_smooth_size < 0:
            raise ValueError('The smooth size of the gradient should be no less than 0')

        # check the mask of the gradient
        if self.grad_mask is not None and self.grad_mask.shape != self.vp_init.shape:
            raise ValueError('The shape of the mask should be the same as the velocity model')
        
        # check the bound constraint
        if self.bound:
            if self.vp_max <= self.vp_min:
                raise ValueError('The maximum velocity should be larger than the minimum velocity')
  
        
    def __info__(self):
        ''' print information about the optimizer
        '''

        print('\nOptimizer information:')
        print('    Misfit function   : {}'.format(self.misfit_type))
        print('    Inversion method  : {}'.format(self.method))
        print('    Maximum iter      : {}'.format(self.niter_max))
        print('    Maximum update    : {} m/s'.format(self.update_vpmax))
        print('    Maximum vp bound  : {:.2f} m/s'.format(self.vp_max))
        print('    Minimum vp bound  : {:.2f} m/s'.format(self.vp_min))
        print('    Initial vp model  : {:.2f} ~ {:.2f} m/s'.format(self.vp_init.min(), self.vp_init.max()))
        print('    Gradient smooth   : {} grids Gaussian smooth'.format(self.grad_smooth_size))
        if self.grad_mask is not None:
            print('    Gradient mask     : provided by user')
        else:
            print('    Gradient mask     : by default, 10 grids on top')
        print('    Debug option      : {}'.format(self.debug))
