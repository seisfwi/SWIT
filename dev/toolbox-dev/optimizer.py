###############################################################################
#
# SWIT v1.1: Seismic Waveform Inversion Toolbox
#
#   A Python package for seismic waveform inversion
#   By Haipeng Li at USTC & Stanford
#   Email: haipengl@mail.ustc.edu.cn, haipeng@stanford.edu 
#
#   Optimizer class for data preprocessing
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
            optimization method, 'SD', 'CG', 'LBFGS', 'PLBFGS'
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
            1 for no mask
        debug : bool
            whether to output the debug information
    '''

    def __init__(self, vp_init = None, rho_init = None, 
                misfit_type = 'waveform', method = 'SD', niter_max = 20,
                bound = False, vp_max = 5000, vp_min = 1000, update_vpmax = 400, 
                grad_smooth_size = 0, grad_mask = None, debug = False, use_default_mask = False):
        ''' initialize the optimizer class
        '''

        # set up the bound constraint if bound is True
        lb = vp_min * np.ones_like(vp_init).flatten() if bound else None
        ub = vp_max * np.ones_like(vp_init).flatten() if bound else None

        # inherit the SEISCOPE Optimization class
        Optimization.__init__(self, niter_max = niter_max, conv = 1e-6, 
                             method = method, bound = bound, lb = lb, ub = ub, 
                             lbfgs_memory = 5, debug = debug)

        # set other parameters
        self.vp_init = vp_init
        self.rho_init = rho_init
        self.misfit_type = misfit_type
        self.vp_max = vp_max
        self.vp_min = vp_min
        self.update_vpmax = update_vpmax
        self.grad_smooth_size = grad_smooth_size
        self.grad_mask = grad_mask
        self.use_default_mask = use_default_mask

        # check the optimizer parameters
        self.__check__()

        # output the optimizer parameters
        self.__info__()


    def __check__(self):
        ''' check the optimizer parameters
        '''

        # check the initial model
        if self.vp_init is None or self.rho_init is None or self.grad_mask is None:
            raise ValueError('Optimizer: the initial vp and rho models and gradient mask should be provided')

        # check the dimension of the initial models and gradient mask
        if self.vp_init.shape != self.rho_init.shape or self.vp_init.shape != self.grad_mask.shape:
            msg  = 'Optimizer: the shape of the initial vp model is {}\n'.format(self.vp_init.shape)
            msg += 'Optimizer: the shape of the initial rho model is {}\n'.format(self.rho_init.shape)
            err = 'Optimizer: the shape of the gradient mask is {}\n'.format(self.grad_mask.shape)
            err += 'Optimizer: the shape of the initial vp and rho models and gradient mask should be the same'
            raise ValueError(msg + err)

        # check the type of misfit function
        misfits = ['waveform', 'envelope', 'crosscorrelation', 'globalcorrelation', 'hybrid']
        if self.misfit_type.lower() not in misfits:
            msg = 'Optimizer: misfit_type must be one of {}'.format(misfits)
            err = 'Optimizer: unknown simulation type: {}'.format(self.misfit_type)
            raise ValueError(msg + '\n' + err)
    
        # check the optimization method
        methods = ['sd', 'cg', 'lbfgs', 'plbfgs']
        if self.method.lower() not in methods:
            msg = 'Optimizer: method must be one of {}'.format(methods)
            err = 'Optimizer: unknown simulation type: {}'.format(self.method)
            raise ValueError(msg + '\n' + err)

        # check the maximum update value
        if self.update_vpmax <= 0 or self.update_vpmax > 500:
            raise ValueError('Optimizer: the maximum update value should be larger than 0 and is no larger than 500 m/s (empirically)')

        # check the maximum number of iterations
        if self.niter_max <= 0:
            raise ValueError('Optimizer: the maximum number of iterations should be larger than 0')

        # check the smooth radius of the gradient
        if self.grad_smooth_size < 0:
            raise ValueError('Optimizer: the smooth size of the gradient should be no less than 0')
     
        # check the bound constraint
        if self.bound:
            if self.vp_max <= self.vp_min:
                raise ValueError('Optimizer: the maximum velocity should be larger than the minimum velocity')
            
            if self.vp_max < self.vp_init.max() or self.vp_min > self.vp_init.min():
                msg  = '\nOptimizer: the range of initial model is {:.2f} ~ {:.2f} m/s\n'.format(self.vp_init.min(), self.vp_init.max())
                msg += 'Optimizer: the range of bound constraint is {:.2f} ~ {:.2f} m/s\n'.format(self.vp_min, self.vp_max)
                err = 'Optimizer: the bound constraint should be larger than the range of initial model\n'
                raise ValueError(msg + err)


        
    def __info__(self):
        ''' print information about the optimizer
        '''

        print('Optimizer information:')
        print('    Misfit function   : {}'.format(self.misfit_type))
        print('    Inversion method  : {}'.format(self.method))
        print('    Maximum iter      : {}'.format(self.niter_max))
        print('    Maximum update    : {} m/s'.format(self.update_vpmax))
        if self.bound:
            print('    Bound constraint  : {:.2f} ~ {:.2f} m/s'.format(self.vp_min, self.vp_max))
        else: 
            print('    Bound constraint  : {}'.format(self.bound))
        print('    Initial vp model  : {:.2f} ~ {:.2f} m/s'.format(self.vp_init.min(), self.vp_init.max()))
        print('    Gradient smooth   : {} grids Gaussian smooth'.format(self.grad_smooth_size))
        if self.use_default_mask:
            print('    Gradient mask     : default damping mask on top of the model (10 grids)')
        else:
            print('    Gradient mask     : provided by user')
        print('    Debug option      : {}'.format(self.debug))
