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

# import the SEISCOPE optimization toolbox
from optimization.optimization import Optimization


class Optimizer(Optimization):
    ''' Optimizer class describes the optimization algorithm
    '''

    def __init__(self, vp_init = None, rho_init = None, misfit_type = 'waveform', optim_scheme = 'CG', 
                max_update_val = 400, grad_smooth_size = 0, grad_mask = None,
                niter_max = 20, method = 'SD', 
                bound = False, lb = None, ub = None, 
                lbfgs_memory = 5, debug = False):
  

        # initialize the Optimization class
        Optimization.__init__(self, niter_max = niter_max, conv = 1e-6, method = method,
        bound = bound, lb = lb, ub = ub, lbfgs_memory = lbfgs_memory, debug = debug)

        # set the other parameters
        self.vp_init = vp_init
        self.rho_init = rho_init
        self.misfit_type = misfit_type
        self.optim_scheme = optim_scheme
        self.max_update_val = max_update_val
        self.grad_smooth_size = grad_smooth_size
        self.grad_mask = grad_mask

        # check the optimizer parameters
        self.__check__()

        # output the optimizer parameters
        self.__info__()


    def __check__(self):
        '''
            check the optimizer parameters
        '''

        # check the initial model
        if self.vp_init is None:
            raise ValueError('The initial vp model should be provided')
        if self.rho_init is None:
            raise ValueError('The initial rho model should be provided')
        if self.vp_init.shape != self.rho_init.shape:
            raise ValueError('The shape of the initial model vp and rho should be the same')

        # check the type of misfit function
        supported_misfit = ['waveform', 'envelope', 'traveltime', 'clobalcorrelation']
        if self.misfit_type.lower() not in supported_misfit:
            msg = 'Only support the following misfit types: {}'.format(supported_misfit)
            err = 'The misfit type {} is not supported.'.format(self.misfit_type)
            raise ValueError(msg + '\n' + err)
    
        # check the optimization scheme
        supported_scheme = ['sd', 'cg', 'lbfgs', 'plbfgs', 'trn', 'ptrn']

        if self.optim_scheme.lower() not in supported_scheme:
            msg = 'Only support the following optimization schemes: {}'.format(supported_scheme)
            err = 'The optimization scheme {} is not supported.'.format(self.optim_scheme)
            raise ValueError(msg + '\n' + err)

        # check the maximum update value
        if self.max_update_val <= 0 or self.max_update_val > 500:
            raise ValueError('The maximum update value should be larger than 0 and is no larger than 500 m/s (empirically)')

        #  check the maximum number of iterations
        if self.niter_max <= 0:
            raise ValueError('The maximum number of iterations should be larger than 0')

        # check the smooth radius of the gradient
        if self.grad_smooth_size < 0:
            raise ValueError('The smooth size of the gradient should be no less than 0')

        # check the mask of the gradient
        if self.grad_mask is not None:
            raise NotImplementedError('The mask of the gradient is not implemented yet')
        
        if self.grad_mask is not None and self.grad_mask.shape != self.vp_init.shape:
            raise ValueError('The shape of the mask should be the same as the initial model')
        
        if self.bound:
            if self.lb is None:
                raise ValueError('The lower bound of the model should be provided')
            if self.ub is None:
                raise ValueError('The upper bound of the model should be provided')
            
            if self.lb.max() > self.ub.min():
                raise ValueError('The lower bound should be no larger than the upper bound of the model')
            
            if self.lb.shape != self.vp_init.shape:
                raise ValueError('The shape of the lower bound should be the same as the initial model')
            
            if self.ub.shape != self.vp_init.shape:
                raise ValueError('The shape of the upper bound should be the same as the initial model')

        
    def __info__(self):
        ''' print information about the optimizer
        '''

        print('Optimizer information:')
        print('    Initial model range: {:.2f} ~ {:.2f} m/s'.format(self.vp_init.min(), self.vp_init.max()))
        print('    Misfit type: {}'.format(self.misfit_type))
        print('    Optimization scheme: {}'.format(self.optim_scheme))
        print('    Maximum update value: {:.2f} m/s'.format(self.max_update_val))
        print('    Gradient smooth size: {}'.format(self.grad_smooth_size))
        if self.grad_mask is not None:
            print('    Gradient mask is applied')
        else:
            print('    Gradient mask is not applied')
        print('    Maximum number of iterations: {}'.format(self.niter_max))
        print('    Convergence criteria: {:.2e}'.format(self.conv))
        print('    Optimization method: {}'.format(self.method))
        print('    Bound constraints option: {}'.format(self.bound))
        if self.bound:
            print('    Lower boundary range: {:.2f} ~ {:.2f} m/s'.format(self.lb.min(), self.lb.max()))
            print('    Upper boundary range: {:.2f} ~ {:.2f} m/s'.format(self.ub.min(), self.ub.max()))
        print('    Debug option: {}'.format(self.debug))
