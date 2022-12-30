###############################################################################
#
# SWIT v1.1: Seismic Waveform Inversion Toolbox
#
#   A Python package for seismic waveform inversion
#   Developed by Haipeng Li at USTC, updated on 2022-12-21 at Stanford
#   haipengl@mail.ustc.edu.cn, haipeng@stanford.edu
#
#   Optimizer module: a self-contained optimization package
# 
#   This code is reimplemented by Haipeng Li in Python based on the SEISCOPE 
#   optimization toolbox (https://seiscope2.osug.fr/SEISCOPE-OPTIMIZATION-TOOLBOX-274). 
#   The original code is written by Ludovic Métivier and Romain Modrak: 
#
#   Métivier, L., & Brossier, R. (2016). The SEISCOPE optimization toolbox: A large-scale 
#   nonlinear optimization library based on reverse communication. Geophysics, 81(2), F1-F15.
# 
#   Note that some of the codes and default parameters are modified to fit the SWIT package.
#
###############################################################################

import numpy as np


class Optimization(object):
    '''
        Optimization based on Fortran subroutines from SEISCOPE OPTIMIZATION TOOLBOX
    '''

    def __init__(self, niter_max = 20, conv = 1e-8, method = 'SD', bound = False, lb = None, ub = None, lbfgs_memory = 20, debug = False):
        '''
            Initialize the optimizer class

            Input:
                niter_max: maximum number of iterations
                conv: convergence criteria
                method: optimization method
                bound: bound constraints option
                lb: lower boundary
                ub: upper boundary
                lbfgs_memory: LBFGS memory
                debug: debug option FOR USER
        '''

        self.niter_max = niter_max
        self.conv = conv
        self.method= method
        self.bound = bound
        self.lb = lb
        self.ub = ub
        self.debug = debug

        # Default parameters
        self.threshold = 100.           # Tolerance on bound constraints satisfaction
        self.FLAG = 'INIT'

        # Initialize linesearch parameters by default
        self.m1 = 0.0                   # Wolfe conditions parameter 1, 0.0 for no Wolfe conditions due to the scale of the misfit
        self.m2 = 0.95                  # Wolfe conditions parameter 2
        self.mult_factor = 5            # bracketting parameter
        self.nls_max = 10               # max number of linesearch
        self.cpt_ls = 0                 # linesearch counter
        self.first_ls = True            # first linesearch flag
        self.alpha = 1.                 # first value for the linesearch steplength

        # Quasi-Newton l-BFGS method
        if self.method in ['LBFGS', 'PLBFGS']:
            self.l = lbfgs_memory 

        # Truncated Newton method
        if self.method in ['TRN', 'PTRN']:
            self.niter_max_CG=5         # maximum number of inner conjugate gradient iterations


    def iterate(self, x, fcost, grad, grad_preco):
        ''' perform nonlinear inversion
        '''

        # copy the input arrays and convert them to 1D arrays
        x = np.copy(x.flatten())
        grad = np.copy(grad.flatten())
        grad_preco = np.copy(grad_preco.flatten())

        # Preconditioned Steepest Descent: PSTD
        if self.method == 'SD':
            x = self.PSTD(x, fcost, grad, grad_preco)

        # Preconditioned nonlinear conjugate gradient: PNLCG
        elif self.method == 'CG':
            x = self.PNLCG(x, fcost, grad, grad_preco)

        # Quasi-Newton l-BFGS method: LBFGS
        elif self.method == 'LBFGS':
            x = self.LBFGS(x, fcost, grad)

        # Preconditioned Quasi-Newton l-BFGS method: PLBFGS
        elif self.method == 'PLBFGS':
            x = self.PLBFGS(x, fcost, grad, grad_preco)
        
        # Truncated Newton method: TRN
        elif self.method == 'TRN':
            x = self.TRN(x, fcost, grad)

        # Preconditioned Truncated Newton method: TRN
        elif self.method == 'PTRN':
            x = self.PTRN(x, fcost, grad, grad_preco)

        else:
            raise AssertionError('method: SD, CG, LBFGS, PLBFGS, TRN, PTRN')

        return x


    #--------------------------------------------------------------------#
    #            Preconditioned Steepest Descent: PSTD                   # 
    #--------------------------------------------------------------------#
    def PSTD(self, x, fcost, grad, grad_preco):
        ''' Preconditioned Steepest Decent method
        '''

        # Initialize the linesearch process
        if(self.FLAG == 'INIT'):
            
            self.init_PSTD(x, fcost, grad, grad_preco)
            x = self.linesearch(x, fcost, grad)

            self.print_info(fcost)
            self.FLAG = 'GRAD'
            self.nfwd_pb += 1
        
        # call the linesearch process
        else:
            x = self.linesearch(x, fcost, grad)

            if(self.task == 'NEW_STEP'):
                self.cpt_iter += 1
                
                # test for convergence
                if(self.std_test_conv(fcost)):
                    self.FLAG = 'CONV'
                    self.print_info(fcost)
                
                # if a NEW_STEP is taken, compute a new descent direction using 
                # current descent, gradient and preconditioned gradient
                else:
                    self.FLAG = 'NSTE'
                    self.grad = grad
                    self.descent = -1. * grad_preco
                    self.print_info(fcost)

            # if the linesearch needs a new gradient then ask the user to provide it
            elif(self.task == 'NEW_GRAD'):
                self.FLAG = 'GRAD'
                self.nfwd_pb = self.nfwd_pb+1
            
            # if the linesearch has failed, inform the user
            elif(self.task == 'FAILURE!'):
                self.FLAG = 'FAIL'
                self.print_info(fcost)
        return x


    def init_PSTD(self, x, fcost, grad, grad_preco):
        '''  Initilize linesearch parameters
        '''

        # set counters
        self.cpt_iter = 0
        self.f0 = fcost
        self.nfwd_pb = 0

        # linesearch parameters
        self.fk = fcost

        # memory allocations
        self.xk = np.zeros_like(x)
        self.grad = np.zeros_like(x)
        self.descent = np.zeros_like(x)

        self.xk[:] = x[:]
        self.grad[:] = grad[:]

        # first descent direction
        self.descent[:] = -1. * grad_preco[:]


    #--------------------------------------------------------------------#
    #       Preconditioned nonlinear conjugate gradient: PNLCG           # 
    #--------------------------------------------------------------------#
    def PNLCG(self, x, fcost, grad, grad_preco):
        ''' Preconditioned nonlinear conjugate gradient
        '''
    
        # initialize the linesearch process
        if(self.FLAG == 'INIT'):

            self.init_PNLCG(x, fcost, grad, grad_preco)
            x = self.linesearch(x, fcost, grad)

            self.print_info(fcost)
            self.FLAG = 'GRAD'
            self.nfwd_pb += 1

            # Store current gradient before the user compute the new one
            self.grad_prev[:] = grad[:]
        
        # call the linesearch process
        else:
            x = self.linesearch(x, fcost, grad)

            if(self.task == 'NEW_STEP'): 
                self.cpt_iter += 1

                # test for convergence
                if(self.std_test_conv(fcost)):
                    self.FLAG = 'CONV'
                    self.grad[:] = grad[:]
                    self.print_info(fcost)
                
                # if a NEW_STEP is taken, compute a new descent direction using 
                # current descent, gradient and preconditioned gradient
                else:
                    self.FLAG = 'NSTE'
                    self.descent_PNLCG(grad, grad_preco)                    
                    self.grad[:] = grad[:]
                    self.print_info(fcost)

            # if the linesearch needs a new gradient then ask the user to provide it
            elif(self.task == 'NEW_GRAD'):
                self.FLAG = 'GRAD'
                self.nfwd_pb += 1
                
                # Store current gradient before the user compute the new one
                self.grad_prev[:] = grad[:]

            # if the linesearch has failed, inform the user
            elif(self.task == 'FAILURE!'):
                self.FLAG ='FAIL'
                self.grad[:] = grad[:]
                self.print_info(fcost)

        return x


    def init_PNLCG(self, x, fcost, grad, grad_preco):
        '''  Initilize linesearch parameters
        '''

        # set counters
        self.cpt_iter = 0
        self.f0 = fcost
        self.nfwd_pb = 0

        # linesearch parameters
        self.fk=fcost

        # memory allocations
        self.grad_prev = np.zeros_like(x)
        self.descent_prev = np.zeros_like(x)
        self.xk = np.zeros_like(x)
        self.grad = np.zeros_like(x)
        self.descent = np.zeros_like(x)

        self.xk[:] = x[:]
        self.grad[:] = grad[:]
            
        # first descent direction
        self.descent[:] = -1. * grad_preco[:]


    def descent_PNLCG(self, grad, grad_preco):
        ''' Descent direction
        '''

        # Storing old descent direction
        self.descent_prev[:] = self.descent[:]
        
        # The original PNLCG method (Dai and Yuan, 1999) in the SEISCOPE optimization toolbox
        # # Computation of beta (Dai and Yuan, 1999)
        # gkpgk = scalL2(grad, grad_preco)
        # skpk  = scalL2(grad - self.grad_prev, self.descent_prev)
        # if skpk != 0.:
        #     beta = gkpgk / skpk
        # else:
        #     beta = 0.

        # # Safeguard (may be useful in some cases)
        # if((beta >= 1e5) or (beta <= -1e5)):
        #     beta = 0.

        # Computation of beta (Polak and Ribière 1969), modified by Haipeng Li
        gkpgk = scalL2(grad, grad - self.grad_prev)
        skpk  = scalL2(self.grad_prev, self.grad_prev)
        if skpk != 0.:
            beta = gkpgk / skpk
        else:
            beta = 0.

        # Safeguard for beta (may be useful in some cases, empirical values)
        if beta > 2 or beta < 0:
            # print('Optimization: reset to descent direction in NLCG due to illegal beta = {:.2f}'.format(beta))
            beta = 0.
            
        # Computation of the descent direction
        self.descent[:] = -1. * grad_preco[:] + beta * self.descent_prev[:]
        

    #--------------------------------------------------------------------#
    #                Quasi-Newton l-BFGS method: LBFGS                   # 
    #--------------------------------------------------------------------#
    def LBFGS(self, x, fcost, grad):
        ''' Quasi-Newton l-BFGS method
        '''
    
        if(self.FLAG == 'INIT'):
            # initialize the linesearch process
            self.init_LBFGS(x, fcost, grad)
            x = self.linesearch(x, fcost, grad)
            self.print_info(fcost)
            self.FLAG = 'GRAD'
            self.nfwd_pb += 1
        else:
            # else call the linesearch process
            x = self.linesearch(x, fcost, grad)
            if(self.task == 'NEW_STEP'):
                self.cpt_iter += 1
                # test for convergence
                if(self.std_test_conv(fcost)):
                    self.FLAG = 'CONV'
                    self.grad[:] = grad[:]
                    self.print_info(fcost)
                else:
                    self.FLAG = 'NSTE'
                    # if a NEW_STEP is taken, compute a new descent direction using current gradient
                    # and l-BFGS approximation of the inverse Hessian preconditioned gradient
           
                    # LBFGS update        
                    self.update_LBFGS(x, grad)
                    
                    # Computation of the new descent direction
                    self.descent_LBFGS(grad)
                    
                    # LBFGS store
                    self.save_LBFGS(x, grad)
                    
                    # print info on current iteration
                    self.grad[:] = grad[:]
                    self.print_info(fcost)
                    
            elif(self.task == 'NEW_GRAD'):
                # if the linesearch needs a new gradient then ask the user to provide it
                self.FLAG = 'GRAD'
                self.nfwd_pb += 1

            elif(self.task == 'FAILURE!'):
                # if the linesearch has failed, inform the user       #
                self.FLAG='FAIL'
                # print info on current iteration
                self.grad[:] = grad[:]
                self.print_info(fcost)
        return x


    def init_LBFGS(self, x, fcost, grad):
        '''  Initilize linesearch parameters
        '''
        
        # set counters
        n = np.size(x)
        self.cpt_iter = 0
        self.f0 = fcost
        self.nfwd_pb = 0  
        self.sk = np.zeros((n, self.l))
        self.yk = np.zeros((n, self.l))
        self.cpt_lbfgs = 1
    
        # initialize linesearch parameters
        self.fk = fcost

        # memory allocations
        self.xk = np.zeros_like(x)
        self.grad = np.zeros_like(x)
        self.descent = np.zeros_like(x)

        self.xk[:] = x[:]
        self.grad[:] = grad[:]

        # first descent direction
        self.descent[:] = -1. * grad[:]
        
        # LBFGS save
        self.save_LBFGS(x, grad)
        
    
    def save_LBFGS(self, x, grad):
        ''' save values for LBFGS
        '''
        
        if(self.cpt_lbfgs <= self.l):
            # if the number of stored pairs does not exceed the maximum value, then save x and grad
            self.sk[:, self.cpt_lbfgs-1] = x[:]
            self.yk[:, self.cpt_lbfgs-1] = grad[:]
        else:
            # otherwise, erase the oldest pair and save the new one (shift)                         
            for i in range (self.l-1):
                self.sk[:, i] = self.sk[:, i+1]
                self.yk[:, i] = self.yk[:, i+1]
        
            self.sk[:, self.l-1] = x[:]
            self.yk[:, self.l-1] = grad[:]
        

    def update_LBFGS(self, x, grad):
        ''' update values for LBFGS
        '''

        if(self.cpt_lbfgs <= self.l):
            # if the number of stored pairs does not exceed the maximum value, 
            # then compute a new pair sk yk and update the counter cpt_lbfgs
            ii = self.cpt_lbfgs - 1
            self.sk[:, ii] = x    - self.sk[:, ii]
            self.yk[:, ii] = grad - self.yk[:, ii]
            self.cpt_lbfgs += 1
        else:
            # otherwise, simply update the lth pair
            self.sk[:,self.l-1] = x[:]    - self.sk[:,self.l-1]
            self.yk[:,self.l-1] = grad[:] - self.yk[:,self.l-1]


    def descent_LBFGS(self, grad):
        ''' Descent direction
        '''
        
        n = np.size(grad)
        borne_i = self.cpt_lbfgs - 1

        # Safeguard
        norml2_sk = normL2(self.sk[:,borne_i-1])
        norml2_yk = normL2(self.yk[:,borne_i-1])
        if( (norml2_sk == 0.) or (norml2_yk == 0.)):
            self.descent = -1. * grad
        else:
            # First phase of the recursion loop
            alpha = np.zeros(self.cpt_lbfgs)
            rho   = np.zeros(self.cpt_lbfgs)
            q     = np.zeros(n)
            q[:]  = grad[:]

            for i in range(borne_i):
                ii = borne_i - i - 1
                rho[ii] = 1. / scalL2(self.yk[:, ii], self.sk[:,ii])
                alpha[ii] = rho[ii] * scalL2(self.sk[:,ii], q)
                q[:] = q[:] - alpha[ii] * self.yk[:,ii]
        
            gamma_num = scalL2(self.sk[:,borne_i-1], self.yk[:,borne_i-1])
            gamma_den = normL2(self.yk[:,borne_i-1])

            # Scaling by gamma
            gamma = gamma_num / (gamma_den * gamma_den)
            self.descent[:] = gamma * q[:]

            # Second phase of the recursion loop
            for i in range (borne_i):
                beta = rho[i] * scalL2(self.yk[:,i], self.descent)
                self.descent[:] += (alpha[i] - beta) * self.sk[:, i]

            self.descent *= -1.


    #--------------------------------------------------------------------#
    #        Quasi-Newton preconditioned l-BFGS method: PLBFGS           # 
    #--------------------------------------------------------------------#
    def PLBFGS(self, x, fcost, grad, grad_preco):
        ''' Quasi-Newton preconditioned l-BFGS method
        '''

        if(self.FLAG == 'INIT'):
            # initialize the linesearch process
            self.init_PLBFGS(x, fcost, grad, grad_preco)
            x = self.linesearch(x, fcost, grad)
            self.print_info(fcost)
            self.FLAG = 'GRAD'
            self.nfwd_pb += 1

        elif(self.FLAG == 'PREC'):
            # if self.FLAG is PREC, we return from a call to the user preconditioner,
            # we have to finish the computation of the descent direction
            self.descent2_PLBFGS()

            # LBFGS save
            self.save_LBFGS(x, grad)
            self.cpt_iter += 1
            
            # before continuing we test for convergence
            if(self.std_test_conv(fcost)):
                self.FLAG = 'CONV'
                # print info on current iteration
                self.grad[:] = grad[:]
                self.print_info(fcost)
            else:
                self.FLAG = 'NSTE'
                # print info on current iteration
                self.grad[:] = grad[:]
                self.print_info(fcost)
        else:
            # else call the linesearch process
            x = self.linesearch(x, fcost, grad)
            if(self.task == 'NEW_STEP'): 
                # LBFGS update
                self.update_LBFGS(x, grad)
                # Start the computation of the new descent direction
                self.descent1_PLBFGS(grad)
                # Set self.FLAG to PREC for asking user to perform preconditioning 
                self.FLAG = 'PREC'                        
            elif(self.task == 'NEW_GRAD'): 
                # if the linesearch needs a new gradient: ask the user to provide it
                self.FLAG = 'GRAD'
                self.nfwd_pb += 1
            elif(self.task == 'FAILURE!'):
                # if the linesearch has failed, inform the user
                self.FLAG = 'FAIL'
                # print info on current iteration
                self.grad[:] = grad[:]
                self.print_info(fcost)
        return x

    def init_PLBFGS(self, x, fcost, grad, grad_preco):
        '''  Initilize linesearch parameters
        '''

        # set counters
        self.cpt_iter = 0
        self.f0 = fcost
        self.nfwd_pb = 0
        self.cpt_lbfgs = 1
        
        # initialize linesearch parameters
        self.fk = fcost

        # memory allocations
        n = np.size(x)

        self.xk = np.zeros_like(x)
        self.grad = np.zeros_like(x)
        self.descent = np.zeros_like(x)
        self.sk = np.zeros((n, self.l))
        self.yk = np.zeros((n, self.l))
        self.q_plb = np.zeros_like(x)

        self.xk[:] = x[:]
        self.grad[:] = grad[:]
        
        # first descent direction
        self.descent[:] = -1. * grad_preco[:]
        self.save_LBFGS(x, grad)


    # First loop
    def descent1_PLBFGS(self, grad):
        ''' First loop for PLBFGS
        '''

        borne_i = self.cpt_lbfgs - 1
        self.alpha_plb = np.zeros(self.cpt_lbfgs)
        self.rho_plb = np.zeros(self.cpt_lbfgs) 
        self.q_plb[:] = grad[:]

        for i in range(borne_i):
            ii = borne_i-i-1
            self.rho_plb[ii]   = 1. / scalL2(self.yk[:, ii], self.sk[:, ii])
            self.alpha_plb[ii] = self.rho_plb[ii] * scalL2(self.sk[:, ii], self.q_plb[:]) 
            self.q_plb[:]     -= self.alpha_plb[ii] * self.yk[:,ii]

    # Second loop
    def descent2_PLBFGS(self):
        ''' Second loop for PLBFGS
        '''
    
        borne_i = self.cpt_lbfgs-1  
        gamma_num = scalL2(self.sk[:, borne_i-1], self.yk[:, borne_i-1])
        gamma_den = normL2(self.yk[:, borne_i-1])
        gamma = gamma_num / (gamma_den**2)
        self.descent[:] = gamma * self.q_plb[:]
        for i in range(borne_i):
            beta = self.rho_plb[i] * scalL2(self.yk[:, i], self.descent[:])
            self.descent[:] += (self.alpha_plb[i]-beta) * self.sk[:,i]
        self.descent[:] *= -1.


    #--------------------------------------------------------------------#
    #                    Truncated Newton method: TRN                    # 
    #--------------------------------------------------------------------#
    def TRN(self, x, fcost, grad):
        ''' Truncated Newton method
        '''
    
        if(self.FLAG=='INIT'):
            # if FLAG is INIT, call the dedicated initialization def to allocate data structure self
            self.init_TRN( x, fcost, grad)
            self.print_info_TRN(fcost)
            self.comm = 'DESC'
            self.CG_phase = 'INIT'
            self.nfwd_pb = self.nfwd_pb+1
            self.conv_CG = False
            self.FLAG = 'NONE'

        if(self.comm == 'DESC'):
            # if self.comm is DES, the selfizer is computing a descent direction through the conjugate gradient
            self.descent_TRN(grad)
            
            if(self.conv_CG):
                # if the conjugate gradient has converged go to next phase: linesearch in the descent direction
                self.comm = 'NSTE'
                self.CG_phase = 'INIT'
                self.FLAG = 'NONE'
            else:
                # else perform a new iteration of conjugate gradient and ask the user to compute a Hessian-vector product
                self.FLAG = 'HESS'
                self.nhess = self.nhess+1
            
        elif(self.comm == 'NSTE'):
            # if self.comm is NSTE, the selfizer is looking for a new step in the descent direction
            x = self.linesearch(x, fcost, grad)
            if(self.task=='NEW_STEP'):          
                # if self.task is 'NEW_STEP, the linesearch process has found the new step
                self.cpt_iter += 1
                # Save the previous gradient norm 
                self.norm_grad_m1=self.norm_grad
                # Compute the new gradient norm 
                self.norm_grad = normL2(grad)
                # Print info on current nonlinear iteration
                self.print_info_TRN(fcost)
                # Test for convergence
                if(self.std_test_conv(fcost)):
                    self.FLAG = 'CONV'
                    self.print_info_TRN(fcost)
                else:
                    # Flags for the computation of the new descent direction
                    self.FLAG = 'NSTE'
                    self.comm = 'DESC'
                    # Update forcing term self.eta following the Eisenstat and Walker formula
                    self.forcing_term_TRN(grad)

            elif(self.task == 'NEW_GRAD'): # STILL SEARCHING THE STEP 
                # if self.task is 'NEW_GRAD, the linesearch process is continuing, the gradient at the current point is required
                self.FLAG = 'GRAD'
                self.nfwd_pb += 1

            elif(self.task =='FAILURE!'):
                # if self.task is 'FAILURE, the linesearch process has failed, the iterations are stopped
                self.FLAG = 'FAIL'
                self.print_info_TRN(fcost)

        return x

    def init_TRN(self, x, fcost, grad):
        '''  Initilize linesearch parameters
        '''
        
        # set counters
        self.cpt_iter = 0
        self.f0 = fcost   
        self.nfwd_pb = 0    
        self.nhess = 0
        self.eta = 0.9   
        #self.eta=0.1   
        #self.eta=1e-6   
            
        # initialize linesearch parameters
        self.fk = fcost
        self.cpt_iter_CG = 0
        
        # memory allocations
        self.xk = np.zeros_like(x)
        self.grad = np.zeros_like(x)
        self.descent = np.zeros_like(x)
        self.descent_prev = np.zeros_like(x)
        self.residual = np.zeros_like(x)
        self.d = np.zeros_like(x)
        self.Hd = np.zeros_like(x)
        self.eisenvect = np.zeros_like(x)

        self.xk[:] = x[:]
        self.grad[:] = grad[:]

        # norm of the first gradient
        self.norm_grad = normL2(grad)


    def descent_TRN(self, grad):
        ''' Descent direction
        '''

        if(self.CG_phase=='INIT'):
            # if self.CG_phase is INIT, initialize the conjugate gradient process
            self.residual[:] = grad[:]
            self.d[:] = -1.*self.residual[:]
            self.Hd[:] = 0.
            self.descent[:] = 0.
            self.qk_CG = 0.
            self.hessian_term = 0.
            self.norm_residual = normL2(self.residual)
            self.conv_CG = False
            self.cpt_iter_CG = 0
            self.print_info_TRN(0.0)
            self.CG_phase='IRUN'
        else:
            # else perform one conjugate gradient iteration
            dHd = scalL2(self.d, self.Hd)
            if(dHd < 0.):
                # if dHd < 0, detection of a negative eigenvalue of the Hessian operator, stop the process
                self.conv_CG = True
                print('Negative curvature')
                if(self.cpt_iter_CG == 0):
                    # if this is the first iteration, :return the opposite of the gradient as descent direction
                    self.descent[:] = self.d[:]

                    #-----------------------------------------------------#
                    # if the debug option is activated, compute the       #
                    # quadratic function minimized during the conjugate   #
                    # gradient process (check is this value decresae      #
                    # throughout the CG iterations )                      #
                    #-----------------------------------------------------#     
                    if(self.debug):
                        mgrad = np.zeros_like(grad)
                        mgrad[:] = -1.*grad[:]
                        self.norm_residual = normL2(self.residual)
                        alpha = (self.norm_residual**2) / dHd
                        self.qkm1_CG[:] = self.qk_CG[:]
                        grad_term = scalL2(self.descent,mgrad)
                        self.hessian_term = self.hessian_term+(alpha**2)*dHd
                        self.qk_CG = -grad_term + 0.5*self.hessian_term
            else:
                # if dHd > 0, : perform one conjugate gradient iteration
            
                # Update descent direction
                self.norm_residual = normL2(self.residual)
                alpha = (self.norm_residual**2)/dHd
                self.descent_prev[:] = self.descent[:]
                self.descent[:] = self.descent[:] + alpha*self.d[:]
                self.residual[:] = self.residual[:] + alpha*self.Hd[:]
            
                # Update CG direction
                norm_residual_prev = self.norm_residual
                self.norm_residual = normL2(self.residual)
                beta = (self.norm_residual**2)/(norm_residual_prev**2)
                self.d[:] = -1.*self.residual[:]+beta*self.d[:]
            
                # Update iteration counter 
                self.cpt_iter_CG = self.cpt_iter_CG+1
                #-----------------------------------------------------#
                # if the debug option is activated, compute the       #
                # quadratic function minimized during the conjugate   #
                # gradient process (check is this value decresae      #
                # throughout the CG iterations )                      #
                #-----------------------------------------------------#        
                if(self.debug):
                    mgrad = np.zeros_like(grad)
                    mgrad[:] = -1.*grad[:]
                    self.qkm1_CG = self.qk_CG
                    grad_term = scalL2(self.descent,mgrad)
                    descent_scal_Hd = scalL2(self.descent_prev,self.Hd)
                    self.hessian_term = self.hessian_term + (alpha**2) * dHd + 2. * alpha * descent_scal_Hd
                    self.qk_CG = -grad_term + 0.5*self.hessian_term
                
                # Check if the Eisenstat stopping critertion is satisfied
                self.conv_CG=( (self.norm_residual<=(self.eta*self.norm_grad)) or 
                    (self.cpt_iter_CG >= self.niter_max_CG))

                # Print information on the current CG iteration
                self.print_info_TRN(0.)


    def forcing_term_TRN(self, grad):
        ''' forcing_term_TRN
        '''
    
        # Computation of the forcing term self.eta following the formula
        eta_save = self.eta
        self.eisenvect[:] = grad[:] - self.residual[:]
        norm_eisenvect = normL2(self.eisenvect)
        self.eta = norm_eisenvect / self.norm_grad_m1
        
        # Additional safeguard if self.eta is too large
        eta_save_power = eta_save**((1.+np.sqrt(5.))/2.)
        if(eta_save_power > 0.1):
            self.eta = max(self.eta,eta_save_power)
        
        if(self.eta > 1.):
            self.eta = 0.9


    #--------------------------------------------------------------------#
    #         Preconditioned truncated Newton method: PTRN               # 
    #--------------------------------------------------------------------#
    def PTRN(self, x, fcost, grad, grad_preco):
        ''' Preconditioned truncated Newton method
        '''
    
        if(self.FLAG=='INIT'):
            # if FLAG is INIT, call the dedicated initialization def to allocate data structure self

            self.init_PTRN(x, fcost, grad)
            self.print_info_PTRN(fcost)
            self.comm = 'DES1'
            self.CG_phase = 'INIT'
            self.nfwd_pb = self.nfwd_pb+1
            self.conv_CG = False
            self.FLAG = 'NONE'
        
        if(self.comm == 'DES1'):
            # if self.comm is DES1, the selfizer starts the computation of a descent direction through the conjugate gradient
            if(self.CG_phase=='INIT'):
                # if self.CG_phase is INIT, initialization of the conjugate gradient process
                self.descent_PTRN0(grad, grad_preco)
                self.CG_phase = 'IRUN'
                self.comm = 'DES1'
                self.FLAG = 'HESS'
                self.nhess += 1
            elif(self.CG_phase == 'IRUN'):
                # if self.CG_phase is IRUN, iterate the conjugate gradient process (first part)
                
                self.descent_PTRN1(grad)
                if(self.conv_CG): 
                    #-----------------------------------------------------#
                    # if the conjugate gradient has already converged     # 
                    # (detection of a negative curvature), go to next     #
                    # phase: linesearch in the descent direction          #
                    #-----------------------------------------------------#
                    self.comm = 'NSTE'
                    self.CG_phase = 'INIT'        
                    self.FLAG = 'NONE'                   
                else:
                    # if the conjugate gradient has not converged prematurly,: ask the user to apply the preconditioner
                    self.FLAG = 'PREC'
                    self.comm = 'DES2'
                    
        elif(self.comm == 'DES2'):
            # if self.comm is DES2, the selfizer finish the current iteration of the conjugate gradient process
            self.descent_PTRN2(grad)

            if(self.conv_CG):
                # if the conjugate gradient has converged go to next phase: linesearch in the descent direction
                self.comm = 'NSTE'
                self.CG_phase = 'INIT'        
                self.FLAG = 'NONE'                
            else:
                # else start a new iteration of conjugate gradient and ask the user to compute a Hessian-vector product
                self.comm = 'DES1'
                self.FLAG = 'HESS'
                self.nhess = self.nhess + 1

        elif(self.comm == 'NSTE'):
            # if self.comm is NSTE, a descent direction has been computed, and a linesearch must be performed in this direction
            x = self.linesearch(x, fcost, grad)

            if(self.task=='NEW_STEP'):
                # if self.task is 'NEW_STEP, the linesearch process has found the new step
                self.cpt_iter += 1
                # Save the previous gradient norm 
                self.norm_grad_m1 = self.norm_grad
                # Computation of the new gradient norm 
                self.norm_grad = normL2(grad)
                # Print infor on current nonlinear iteration
                self.print_info_PTRN(fcost)        
                # Test for convergence
                if(self.std_test_conv(fcost)):
                    self.FLAG = 'CONV'
                    self.print_info_PTRN(fcost)
                else:
                    self.FLAG = 'NSTE'
                    # Flag for the computation of the new descent direction
                    self.comm = 'DES1'
                    # Update forcing term self.eta following the Eisenstat and Walker formula
                    self.forcing_term_TRN(grad)
                    
            elif(self.task == 'NEW_GRAD'): 
                # if self.task is 'NEW_GRAD, the linesearch process is continuing, the gradient at the current point is required
                self.FLAG = 'GRAD'         
                self.nfwd_pb += 1
            elif(self.task == 'FAILURE!'):        
                # if self.task is 'FAILURE, the linesearch process has failed, the iterations are stopped
                self.FLAG = 'FAIL'
                self.print_info_PTRN(fcost)

        return x

    def init_PTRN(self, x, fcost, grad):
        '''  Initilize linesearch parameters
        '''
        
        # set counters
        self.cpt_iter=0
        self.f0=fcost   
        self.nfwd_pb=0    
        self.nhess=0
        self.eta=0.9   
        #self.eta=0.1   
        #self.eta=1e-6   
            
        # initialize linesearch parameters
        self.fk=fcost
        self.cpt_iter_CG = 0
        
        # memory allocations
        self.xk = np.zeros_like(x)
        self.grad = np.zeros_like(x)
        self.descent = np.zeros_like(x)
        self.descent_prev = np.zeros_like(x)
        self.residual = np.zeros_like(x)
        self.residual_preco = np.zeros_like(x)
        self.d = np.zeros_like(x)
        self.Hd = np.zeros_like(x)
        self.eisenvect = np.zeros_like(x)

        self.xk[:] = x[:]
        self.grad[:] = grad[:]

        # norm of the first gradient
        self.norm_grad = normL2(grad)


    def descent_PTRN0(self, grad, grad_preco):
        ''' Descent direction
        '''
    
        # Initialization of the conjugate gradient process
        self.residual[:] = grad[:]
        self.residual_preco[:] = grad_preco[:]
        self.d[:] = -1. * self.residual_preco[:]
        self.Hd[:] = 0.
        self.descent[:] = 0.
        self.qk_CG = 0.
        self.hessian_term = 0.
        self.norm_residual = normL2(self.residual)
        self.conv_CG = False
        self.cpt_iter_CG = 0
        self.print_info_PTRN(0.0)

    def descent_PTRN1(self, grad):
        ''' descent_PTRN1
        '''
    
        # start one conjugate gradient iteration
        self.dHd = scalL2(self.d, self.Hd)
        if(self.dHd < 0.):
            # if dHd < 0, detection of a negative eigenvalue of the Hessian operator, stop the process
            self.conv_CG = True
            print('Negative curvature')
            if(self.cpt_iter_CG == 0):
                #-----------------------------------------------------#
                # if this is the first iteration,:return the          #
                # opposite of the preconditioned gradient as descent  #
                # direction: preconditioned steepest descent direction#
                #-----------------------------------------------------#   
                self.descent[:] = self.d[:]
                #-----------------------------------------------------#
                # if the debug option is activated, compute the       #
                # quadratic function minimized during the conjugate   #
                # gradient process (check is this value decresae      #
                # throughout the CG iterations )                      #
                #-----------------------------------------------------#
                if(self.debug):
                    self.res_scal_respreco = scalL2(self.residual,self.residual_preco)
                    self.alpha_CG = self.res_scal_respreco / self.dHd
                    self.qkm1_CG = self.qk_CG
                    mgrad = np.zeros_like(grad)
                    mgrad[:] = -1. * grad[:]
                    grad_term = scalL2(self.descent,mgrad)
                    self.hessian_term = self.hessian_term + (self.alpha_CG**2) * self.dHd
                    self.qk_CG = -grad_term + 0.5 * self.hessian_term
        else:
            # if dHd > 0,: start one conjugate gradient iteration

            # Update descent direction
            self.res_scal_respreco = scalL2(self.residual,self.residual_preco)
            self.alpha_CG = self.res_scal_respreco / self.dHd
            self.descent_prev[:] = self.descent[:]
            self.descent[:] = self.descent[:] + self.alpha_CG * self.d[:]
            self.residual[:] = self.residual[:] + self.alpha_CG * self.Hd[:]
            
            # STOP HERE and wait for preconditioning
        

    def descent_PTRN2(self, grad):
        ''' descent_PTRN2
        '''
        
        # continue the current conjugate gradient iteration
        # Update CG direction
        res_scal_respreco_prev = self.res_scal_respreco
        self.res_scal_respreco = scalL2(self.residual, self.residual_preco)
        beta = (self.res_scal_respreco) / (res_scal_respreco_prev)
        self.d[:] = -1. * self.residual_preco[:] + beta * self.d[:]
        
        # Update iteration counter 
        self.cpt_iter_CG += 1

        #-----------------------------------------------------#
        # if the debug option is activated, compute the       #
        # quadratic function minimized during the conjugate   #
        # gradient process (check is this value decresae      #
        # throughout the CG iterations )                      #
        #-----------------------------------------------------#  
        if(self.debug):
            self.qkm1_CG = self.qk_CG
            mgrad = np.zeros_like(grad)
            mgrad[:]=-1.*grad[:]
            grad_term = scalL2(self.descent,mgrad)
            descent_scal_Hd = scalL2(self.descent_prev,self.Hd)
            self.hessian_term = self.hessian_term+(self.alpha_CG**2) * self.dHd + 2. * self.alpha_CG * descent_scal_Hd
            self.qk_CG = -grad_term + 0.5 * self.hessian_term

        # Check if the Eisenstat stopping critertion is satisfied
        self.norm_residual = normL2(self.residual)
        self.conv_CG=((self.norm_residual<=(self.eta*self.norm_grad)) or 
                    (self.cpt_iter_CG>=self.niter_max_CG))
        
        # Print information on the current CG iteration
        self.print_info_PTRN(0.)


    #--------------------------------------------------------------------#
    #                               linesearch                           # 
    #--------------------------------------------------------------------#
    def linesearch(self, x, fcost, grad):
        ''' Linesearch
        '''
        if(self.first_ls):
            # FIRST LINESEARCH: initialization step
            self.fk = fcost
            self.q0 = scalL2(grad, self.descent)

            # Set the search interval bounds to 0
            self.alpha_L = 0.
            self.alpha_R = 0.
            self.task = 'NEW_GRAD'
            self.first_ls = False
            self.xk[:] = x[:]
            x[:] = self.xk[:] + self.alpha * self.descent[:]

            # If bounds activated, project x into the feasible range
            if(self.bound == 1):
                x = self.project(x)
            self.cpt_ls=0

        elif (self.cpt_ls >= self.nls_max and fcost < self.fk):
            # if the number of linesearch iteration outreaches the maximum allowed
            # but a decrease of the misfit is produced then accept the steplength

            self.task = 'NEW_STEP'
            self.first_ls = True
            
            # Compute new x in the descent direction
            x[:] = self.xk[:] + self.alpha * self.descent[:]
            
            # If bounds activated, project x into the feasible range
            if(self.bound == 1):
                x = self.project(x)

        elif(self.cpt_ls >= self.nls_max):
            # if the number of linesearch iteration outreaches the maximum allowed
            # without decreasing the misfit then the linesearch has failed
            self.task='FAILURE!'
        else:
            # If not initialization step and number of linesearch iteration ok
            # then perform one linesearch iteration                     
            
            self.q = scalL2(grad, self.descent)

            if (fcost <= (self.fk+(self.m1*self.alpha*self.q0))) and (self.q >= (self.m2*self.q0)):
                #--------------------------------------------------------------------#
                # First test if the Wolfe conditions are satisfied with              #     
                # current steplength, if this is the case, linesearch                # 
                # ends here with success                                             #
                #--------------------------------------------------------------------#
                self.task = 'NEW_STEP'
                self.first_ls = True
                if(self.debug):
                    print('fcost :', fcost)
                    print('optimize.f0    :', self.f0)
                    print('optimize.fk    :', self.fk)
                    print('optimize.alpha :', self.alpha)
                    print('optimize.q     :', self.q)
                    print('optimize.q0    :', self.q0)
                    print('m1             :', self.m1)
                    print('cpt_ls is      : ', self.cpt_ls)

            elif (fcost > (self.fk + (self.m1 * self.alpha * self.q0))):
                #--------------------------------------------------------------------#
                # If the first condition is not satisfied then shrink the            #
                # search interval                                                    #
                #--------------------------------------------------------------------#
                if(self.debug):
                    print('failure 1')
                    print('fcost          :', fcost)
                    print('optimize.fk    :', self.fk)
                    print('optimize.alpha :', self.alpha)
                    print('optimize.q0    :', self.q0)
                    print('m1             :', self.m1)
                    print('cpt_ls is      : ', self.cpt_ls)
                
                self.alpha_R = self.alpha
                new_alpha = (self.alpha_L+self.alpha_R) / 2.
                self.alpha = new_alpha
                self.task = 'NEW_GRAD'
                self.cpt_ls = self.cpt_ls+1

            elif( (fcost <= (self.fk + (self.m1 * self.alpha * self.q0))) and (self.q < (self.m2 * self.q0))):
                #--------------------------------------------------------------------#
                # If the second condition is not satisfied then shrink the           #
                # search interval unless the right bound of the search interval      #
                # as not yet been defined                                            #
                #--------------------------------------------------------------------#
                if(self.debug):
                    print('failure 2')
                    print('fcost          :', fcost)
                    print('optimize.fk    :', self.fk)
                    print('optimize.alpha :', self.alpha)
                    print('optimize.q0    :', self.q0)
                    print('optimize.q     :', self.q)
                    print('m1             :', self.m1)
                    print('m2             :', self.m2)
                    print('cpt_ls is      : ', self.cpt_ls)

                self.alpha_L = self.alpha

                if(self.alpha_R != 0.):
                    new_alpha = (self.alpha_L + self.alpha_R) / 2.
                else:
                    new_alpha = self.mult_factor * self.alpha
                
                self.alpha = new_alpha
                self.task = 'NEW_GRAD'
                self.cpt_ls = self.cpt_ls + 1
            
            # Compute new x in the descent direction
            x[:] = self.xk[:] + self.alpha * self.descent[:]

            # If bounds activated, project x into the feasible range
            if(self.bound == 1):
                x = self.project(x)

        return x


    def project(self, x):
        '''
        project the vector x into the box defined by the bound constraints optim.lb and optim.ub 
        '''

        n = np.size(x)
        for i in range(n):
            if(x[i] > self.ub[i]):
                x[i] = self.ub[i] - self.threshold

            if(x[i] < self.lb[i]):
                x[i] = self.lb[i] + self.threshold

        return x

     
    def std_test_conv(self, fcost):
        '''
        This routine implements a simple convergence test based on the relative cost decrease. 
        If the current relative cost is lower than a value set by the user in optim.conv
        then test_conv is returned equal to True, otherwise, it is returned equal to False
        '''
    
        test_conv = ((fcost/self.f0 < self.conv) or (self.cpt_iter >= self.niter_max))

        return test_conv


    def print_info(self, fcost):
        ''' print  information
        '''

        ng = normL2(self.grad)

        if(self.FLAG == 'INIT'):
            f = open('iterate_%s.log'%(self.method), 'w')
            if(self.method == 'SD'):
                f.write('**********************************************************************\n')
                f.write('         STEEEPEST DESCENT ALGORITHM         \n')
                f.write('**********************************************************************\n')
            elif(self.method == 'CG'):
                f.write('**********************************************************************\n')
                f.write('         NONLINEAR CONJUGATE GRADIENT ALGORITHM         \n')
                f.write('**********************************************************************\n')
            elif(self.method == 'LBFGS'):
                f.write('**********************************************************************\n')
                f.write('             l-BFGS ALGORITHM                \n')
                f.write('**********************************************************************\n')
            elif(self.method == 'PLBFGS'):
                f.write('**********************************************************************\n')
                f.write('             PRECONDITIONED l-BFGS ALGORITHM                \n')
                f.write('**********************************************************************\n')

            f.write('     Convergence criterion  : %10.2e\n' %self.conv)
            f.write('     Niter_max              : %10d\n'   %self.niter_max)
            f.write('     Initial cost is        : %10.2e\n' %self.f0)
            f.write('     Initial norm_grad is   : %10.2e\n' %ng)
            f.write('**********************************************************************\n')
            f.write('   Niter      fk         ||gk||       fk/f0        alpha        nls      ngrad    \n')
            f.write('%6d %12.2e %12.2e %12.2e %12.2e %8d %8d\n'
            %(self.cpt_iter, fcost, ng, fcost/self.f0, self.alpha, self.cpt_ls, self.nfwd_pb))
            f.close()

        elif(self.FLAG == 'CONV'):
            f = open('iterate_%s.log'%(self.method), 'a')
            f.write('%6d %12.2e %12.2e %12.2e %12.2e %8d %8d\n'
            %(self.cpt_iter, fcost, ng, fcost/self.f0, self.alpha, self.cpt_ls, self.nfwd_pb))
            f.write('**********************************************************************\n')

            if(self.cpt_iter >= self.niter_max):
                f.write('  STOP: MAXIMUM NUMBER OF ITERATION REACHED    \n')
            else:
                f.write('  STOP: CONVERGENCE CRITERION SATISFIED        \n')
                f.write('**********************************************************************\n')
                f.close()

        elif(self.FLAG == 'FAIL'):
            f = open('iterate_%s.log'%(self.method), 'a')
            f.write('**********************************************************************\n')
            f.write('  STOP: LINESEARCH FAILURE    \n')
            f.write('**********************************************************************\n')
            f.close()
        else:
            f = open('iterate_%s.log'%(self.method), 'a')
            f.write('%6d %12.2e %12.2e %12.2e %12.2e %8d %8d\n'
            %(self.cpt_iter, fcost, ng, fcost/self.f0, self.alpha, self.cpt_ls, self.nfwd_pb))
            f.close()


    def print_info_TRN(self, fcost):
        
        if(self.FLAG=='INIT') :
            f = open('iterate_TRN.log', 'w')
            f.write('******************************************************************************************\n')
            f.write('                                 TRUNCATED NEWTON ALGORITHM                               \n')
            f.write('******************************************************************************************\n')
            f.write('     Convergence criterion  : %10.2e \n' %self.conv)
            f.write('     Niter_max              : %7d \n'    %self.niter_max)
            f.write('     Initial cost is        : %10.2e \n' %self.f0)
            f.write('     Initial norm_grad is   : %10.2e \n' %self.norm_grad)
            f.write('     Maximum CG iter        : %7d \n'    %self.niter_max_CG)
            f.write('******************************************************************************************\n')
            f.write('   Niter      fk         ||gk||       fk/f0          alpha      nls     nit_CG      eta       ngrad     nhess \n')
            f.write('%6d %12.2e %12.2e %12.2e %12.2e %8d %8d %12.2e %8d %8d \n'%(
            self.cpt_iter, fcost, self.norm_grad, fcost/self.f0, self.alpha, self.cpt_ls, self.cpt_iter_CG, self.eta, self.nfwd_pb, self.nhess))
            f.close()

            f = open('iterate_TRN_CG.log', 'w')
            f.write('******************************************************************************************\n')
            f.write('                                 TRUNCATED NEWTON ALGORITHM                               \n')
            f.write('                                      INNER CG HISTORY                                    \n')
            f.write('******************************************************************************************\n')
            f.write('     Convergence criterion  : %10.2e \n' %self.conv)
            f.write('     Niter_max              : %7d \n'    %self.niter_max)
            f.write('     Initial cost is        : %10.2e \n' %self.f0)
            f.write('     Initial norm_grad is   : %10.2e \n' %self.norm_grad)
            f.write('     Maximum CG iter        : %7d \n'    %self.niter_max_CG)
            f.write('******************************************************************************************\n')
            f.close()
        elif(self.FLAG=='CONV'):
            f = open('iterate_TRN.log', 'a')
            f.write('**********************************************************************\n')
            if(self.cpt_iter == self.niter_max):
                f.write('  STOP: MAXIMUM NUMBER OF ITERATION REACHED    \n')
            else:
                f.write('  STOP: CONVERGENCE CRITERION SATISFIED        \n')
            
            f.write('**********************************************************************\n')
            f.close()

        elif(self.FLAG=='FAIL'):
            f = open('iterate_TRN.log', 'a')
            f.write('**********************************************************************\n')
            f.write('  STOP: LINESEARCH FAILURE    \n')
            f.write('**********************************************************************\n')
            f.close()
        elif(self.comm=='DESC'):
            f = open('iterate_TRN_CG.log', 'a')
            if(self.CG_phase=='INIT'):
                f.write('-------------------------------------------------------------------------------------------------\n')
                f.write(' NONLINEAR ITERATION %4d  ETA IS : %12.2e \n'%(self.cpt_iter, self.eta))
                f.write('-------------------------------------------------------------------------------------------------\n')
                f.write('  Iter_CG        qk       norm_res      norm_res/||gk||  \n')
                f.write('%8d %12.2e %12.2e %12.2e \n'%(self.cpt_iter_CG, self.qk_CG, self.norm_residual, self.norm_residual/self.norm_grad))
            else:
                f.write('%8d %12.2e %12.2e %12.2e \n'%(self.cpt_iter_CG, self.qk_CG, self.norm_residual, self.norm_residual/self.norm_grad))
            f.close()            
        else:
            f = open('iterate_TRN.log', 'a')
            f.write('%6d %12.2e %12.2e %12.2e %12.2e %8d %8d %12.2e %8d %8d \n'%(
            self.cpt_iter, fcost, self.norm_grad, fcost/self.f0, self.alpha, self.cpt_ls, self.cpt_iter_CG, self.eta, self.nfwd_pb, self.nhess))
            f.close()


    def print_info_PTRN(optim, fcost, FLAG):
        ''' Print information
        '''

        if(FLAG=='INIT'):
            f = open('iterate_PTRN.log', 'w')
            f.write('******************************************************************************************\n')
            f.write('                      PRECONDITIONED TRUNCATED NEWTON ALGORITHM                           \n')
            f.write('******************************************************************************************\n')
            f.write('     Convergence criterion  : %10.2e \n'%optim.conv)
            f.write('     Niter_max              : %7d \n'%optim.niter_max)
            f.write('     Initial cost is        : %10.2e \n'%optim.f0)
            f.write('     Initial norm_grad is   : %10.2e \n'%optim.norm_grad)
            f.write('     Maximum CG iter        : %7d \n'%optim.niter_max_CG)
            f.write('******************************************************************************************\n')
            f.write('   Niter      fk         ||gk||       fk/f0          alpha      nls     nit_CG      eta       ngrad     nhess \n')
            f.write('%6d %12.2e %12.2e %12.2e %12.2e %8d %8d %12.2e %8d %8d \n'%(
            optim.cpt_iter, fcost, optim.norm_grad, fcost/optim.f0, optim.alpha, optim.cpt_ls, optim.cpt_iter_CG, optim.eta, optim.nfwd_pb, optim.nhess))
            f.close()

            f = open('iterate_PTRN_CG.log', 'w')
            f.write('******************************************************************************************\n')
            f.write('                        PRECONDITIONED TRUNCATED NEWTON ALGORITHM                         \n')
            f.write('                                      INNER CG HISTORY                                    \n')
            f.write('******************************************************************************************\n')
            f.write('     Convergence criterion  : %10.2e \n'%optim.conv)
            f.write('     Niter_max              : %7d \n'%optim.niter_max)
            f.write('     Initial cost is        : %10.2e \n'%optim.f0)
            f.write('     Initial norm_grad is   : %10.2e \n'%optim.norm_grad)
            f.write('     Maximum CG iter        : %7d \n'%optim.niter_max_CG)
            f.write('******************************************************************************************\n')
            f.close()

        elif(FLAG=='CONV'):
            f = open('iterate_PTRN.log', 'a')
            f.write('**********************************************************************\n')
            if(optim.cpt_iter==optim.niter_max):
                f.write('  STOP: MAXIMUM NUMBER OF ITERATION REACHED    \n')
            else:
                f.write('  STOP: CONVERGENCE CRITERION SATISFIED        \n')
            f.write('**********************************************************************\n')
            f.close()

        elif(FLAG=='FAIL'):
            f = open('iterate_PTRN.log', 'a')
            f.write('**********************************************************************\n')
            f.write('  STOP: LINESEARCH FAILURE    \n')
            f.write('**********************************************************************\n')
            f.close()     
        elif(optim.comm=='DES1'):
            f = open('iterate_PTRN_CG.log', 'a')
            f.write('-------------------------------------------------------------------------------------------------\n')
            f.write(' NONLINEAR ITERATION %4d  ETA IS : %12.2e \n'%(optim.cpt_iter, optim.eta))
            f.write('-------------------------------------------------------------------------------------------------\n')
            f.write('  Iter_CG        qk       norm_res      norm_res/||gk||  \n')
            f.write('%8d %12.2e %12.2e %12.2e \n'%(optim.cpt_iter_CG, optim.qk_CG, optim.norm_residual, optim.norm_residual/optim.norm_grad))
            f.close()     

        elif(optim.comm=='DES2'):
            f = open('iterate_PTRN_CG.log', 'a')
            f.write('%8d %12.2e %12.2e %12.2e \n'%(optim.cpt_iter_CG, optim.qk_CG, optim.norm_residual, optim.norm_residual/optim.norm_grad))
            f.close()     
        else:
            f = open('iterate_PTRN.log', 'a')
            f.write('%6d %12.2e %12.2e %12.2e %12.2e %8d %8d %12.2e %8d %8d \n'%(
            optim.cpt_iter, fcost, optim.norm_grad, fcost/optim.f0, optim.alpha, optim.cpt_ls, optim.cpt_iter_CG, optim.eta, optim.nfwd_pb, optim.nhess))
            f.close()


    def print_info_PTRN(self, fcost):
        ''' Print information
        '''
        
        if(self.FLAG=='INIT'):
            f = open('iterate_PTRN.log', 'w')
            f.write('******************************************************************************************\n')
            f.write('                      PRECONDITIONED TRUNCATED NEWTON ALGORITHM                           \n')
            f.write('******************************************************************************************\n')
            f.write('     Convergence criterion  : %10.2e \n'%self.conv)
            f.write('     Niter_max              : %7d \n'   %self.niter_max)
            f.write('     Initial cost is        : %10.2e \n'%self.f0)
            f.write('     Initial norm_grad is   : %10.2e \n'%self.norm_grad)
            f.write('     Maximum CG iter        : %7d \n'   %self.niter_max_CG)
            f.write('******************************************************************************************\n')
            f.write('   Niter      fk         ||gk||       fk/f0          alpha      nls     nit_CG      eta       ngrad     nhess \n')
            f.write('%6d %12.2e %12.2e %12.2e %12.2e %8d %8d %12.2e %8d %8d \n'%(
            self.cpt_iter, fcost, self.norm_grad, fcost/self.f0, self.alpha, self.cpt_ls, self.cpt_iter_CG, self.eta, self.nfwd_pb, self.nhess))
            f.close()

            f = open('iterate_PTRN_CG.log', 'w')
            f.write('******************************************************************************************\n')
            f.write('                        PRECONDITIONED TRUNCATED NEWTON ALGORITHM                         \n')
            f.write('                                      INNER CG HISTORY                                    \n')
            f.write('******************************************************************************************\n')
            f.write('     Convergence criterion  : %10.2e \n'%self.conv)
            f.write('     Niter_max              : %7d \n'   %self.niter_max)
            f.write('     Initial cost is        : %10.2e \n'%self.f0)
            f.write('     Initial norm_grad is   : %10.2e \n'%self.norm_grad)
            f.write('     Maximum CG iter        : %7d \n'   %self.niter_max_CG)
            f.write('******************************************************************************************\n')
            f.close()

        elif(self.FLAG=='CONV'):
            f = open('iterate_PTRN.log', 'a')
            f.write('**********************************************************************\n')
            if(self.cpt_iter == self.niter_max):
                f.write('  STOP: MAXIMUM NUMBER OF ITERATION REACHED    \n')
            else:
                f.write('  STOP: CONVERGENCE CRITERION SATISFIED        \n')
            f.write('**********************************************************************\n')
            f.close()

        elif(self.FLAG=='FAIL'):
            f = open('iterate_PTRN.log', 'a')
            f.write('**********************************************************************\n')
            f.write('  STOP: LINESEARCH FAILURE    \n')
            f.write('**********************************************************************\n')
            f.close()     
        elif(self.comm=='DES1'):
            f = open('iterate_PTRN_CG.log', 'a')
            f.write('-------------------------------------------------------------------------------------------------\n')
            f.write(' NONLINEAR ITERATION %4d  ETA IS : %12.2e \n'%(self.cpt_iter, self.eta))
            f.write('-------------------------------------------------------------------------------------------------\n')
            f.write('  Iter_CG        qk       norm_res      norm_res/||gk||  \n')
            f.write('%8d %12.2e %12.2e %12.2e \n'%(self.cpt_iter_CG, self.qk_CG, self.norm_residual, self.norm_residual/self.norm_grad))
            f.close()     

        elif(self.comm=='DES2'):
            f = open('iterate_PTRN_CG.log', 'a')
            f.write('%8d %12.2e %12.2e %12.2e \n'%(self.cpt_iter_CG, self.qk_CG, self.norm_residual, self.norm_residual/self.norm_grad))
            f.close()     
        else:
            f = open('iterate_PTRN.log', 'a')
            f.write('%6d %12.2e %12.2e %12.2e %12.2e %8d %8d %12.2e %8d %8d \n'%(
            self.cpt_iter, fcost, self.norm_grad, fcost/self.f0, self.alpha, self.cpt_ls, self.cpt_iter_CG, self.eta, self.nfwd_pb, self.nhess))
            f.close()

# =================================================================================================
# utility functions
# =================================================================================================

def normL2(x):
    '''The routine normL2 returns the Euclidian norm of a vector x of size n
    
    Args:

    Returns:
    '''

    return np.linalg.norm(x)


def scalL2(x, y):
    '''
    The routine scalL2 returns the scalar product between two vectors x and y
    '''

    return np.inner(x, y)