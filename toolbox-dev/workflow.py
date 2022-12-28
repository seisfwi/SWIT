###############################################################################
#
# SWIT v1.1: Seismic Waveform Inversion Toolbox
#
#   A Python package for seismic waveform inversion
#   Developed by Haipeng Li at USTC, updated on 2022-12-21 at Stanford
#   haipengl@mail.ustc.edu.cn, haipeng@stanford.edu
#
#   Workflow module, defining workflow of FWI and RTM
#
###############################################################################

from multiprocessing import Pool
import os
import numpy as np
from misfit import calculate_adjoint_misfit_is
from tools import load_float


class FWI(object):
    ''' Full waveform inversion workflow
    '''

    def __init__(self, solver, optimizer, processor = None):
        ''' Initialize FWI workflow
        '''

        self.solver = solver
        self.optimizer = optimizer
        self.processor = processor


    def prepare_adjoint_misfit(self):
        ''' Prepare the adjoint source and calculate misfit function using the multiprocessing

        Returns
        -------
        rsd : list of float
            each element is the summed misfit over all traces for one source
        
        Note: the adjoint source will be saved in workpath/config/wavelet/src1_adj.bin ...
        '''

        # get the parameters
        work_path = self.solver.config.path
        pool_num  = self.solver.config.mpi_num
        src_num  = self.solver.source.num
        dt = self.solver.model.dt
        nt = self.solver.model.nt
        misfit_type = self.optimizer.misfit_type

        # initialize the Pool object
        pool = Pool(pool_num)

        # apply the function to all sources
        results = [pool.apply_async(calculate_adjoint_misfit_is, 
        (isrc, work_path, dt, nt, misfit_type, )) for isrc in range(src_num)]
        
        # close the pool and wait for the work to finish
        pool.close()

        # get the misfits for all sources (not the summed misfit)
        fcost_all_src = [p.get() for p in results]

        # join the processes
        pool.join()
        
        # sum the misfit over all sources
        fcost = 0.0
        for _fcost in fcost_all_src:
            fcost += _fcost

        # return the summed misfit over all sources
        return fcost_all_src, fcost



    def calculate_gradient_misfit(self, vp = None, rho = None):
        ''' Calculate gradient of misfit function

        Parameters
        ----------
        vp : flatted 1D array (float32)
            velocity model
        rho : flatted 1D array (float32)
            density model (optional)
        '''
        
        if vp is None and rho is None:
            msg = 'ERROR: no model provided for gradient calculation'
            raise ValueError(msg)

        # set model
        self.solver.set_model(vp = vp, rho = rho)

        # generate sythetic data (save boundary data for wavefield reconstruction)
        self.solver.run(simu_type = 'forward', simu_tag = 'syn', save_boundary = True)

        # process synthetic data
        # self.processor.run()
        
        # prepare adjoint source and calculate misfit
        fcost_all_src, fcost = self.prepare_adjoint_misfit()

        # calculate gradient via adjoint modeling
        self.solver.run(simu_type = 'gradient', simu_tag = 'syn')

        # process gradient
        grad = self.postprocess_gradient()

        # save misfit and gradient for each iteration

        return fcost_all_src, fcost, grad



    def postprocess_gradient(self):
        ''' Post-process gradient

            1. read gradient from binary files and sum over all sources
            2. apply mask 
            3. apply smoothness
            4. apply normalization with proper scaling factor
        '''
        # get the parameters
        nx = self.solver.model.nx
        nz = self.solver.model.nz 
        src_num  = self.solver.source.num

        # read gradient
        grad = np.zeros((nx, nz))
        for isrc in range(src_num):
            path = os.path.join(self.solver.config.path, 'data/syn/src{}/vp_gradient.bin'.format(isrc+1))
            grad += load_float(path).reshape(nx, nz)
        grad = grad/grad.max() * self.optimizer.max_update_val

        return grad


    def run(self):
        ''' Run FWI workflow
        '''
        
        # calculate gradient and misfit from initial model
        vp = self.optimizer.vp_init
        rho = self.optimizer.rho_init
        fcost_all_src, fcost, grad = self.calculate_gradient_misfit(vp = vp, rho = rho)
        grad_preco = np.copy(grad)

        count = 0
        # save initial misfit and gradient
        # self.save_results(count, vp, grad)

        print("Initial fcost: {} \n".format(fcost))

        # keep iterate while convergence not reached or linesearch not failed
        while ((self.optimizer.FLAG != 'CONV') and (self.optimizer.FLAG != 'FAIL')):
            
            # flatten 1D array
            vp = vp.flatten()
            grad = grad.flatten()
            grad_preco= grad_preco.flatten()

            # update model and linesearch using the preconditioned gradient
            vp = self.optimizer.iterate(vp, fcost, grad, grad_preco)
 
            if(self.optimizer.FLAG == 'GRAD'):
                
                # compute cost and gradient of the updated model
                count += 1
                fcost_all_src, fcost, grad = self.calculate_gradient_misfit(vp = vp, rho = rho)
                
                # precondition gradient
                grad_preco = np.copy(grad)

                # save iteration history
                # self.save_results(count, vp, grad, fcost_all)

        print('END OF FWI')
        print('See the convergence history in iterate_{}.dat'.format(self.optimizer.scheme))
        


    # def save_results(self, count, vp, grad, fcost_all):
    #     ''' Save results during FWI iterations
    #     '''
        
    #     np.save(self.solver.path + 'outputs/grad_ite_{:04d}.bin'.format(count), grad.flatten())
    #     np.save(self.solver.path + 'outputs/vp_ite_{:04d}.bin'.format(count), vp.flatten())
    #     np.save(self.solver.path + 'outputs/misfit_iter_{:04d}.bin'.format(count), fcost_all)



class RTM(object):
    ''' Reverse time migration (RTM) class
    '''

    def __init__(self, solver, processor):
        ''' Initialize RTM class

        Parameters
        ----------
        solver : Solver object
            solver object
        processor : Processor object
            processor object
        '''
        self.solver = solver
        self.processor = processor
    
    
    def run(self):
        ''' Run RTM workflow
        '''
        
        raise NotImplementedError