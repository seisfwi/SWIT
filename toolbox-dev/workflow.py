###############################################################################
#
# SWIT v1.1: Seismic Waveform Inversion Toolbox
#
#   A Python package for seismic waveform inversion
#   By Haipeng Li at USTC & Stanford
#   Email: haipengl@mail.ustc.edu.cn, haipeng@stanford.edu 
#
#   Workflow mudule defines the implementation of FWI and RTM workflows
#
###############################################################################

import os
import time
from multiprocessing import Pool

import numpy as np
from misfit import calculate_adjoint_misfit_is
from plot import plot_misfit, plot_model
from tools import smooth2d
from utils import preconditioner


class FWI(object):
    ''' Full Waveform Inversion workflow

    Parameters
    ----------
        solver : object
            solver object
        optimizer : object
            optimizer object
        preprocessor : object
            preprocessor object
    '''

    def __init__(self, solver, optimizer, preprocessor):
        ''' Initialize FWI workflow
        '''
        # basic parameters
        self.solver = solver
        self.optimizer = optimizer
        self.preprocessor = preprocessor

        # build directories
        self.__build_dir()

        # check the existence of obs data
        self.__check_obs_data()
        

    def run(self):
        ''' Run FWI workflow
        '''

        # start the timer
        start_time = time.time()

        # preprocess the obs data
        self.preprocessor.run(data_path = self.solver.system.path + 'data/obs/', 
                            src_num = self.solver.source.num, 
                            mpi_num = self.solver.system.mpi_num, 
                            nt = self.solver.model.nt, 
                            dt = self.solver.model.dt,  
                            src_coord = self.solver.source.coord, 
                            rec_coord = self.solver.receiver.coord)

        # compute cost and gradient for the initial model
        vp  = self.optimizer.vp_init
        rho = self.optimizer.rho_init # rho is not updated in FWI workflow
        fcost, fcost_all, grad_preco = self.objective_function(vp = vp, rho = rho)

        # iterate until convergence or linesearch failure
        while ((self.optimizer.FLAG != 'CONV') and (self.optimizer.FLAG != 'FAIL')):

            # update model and linesearch using the preconditioned gradient
            vp_pre = np.copy(vp)
            vp = self.optimizer.iterate(vp, fcost, grad_preco, np.copy(grad_preco))

            # calculate gradient and misfit from updated model
            if(self.optimizer.FLAG == 'GRAD'):
                # print information
                self.print_iterate_info(fcost, vp, vp_pre)

                # compute cost and gradient for the updated model
                fcost, fcost_all, grad_preco = self.objective_function(vp = vp, rho = rho)
                
                # save the iteration history
                self.save_results(vp, grad_preco, fcost_all)

                # plot the model, gradient and misfit
                self.plot_results(vp, grad_preco)

        # print the end information
        self.print_end_info(start_time)


    def objective_function(self, vp = None, rho = None):
        ''' Objective function for FWI

        Parameters
        ----------
        vp : flatted 1D array (float32)
            velocity model
        rho : flatted 1D array (float32)
            density model (optional)
        '''
        
        if vp is None and rho is None:
            msg = 'FWI workflow ERROR: no model provided for gradient calculation'
            raise ValueError(msg)

        # reset model to solver
        self.solver.set_model(vp = vp, rho = rho)

        # model the syn data (save boundary for wavefield reconstruction)
        self.solver.run(simu_type = 'forward', simu_tag = 'syn', save_boundary = True)

        # preprocess the syn data
        self.preprocessor.run(data_path = self.solver.system.path + 'data/syn/', 
                            src_num = self.solver.source.num, 
                            mpi_num = self.solver.system.mpi_num, 
                            nt = self.solver.model.nt, 
                            dt = self.solver.model.dt,  
                            src_coord = self.solver.source.coord, 
                            rec_coord = self.solver.receiver.coord)
        
        # calculate adjoint source and calculate misfit
        fcost, fcost_all = self.calculate_adjoint_misfit()

        # calculate the gradient via adjoint modeling
        grad, for_illum, adj_illum = self.solver.run(simu_type = 'gradient', simu_tag = 'syn')

        # postprocess the gradient
        grad_preco = self.postprocess_gradient(grad, for_illum, adj_illum)

        # save misfit and gradient for each iteration
        return fcost, fcost_all, grad_preco


    def calculate_adjoint_misfit(self):
        ''' Prepare the adjoint source and calculate misfit function

        Returns
        -------
            rsd : list of float
                each element is the summed misfit over all traces for one source
            
            Note: the adjoint source is saved in config/wavelet/src1_adj.bin, ...
        '''

        # initialize the Pool object
        pool = Pool(self.solver.system.mpi_num)

        # TODO: add the MPI support for the calculation of misfit
        # apply to all sources
        results = [pool.apply_async(calculate_adjoint_misfit_is, (isrc, 
                            self.solver.system.path, 
                            self.solver.model.dt, 
                            self.solver.model.nt, 
                            self.optimizer.misfit_type, ))
                for isrc in range(self.solver.source.num)]

        # close the pool and wait for the work to finish
        pool.close()

        # get the misfits for all sources (not the summed misfit)
        fcost_all = np.array([p.get() for p in results])

        # block at this line until all processes are done
        pool.join()
    
        # return the summed misfit over all sources
        return np.sum(fcost_all), fcost_all


    def postprocess_gradient(self, grad, for_illum, adj_illum):
        ''' Postprocess gradient. 
            Note that the sequence of operations is important here.

        Parameters
        ----------
            grad : 2D array
                gradient of misfit function
            for_illum : 2D array
                forward illumination
            adj_illum : 2D array
                adjoint illumination

        Returns
        -------
        grad : 2D array
            preconditioned gradient of misfit function
        '''

        # mask the gradient
        grad = grad * self.optimizer.grad_mask
   
        # preconditioning the gradient
        grad = grad / np.power(preconditioner(for_illum, adj_illum), 1.)

        # smooth the gradient if needed
        if self.optimizer.grad_smooth_size > 0:
            grad = smooth2d(grad, span = self.optimizer.grad_smooth_size)

        # scale the gradient with the maximum allowed update value
        grad = grad / abs(grad).max() * self.optimizer.update_vpmax

        return grad


    def print_iterate_info(self, fcost, vp, vp_pre):
        ''' Print information during FWI iterations
        '''
        # calculate the maximum update
        max_update = np.max(np.abs(vp-vp_pre.flatten()))

        # print the information
        print('Iteration: {} \t fcost: {:.4e} \t step: {:7.4f} \t max update: {:7.2f} m/s'.format(
                self.optimizer.cpt_iter+1, fcost, self.optimizer.alpha, max_update))
    
    
    def print_end_info(self, start_time):
        ''' Print end information
        '''
        # print the end information
        hours, rem = divmod(time.time()-start_time, 3600)
        minutes, seconds = divmod(rem, 60)
        print('\nFWI workflow: finished in {:0>2}h {:0>2}m {:.0f}s'.format(int(hours), int(minutes), seconds))
        print('FWI workflow: see convergence history in iterate_{}.log\n'.format(self.optimizer.method))


    def save_results(self, vp, grad, fcost_all):
        ''' Save results during FWI iterations
        '''

        iter = self.optimizer.cpt_iter+1
        np.save(self.solver.system.path + 'fwi/grad/grad_it_{:04d}.npy'.format(iter), grad)
        np.save(self.solver.system.path + 'fwi/model/vp_it_{:04d}.npy'.format(iter), vp)
        np.save(self.solver.system.path + 'fwi/misfit/fcost_all_it_{:04d}.npy'.format(iter), fcost_all)


    def plot_results(self, vp, grad):
        ''' Plot results during FWI iterations
        '''
        # plot model
        plot_model(self.solver.model.x, 
            self.solver.model.z,
            vp.reshape(self.solver.model.nx, self.solver.model.nz).T, 
            self.optimizer.vp_min, 
            self.optimizer.vp_max,
            os.path.join(self.solver.system.path, 'fwi/figures/vp_it_{:04d}.png'.format(self.optimizer.cpt_iter+1)), 
            'vp', 
            figaspect = 1, 
            colormap = 'jet')

        # plot gradient
        plot_model(self.solver.model.x, 
            self.solver.model.z,
            grad.reshape(self.solver.model.nx, self.solver.model.nz).T, 
            -self.optimizer.update_vpmax, 
            self.optimizer.update_vpmax,
            os.path.join(self.solver.system.path, 'fwi/figures/grad_it_{:04d}.png'.format(self.optimizer.cpt_iter+1)), 
            'grad', 
            figaspect = 1, 
            colormap = 'seismic')

        # plot gradient mask on first iteration
        if self.optimizer.cpt_iter == 0:
            plot_model(self.solver.model.x, 
                self.solver.model.z,
                self.optimizer.grad_mask.reshape(self.solver.model.nx, self.solver.model.nz).T, 
                0, 
                1,
                os.path.join(self.solver.system.path, 'fwi/figures/grad_mask.png'), 
                'grad mask', 
                figaspect = 1, 
                colormap = 'gray')

        # plot_misfit
        plot_misfit(self.solver.system.path, 
            self.optimizer.method, 
            self.optimizer.cpt_iter+1, 
            self.optimizer.niter_max, 
            self.solver.source.num)

        # plot waveform comparison
        


    def __check_obs_data(self):
        ''' Check the existence of observed data
        '''

        for isrc in range(self.solver.source.num):
            sg_file = self.solver.system.path + 'data/obs/src' + str(isrc+1) + '/sg'
            
            if (not os.path.exists(sg_file + '.segy') and 
                not os.path.exists(sg_file + '.su') and 
                not os.path.exists(sg_file + '.bin')) :
                msg = 'FWI workflow ERROR: observed data are not found: {}.segy (.su or .bin)'.format(sg_file)
                raise ValueError(msg)

        print('FWI workflow: find observed data  in {}data/obs/'.format(self.solver.system.path))
        print('FWI workflow: start iteration ...\n')


    def __build_dir(self):
        ''' Build directories for FWI workflow and clean up the previous results if any
        '''

        # print the working path
        path = self.solver.system.path
        print('\nFWI workflow: the working path is in {}'.format(path + 'fwi'))

        # build required directories and clean up the previous results if any
        folders = [path + 'fwi', 
                   path + 'fwi/grad',
                   path + 'fwi/model',
                   path + 'fwi/misfit',
                   path + 'fwi/waveform',
                   path + 'fwi/figures',]

        for folder in folders:
            if os.path.exists(folder):
                os.system('rm -rf ' + folder)
                print('FWI workflow: clean previous data in {}'.format(path + 'fwi'))
            os.makedirs(folder)


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


class Configuration(object):
    ''' Configuration class converts dictionary to object attributes
    '''
    def __init__(self, dict):
        for _, value in dict.items():
            for k, v in value.items():
                setattr(self, k, v)
