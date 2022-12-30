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

# import modules
import os
import time
import numpy as np
from multiprocessing import Pool

# import SWIT modules
from misfit import calculate_adjoint_misfit_is
from tools import load_float, smooth2d
from utils import generate_preconditioner, generate_mask
from plot import plot_model, plot_misfit


class FWI(object):
    ''' Full waveform inversion workflow

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

        self.solver = solver
        self.optimizer = optimizer
        self.preprocessor = preprocessor

        # build directories
        self.__build_dir()

        # check the existence of obs data
        self.__check_obs_data()
        

    def __check_obs_data(self):
        ''' Check the existence of observed data
        '''

        for isrc in range(self.solver.source.num):
            sg_file = self.solver.config.path + 'data/obs/src' + str(isrc+1) + '/sg'
            
            if (not os.path.exists(sg_file + '.segy') and 
                not os.path.exists(sg_file + '.su') and 
                not os.path.exists(sg_file + '.bin')) :
                msg = 'FWI workflow ERROR: observed data are not found: {}.segy (.su or .bin)'.format(sg_file)
                raise ValueError(msg)

        print('FWI workflow: find observed data  in {}data/obs/'.format(self.solver.config.path))
        print('FWI workflow: start iteration ...\n')


    def __build_dir(self):
        ''' Build directories for FWI workflow and clean up the previous results if any
        '''

        path = self.solver.config.path
        folders = [path + 'fwi', 
                   path + 'fwi/grad',
                   path + 'fwi/model',
                   path + 'fwi/misfit',
                   path + 'fwi/waveform',
                   path + 'fwi/figures',]

        # print message
        print('\n')
        # clean up the previous results and build directories
        for folder in folders:
            if os.path.exists(folder):
                os.system('rm -rf ' + folder)
                print('FWI workflow: clean previous data in {}'.format(path + 'fwi'))
            os.makedirs(folder)
        print('FWI workflow: build working paths in {}'.format(path + 'fwi'))


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

        # block at this line until all processes are done
        pool.join()
        
        # sum the misfit over all sources
        fcost = 0.0
        for _fcost in fcost_all_src:
            fcost += _fcost

        # return the summed misfit over all sources
        return np.array(fcost_all_src), fcost


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
        self.preprocessor.run(data_path = self.solver.config.path + 'data/syn/', 
                            src_num = self.solver.source.num, 
                            mpi_num = self.solver.config.mpi_num, 
                            nt = self.solver.model.nt, 
                            dt = self.solver.model.dt,  
                            src_coord = self.solver.source.coord, 
                            rec_coord = self.solver.receiver.coord)
        
        # prepare adjoint source and calculate misfit
        fcost_all_src, fcost = self.prepare_adjoint_misfit()

        # calculate gradient via adjoint modeling
        self.solver.run(simu_type = 'gradient', simu_tag = 'syn')

        # process gradient
        grad = self.postprocess_gradient()

        # save misfit and gradient for each iteration
        return fcost_all_src, fcost, grad


    def postprocess_gradient(self):
        ''' Postprocess gradient
            1. read gradient from binary files and sum over all sources
            2. apply mask
            3. apply preconditioning
            4. apply smoothness
            5. apply normalization with proper scaling factor

        Returns
        -------
        grad : 2D array (float32)
            gradient of misfit function
        '''
        # get the parameters
        nx = self.solver.model.nx
        nz = self.solver.model.nz 
        src_num  = self.solver.source.num
        acquisition_type = self.solver.model.acquisition_type

        # load gradient and wavefield illumaition from binary files
        grad = np.zeros((nx, nz))
        for_illum = np.zeros((nx, nz))
        adj_illum = np.zeros((nx, nz))

        # sum over all sources
        for isrc in range(src_num):
            path = os.path.join(self.solver.config.path, 'data/syn/src{}/'.format(isrc+1))
            grad += load_float(path + 'vp_gradient.bin').reshape(nx, nz)
            for_illum += load_float(path + 'forward_illumination.bin').reshape(nx, nz)
            adj_illum += load_float(path + 'adjoint_illumination.bin').reshape(nx, nz)

        # generate a default mask or use the provided mask
        if self.optimizer.grad_mask is None:
            mask = generate_mask(nx, nz, acquisition_type, threshold = 0.05, mask_size = 10)
        else:
            mask = self.optimizer.grad_mask
        grad *= mask

        # TODO: generete the mask only at the beginning of the inversion
        # apply the preconditioning, which is the approximated inverse Hessian
        precond = generate_preconditioner(for_illum, adj_illum)
        grad = grad / np.power(precond, 1.)

        # apply smoothness to the gradient, if smoothness is provided
        if self.optimizer.grad_smooth_size > 0:
            grad = smooth2d(grad, span=self.optimizer.grad_smooth_size)

        # scale the gradient properly
        grad = grad / abs(grad).max() * self.optimizer.update_vpmax

        return grad


    def run(self):
        ''' Run FWI workflow
        '''

        # start the timer
        start_time = time.time()

        # preprocess the obs data
        self.preprocessor.run(data_path = self.solver.config.path + 'data/obs/', 
                            src_num = self.solver.source.num, 
                            mpi_num = self.solver.config.mpi_num, 
                            nt = self.solver.model.nt, 
                            dt = self.solver.model.dt,  
                            src_coord = self.solver.source.coord, 
                            rec_coord = self.solver.receiver.coord)

        # calculate gradient and misfit from initial model
        vp = self.optimizer.vp_init
        rho = self.optimizer.rho_init
        fcost_all_src, fcost, grad = self.calculate_gradient_misfit(vp = vp, rho = rho)

        # keep iterate while convergence not reached or linesearch not failed
        while ((self.optimizer.FLAG != 'CONV') and (self.optimizer.FLAG != 'FAIL')):
            # TODO:  plot initial model
            # TODO:  simplicify the code below
            # update model and linesearch using the preconditioned gradient
            grad_preco = np.copy(grad)
            vp_old = np.copy(vp).flatten()
            vp = self.optimizer.iterate(vp, fcost, grad, grad_preco)

            # calculate gradient and misfit from updated model
            if(self.optimizer.FLAG == 'GRAD'):
                # print the iteration information
                print('Iteration: {} \t fcost: {:.4e} \t step: {:7.4f} \t max update: {:7.2f} m/s'.format(
                    self.optimizer.cpt_iter+1, fcost, self.optimizer.alpha, np.max(np.abs(vp-vp_old))))

                # compute cost and gradient
                fcost_all_src, fcost, grad = self.calculate_gradient_misfit(vp = vp, rho = rho)
                
                # save the iteration history
                self.save_results(vp, grad, fcost_all_src)

                # plot the model, gradient and misfit
                self.plot_results(vp, grad)

        # print the end information
        hours, rem = divmod(time.time()-start_time, 3600)
        minutes, seconds = divmod(rem, 60)
        print('\nFWI workflow: job is finished.')
        print('FWI workflow: running time is for {:0>2}h {:0>2}m {:.0f}s'.format(int(hours), int(minutes), seconds))
        print('FWI workflow: see convergence history in iterate_{}.log\n'.format(self.optimizer.method))


    def save_results(self, vp, grad, fcost_all):
        ''' Save results during FWI iterations
        '''
        np.save(self.solver.config.path + 'fwi/grad/grad_it_{:04d}.npy'.format(self.optimizer.cpt_iter+1), grad)
        np.save(self.solver.config.path + 'fwi/model/vp_it_{:04d}.npy'.format(self.optimizer.cpt_iter+1), vp)
        np.save(self.solver.config.path + 'fwi/misfit/fcost_all_it_{:04d}.npy'.format(self.optimizer.cpt_iter+1), fcost_all)


    def plot_results(self, vp, grad):
        ''' Plot results during FWI iterations
        '''
        # plot model
        plot_model(self.solver.model.x, 
            self.solver.model.z,
            vp.reshape(self.solver.model.nx, self.solver.model.nz).T, 
            self.optimizer.vp_min, 
            self.optimizer.vp_max,
            os.path.join(self.solver.config.path, 'fwi/figures/vp_it_{:04d}.png'.format(self.optimizer.cpt_iter+1)), 
            'vp', 
            figaspect = 1, 
            colormap = 'jet')

        # plot gradient
        plot_model(self.solver.model.x, 
            self.solver.model.z,
            grad.reshape(self.solver.model.nx, self.solver.model.nz).T, 
            -self.optimizer.update_vpmax, 
            self.optimizer.update_vpmax,
            os.path.join(self.solver.config.path, 'fwi/figures/grad_it_{:04d}.png'.format(self.optimizer.cpt_iter+1)), 
            'grad', 
            figaspect = 1, 
            colormap = 'seismic')

        # plot gradient mask
        if self.optimizer.grad_mask is not None and self.optimizer.cpt_iter == 0:
            plot_model(self.solver.model.x, 
                self.solver.model.z,
                self.optimizer.grad_mask.reshape(self.solver.model.nx, self.solver.model.nz).T, 
                0, 
                1,
                os.path.join(self.solver.config.path, 'fwi/figures/grad_mask.png'), 
                'grad mask', 
                figaspect = 1, 
                colormap = 'gray')

        # plot_missfit
        plot_misfit(self.solver.config.path, 
            self.optimizer.method, 
            self.optimizer.cpt_iter+1, 
            self.optimizer.niter_max, 
            self.solver.source.num)


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