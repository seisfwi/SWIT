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


import numpy as np


class FWI(object):
    ''' Full waveform inversion workflow
    '''

    def __init__(self, solver, processor, optimizer):
        ''' Initialize FWI workflow
        '''

        self.receiver = solver
        self.processor = processor
        self.optimizer = optimizer


    def calculate_adjoint_source(self):
        ''' Calculate adjoint source
        '''

    def calculate_gradient(self):
        ''' Calculate gradient of objective function
        '''

    def calculate_misfit(self):
        ''' Calculate misfits
        '''

    def calculate_gradient_misfit(self, vp):
        ''' Calculate gradient of misfit function

        Parameters
        ----------
        vp : flatted 1D array (float32)
            velocity model
        rho : flatted 1D array (float32)
            density model (optional)
        '''
        
        # set model
        self.solver.set_model(vp = vp)

        # generate sythetic data (save boundary data for wavefield reconstruction)
        self.solver.run(simu_type = 'forward', simu_tag = 'syn', save_boundary = True)

        # process synthetic data
        self.processor.run()
        
        # calculate adjoint source and save to work_path/config.path/wavefield/src1_adj.bin, src2_adj.bin, ...
        self.calculate_adjoint_source()

        # calculate misfit
        self.calculate_misfit()

        # calculate gradient via adjoint modeling
        self.solver.run(simu_type = 'gradient', simu_tag = 'syn')

        # process gradient
        self.postprocess_gradient()

        return np.sum(self.misfit), self.misfit, self.gradient


    def save_results(self, count, vp, grad, fcost_all):
        ''' Save results during FWI iterations
        '''
        
        np.save(self.solver.path + 'outputs/grad_ite_{:04d}.bin'.format(count), grad.flatten())
        np.save(self.solver.path + 'outputs/vp_ite_{:04d}.bin'.format(count), vp.flatten())
        np.save(self.solver.path + 'outputs/misfit_iter_{:04d}.bin'.format(count), fcost_all)
        

    def run(self):
        ''' Run FWI workflow
        '''

        # set the initial model
        vp = self.optimizer.initial_model['vp']
        self.solver.set_model(vp = vp)

        # calculate gradient and misfit
        fcost, fcost_all, grad = self.calculate_gradient_misfit()

        count = 0
        # save initial misfit and gradient
        self.save_results(count, vp, grad, fcost_all)

        print("Initial fcost: {} \n".format(fcost))

        # keep iterate while convergence not reached or linesearch not failed
        while ((self.optimizer.FLAG != 'CONV') and (self.optimizer.FLAG != 'FAIL')):
            
            vp = self.optimizer.iterate(vp, fcost, grad, grad_preco)
 
            if(self.optimizer.FLAG == 'GRAD'):
                
                # compute cost and gradient of model
                count += 1
                fcost, fcost_all, grad = self.gradient_misfit(vp)
                
                # precondition gradient
                grad_preco = np.copy(grad)

                # save iteration history
                self.save_results(count, vp, grad, fcost_all)


        print('END OF FWI')
        print('See the convergence history in iterate_{}.dat'.format(self.optimizer.scheme))