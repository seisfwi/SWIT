###############################################################################
#
# SWIT v1.1: Seismic Waveform Inversion Toolbox
#
#   A Python package for seismic waveform inversion
#   Developed by Haipeng Li at USTC, updated on 2022-12-21 at Stanford
#   haipengl@mail.ustc.edu.cn, haipeng@stanford.edu
#
#   Workflow module
#
###############################################################################



class SWIT(object):

    def __init__(self, config, model, source, receiver, precessor, optimizer):
        self.config = config
        self.model = model
        self.source = source
        self.receiver = receiver
        self.precessor = precessor
        self.optimizer = optimizer

    def initilize(self):
        
        # initilize solver, precessor and optimizer
        self.solver = solver.Solver(self.config, self.model, self.source, self.receiver)
        self.precessor = precessor.Precessor(self.config, self.precessor)
        self.optimizer = optimizer.Optimizer(self.config, self.optimizer)

        # initilize solver, precessor and optimizer
        self.solver.initilize()
        self.precessor.initilize()
        self.optimizer.initilize()


    def run_FWI(self):
        pass

    def run_RTM(self):
        pass

    