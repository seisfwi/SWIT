###############################################################################
#
# SWIT: Seismic Waveform Inversion Toolbox
#
# by Haipeng Li at USTC, haipengl@mail.ustc.edu.cn
# 
# June, 2021  
#
# Base module
#
###############################################################################

import os
from multiprocessing import cpu_count

import numpy as np


class simulate(object):
    ''' simulate class defining all parameters for seismic wavefield modeling
    '''

    def __init__(self, model, source, receiver, system):
        ''' define all parameters here
        '''
        self.model = model
        self.source = source
        self.receiver = receiver
        self.system = system
        self.initialize()      # initialize


    def __check(self):
        ''' check all parameters are acceptable
        '''
        dx = self.model.dx
        dt = self.model.dt
        f0 = self.source.f0

        vpmin = np.min(self.model.vp)
        vpmax = np.max(self.model.vp)
        dt_req = np.sqrt(3.0/8.0) * dx / vpmax
        dx_req = vpmin / f0 / 10.
        f0_req = vpmin / dx / 10.

        # check model size
        if np.shape(self.model.vp) != (self.model.nx, self.model.nz):
            raise ValueError('Wrong dimensions of Vp, not consistant with nx = %d, nz = %d'%(self.model.nx, self.model.nz))
       
        if np.shape(self.model.rho) != (self.model.nx, self.model.nz):
            raise ValueError('Wrong dimensions of Rho, not consistant with nx = %d, nz = %d'%(self.model.nx, self.model.nz))

        # check path format
        if self.system.homepath[-1] != '/':
            self.system.homepath += '/'

        # Check the stable condition: 4-th order FD: dt <= sqrt(3/8) * dx / vmax
        if dt_req <= dt:
            raise ValueError('modeling stability: dt = %.4f ms > dt_required = %.4f ms: ' % (dt*1000, dt_req*1000))

        # Check the dispersion condition: dx <= vmin/(10*f0)
        if dx_req < dx:
            print('Warning: modeling dispersion, dx = %.2f m > dx_required =  %.2f m' %(dx, dx_req))
            print('Warning: modeling dispersion, f0 = %.2f Hz > f0_required = %.2f Hz' %(f0, f0_req))
  
        if (self.source.xz[:,0].min() < self.model.xx.min() or 
            self.source.xz[:,0].max() > self.model.xx.max() or  
            self.receiver.xz[:,0].min() < self.model.xx.min() or 
            self.receiver.xz[:,0].max() > self.model.xx.max()):
        
            raise ValueError('source or receiver coordinates are out of range')
        
        # Check the system and set the number of CPUs
        self.system.mpiproc = min([self.system.mpiproc, self.source.n, cpu_count() // 2])
        
        # Use one thread in calling scipy to do the filtering
        os.environ["OMP_NUM_THREADS"] = "1" # export OMP_NUM_THREADS=1


    def __builddir(self):
        ''' Build and clean the working folders
        '''

        homepath = self.system.homepath
        # clean the previous data if exists
        folder_list = [homepath+'data/syn',
                       homepath+'data/tempdata',
                       homepath+'parfile',
                       homepath+'outputs',
                       homepath+'figures',
                    ]
        for _, ifolder in enumerate(folder_list):
            if os.path.exists(ifolder):
                os.system('rm  -r %s' % ifolder)

        # creat the working folders
        folder_list = ['-p ' + homepath,                     # Home folder
                       homepath+'data',                      # Data folder
                       homepath+'data/obs/',                 # Observed data
                       homepath+'data/syn/',                 # Synthetic data
                       homepath+'data/tempdata/',            # Tempdata data
                       homepath+'parfile',                   # Parfile
                       homepath+'parfile/forward_parfile/',  # Forward parfile
                       homepath+'parfile/forward_source/',   # Forward source
                       homepath+'parfile/adjoint_parfile/',  # Adjoint parfile
                       homepath+'parfile/adjoint_source/',   # Adjoint source
                       homepath+'parfile/model/',            # Model
                       homepath+'outputs',                   # Outputs
                       homepath+'outputs/velocity/',         # Outputs: velocity
                       homepath+'outputs/gradient/',         # Outputs: gradient
                       homepath+'outputs/direction/',        # Outputs: direction
                       homepath+'outputs/LBFGS_memory/',     # Outputs: LBFGS_memory
                       homepath+'figures',                   # Figures
                       homepath+'figures/model/',            # Figures: gradient and velocity
                       homepath+'figures/waveform/',         # Figures: waveform
                       ]
        for _, ifolder in enumerate(folder_list):
            if not os.path.exists(ifolder):
                os.system('mkdir %s' % ifolder)


    def __screenprint(self):
        ''' Print information to screen
        '''

        # print the primary modeling parameters
        print('*****************************************************')
        print('\n        Seismic Waveform Inversion Toolbox         \n')
        print('*****************************************************\n')
        print('Forward modeling : nx, nz = %d, %d' %(self.model.nx, self.model.nz) )
        print('Forward modeling : dx = %.1f m'             %(self.model.dx))
        print('Forward modeling : dt = %.2f ms, %d steps'  %(self.model.dt * 1000, self.model.nt))
        print('Forward modeling : %d shots run in mpi, %d CPU available'%(self.system.mpiproc, cpu_count() // 2))


    def initialize(self):
        ''' Initialize the simulation
        '''

        self.__check()              # Check simulation parameter
        self.__builddir()           # Build and clean working folder
        self.__screenprint()        # Print to screen


class model(object):
    ''' model parameters for wavefiled simulation (2D or 3D)
    '''

    def __init__(self, nx, nz, dx, dt, nt, fs, pml, vp, rho):

        # basic model
        self.nx = nx
        self.nz = nz
        self.dx = dx
        self.dt = dt
        self.nt = nt
        self.fs = fs                                           
        self.pml = pml
        self.nx_pml = self.nx + self.pml * 2
        self.nz_pml = self.nz + self.pml * (2 - self.fs)
        # coordinate points (x, z), and the time array
        self.xx = np.arange(0, self.nx * self.dx, self.dx)
        self.zz = np.arange(0, self.nz * self.dx, self.dx)
        self.t  = np.linspace(0, self.dt*self.nt, num=self.nt, endpoint=False)
        # velocity, density, etc.
        self.rho = rho                                         # density in kg/m^3
        self.vp = vp                                           # p-wave velocity in m/s
        self.vpmax = vp.max()                                  # maximum vp
        self.vpmin = vp.min()                                  # mininum vp
        # output
        self.savesnap = 0
        self.savestep = 1

        # Notice:
        # the forward modeling code always uses free surface, no matter self.fs = True or False
        #    

class source(object):
    ''' source parameters for forward wavefield simulation (2D or 3D)
    '''

    def __init__(self, f0, n, xz, wavelet):
        self.f0 = f0
        self.n = n
        self.xz = xz
        self.wavelet = wavelet


class receiver(object):
    ''' receiver parameters for forward wavefield simulation (2D or 3D)
    '''

    def __init__(self, n, xz):
        self.n = n
        self.xz = xz


class system(object):
    ''' system setting
    '''

    def __init__(self, homepath, mpiproc, figaspect=1):
        self.homepath = homepath
        self.mpiproc = mpiproc
        self.figaspect = figaspect


# optimize class
class optimize(object):
    ''' FWI optimation parameters
    '''
 
    def __init__(self, misfit_type, scheme, maxiter, step_length, vpmax, vpmin, marine_or_land,
                 grad_mute, grad_smooth,
                 fre_filter, fre_low, fre_high, 
                 mute_late_arrival, mute_late_window, normalize,
                 mute_offset_short, mute_offset_long, 
                 mute_offset_short_dis, mute_offset_long_dis, grad_mask=None):
        ''' Define all parameters
        '''

        # basic inversion parameters
        self.iter = 0
        self.misfit_type = misfit_type
        self.scheme = scheme
        self.maxiter = maxiter
        self.step_length = step_length
        self.vpmax = vpmax
        self.vpmin = vpmin
        self.marine_or_land = marine_or_land
        # gradient preconditioning
        self.grad_mute = grad_mute
        self.grad_smooth = grad_smooth
        self.grad_mask = grad_mask
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
        self.mute_offset_long  = mute_offset_long
        self.mute_offset_short_dis = mute_offset_short_dis           # (units: m)
        self.mute_offset_long_dis  = mute_offset_long_dis            # (units: m)
        
        # set the taper for muting the gradient around the source 
        if self.marine_or_land.lower() in ['marine', 'offshore']:
            self.grad_thred = 0.0
        elif self.marine_or_land.lower() in ['land', 'onshore']:
            self.grad_thred = 0.001
        else:
            raise ValueError('not supported modeling marine_or_land: %s'%(self.marine_or_land))

        # initilize
        self.initialize()


    def __check(self):
        ''' check parameters.
        '''
        
        if self.scheme not in ['NLCG', 'LBFGS']:
            raise ValueError('not supported inversion scheme: %s' % self.scheme)

        if self.misfit_type not in ['Waveform', 'Envelope', 'Traveltime', 'Globalcorrelation', 'RTM']:
            raise ValueError('not supported misfit function: %s' % self.misfit_type)

        if self.misfit_type in ['RTM']:
            self.maxiter = 1

        if ('Max-Trace' not in self.normalize and 
            'L1-Event' not in self.normalize and 
            'L2-Event' not in self.normalize and 
            'L1-Trace' not in self.normalize and 
            'L2-Trace' not in self.normalize and 
            'None' not in self.normalize):
            raise ValueError('not supported normalization:', self.normalize)

        if self.fre_filter not in ['Lowpass', 'Bandpass', 'Highpass', 'None']:
            raise ValueError('not supported frequency filter: %s' %self.fre_filter)

        if self.vpmax < self.vpmin:
            raise ValueError('vpmax=%f m/s is less than vpmin=%f m/s\n' %(self.vpmax, self.vpmin))

        if self.fre_low > self.fre_high:
            raise ValueError('fre_low > fre_high')

        if self.mute_offset_short_dis > self.mute_offset_long_dis:
            raise ValueError('mute_offset_short_dis > mute_offset_long_dis')


    def __screenprint(self):
        ''' Print information to screen.
        '''
        # basic inversion parameter
        print('Inversion scheme : %s' % self.scheme)
        print('Inversion misfit : %s' % self.misfit_type)
        print('Inversion maxiter: %d' % self.maxiter)
        print('Inversion step   : %.3f'      % self.step_length)
        print('Inversion vpmin  : %.1f m/s'  % self.vpmin)
        print('Inversion vpmax  : %.1f m/s'  % self.vpmax)
        print('Gradient  mute   : %d grids on top' % self.grad_mute)
        print('Gradient  smooth : %d grids Gaussian smooth' % self.grad_smooth)

        # filtering 
        if self.fre_filter in ['None']:
            print('Data processing  : no filtering')
        elif self.fre_filter in ['Bandpass']:
            print('Data processing  : %s, %.2f ~ %.2f Hz' %(self.fre_filter, self.fre_low, self.fre_high))
        elif self.fre_filter in ['Lowpass']:
            print('Data processing  : %s, < %.2f Hz' %(self.fre_filter, self.fre_low))
        elif self.fre_filter in ['Highpass']:
            print('Data processing  : %s, > %.2f Hz' %(self.fre_filter, self.fre_high))

        # pick and mute
        if self.mute_late_arrival:
            print('Data processing  : time window, %.2f s after the first break' % self.mute_late_window)
        else:
            print('Data processing  : no time window')

        # mute offset
        if self.mute_offset_short:
                print('Data processing  : mute short offset, %.1f m' % self.mute_offset_short_dis)
        else:
            print('Data processing  : mute short offset, none')
        if self.mute_offset_long:
                print('Data processing  : mute long offset, %.1f m' % self.mute_offset_long_dis)
        else:
            print('Data processing  : mute long offset, none')
        
        # normalize
        if self.normalize in ['None']:
            print('Data processing  : no normalization (keep AVO effect)')
        else:
            #print('Data processing  : normalization, %s' % (', '.join(self.normalize)))
            print('Data processing  : normalization, %s' % (self.normalize))
        print('Data processing  : OMP Threads = %s' %os.environ["OMP_NUM_THREADS"])
        
        print('\nsee more in json-parameter files under parfile folder\n')
        
        print('*****************************************************\n')


    def initialize(self):
        ''' Initialize the optimazation.
        '''

        self.__check()              # Check simulation parameter
        self.__screenprint()
