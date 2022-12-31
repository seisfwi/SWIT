###############################################################################
#
# SWIT v1.1: Seismic Waveform Inversion Toolbox
#
#   A Python package for seismic waveform inversion
#   By Haipeng Li at USTC & Stanford
#   Email: haipengl@mail.ustc.edu.cn, haipeng@stanford.edu 
#
#   Survey class for seismic acquisition, including sub-classes: System, Model, 
#   Source, and Receiver
#
###############################################################################

import os
from multiprocessing import cpu_count

import matplotlib.pyplot as plt
import numpy as np


class System(object):
    ''' System class describes the system configuration

        parameters:
        ----------
            path: str
                path to perform the modeling or inversion
            mpi_num: int
                number of mpi processes
            max_cpu_num: int
                maximum number of CPUs on the PC/cluster
            fig_aspect: float
                aspect ratio of the figure (default: 1.0)
    '''

    def __init__(self, path, mpi_num, max_cpu_num = cpu_count()//2, fig_aspect = 1.0):
        ''' initialize system class
        '''

        # basic parameters
        self.path = path
        self.mpi_num = mpi_num
        self.max_cpu_num = max_cpu_num
        self.fig_aspect = fig_aspect
        
        # check the parameters
        self.__check__()

 
    def __check__(self):
        ''' check the config parameters
        '''
        # check path format
        if self.path[-1] != '/':
            self.path += '/'

        # check the existence of the working path
        if not os.path.exists(self.path):
            os.makedirs(self.path)

        # check the number of mpi processes
        if self.mpi_num > self.max_cpu_num:
            msg  = 'Survey Warning: number of mpi processes {} exceeds the number of CPUs\n'.format(self.mpi_num)
            msg += 'Survey Warning: number of mpi processes is reset to {}.'.format(self.max_cpu_num)
            print(msg)
            self.mpi_num = self.max_cpu_num
    

class Model(object):
    ''' Model class describes the model for wavefield simulation

        parameters:
        ----------
            nx: int
                number of grid points in x direction
            nz: int 
                number of grid points in z direction
            dx: float
                grid spacing in x and z directions
            dt: float   
                time step, in seconds
            nt: int
                number of time steps
            pml: int
                width of pml boundary
            vp: 2d array    
                p-wave velocity model with dimensions of nx * nz, in m/s
            rho: 2d array
                density model with dimensions of nx * nz, in kg/m^3
            topo_coord: 2d array (optional)
                coordinates of topography (x, z), in meters
    '''

    def __init__(self, nx, nz, dx, dt, nt, pml, vp, rho):
        ''' initialize model class
        '''

        # basic parameters
        self.nx = nx
        self.nz = nz
        self.dx = dx
        self.dt = dt
        self.nt = nt
        self.pml = pml
        self.rho = np.copy(rho)
        self.vp = np.copy(vp)

        # pml parameters, free surface at the top
        self.nx_pml = self.nx + self.pml * 2
        self.nz_pml = self.nz + self.pml

        # coordinates of grid points
        self.x = np.arange(0, self.nx * self.dx, self.dx)
        self.z = np.arange(0, self.nz * self.dx, self.dx)

        # time axis
        self.t = np.linspace(0, self.dt*self.nt, num=self.nt, endpoint=False)


        ## For now, we do not consider topography
        # # set the topography
        # if topo_coord is None:
        #     self.topo = np.zeros(self.nx)
        # else:
        # # interpolate the topography and reset the largest value to be zero
        #     self.topo = np.interp(self.x, topo_coord[:,0], topo_coord[:,1])
        #     self.topo = - self.topo + np.max(self.topo)
        
        # default parameters for output wavefield
        self.save_snap = 0
        self.save_step = 10

        # check the parameters
        self.__check__()


    def __check__(self):
        ''' check the model parameters
        '''

        # check the dimensions of models
        if np.shape(self.vp) != np.shape(self.rho) != (self.nx, self.nz):
            raise ValueError('Survey Error: the dimensions of vp/rho are not consistant \with nx = {}, nz = {}'.format(self.nx, self.nz))
       
        # check the pml width
        if self.pml < 20:
            print('Survey Warning: the width of pml boundary is recommended to be larger than 20 grid points')

        # # check the topography
        # if len(self.topo) != self.nx:
        #     raise ValueError('Survey: the length of topography is not consistant with nx = {}'.format(self.nx))

        # if np.max(self.topo) > self.z[-1]:
        #     raise ValueError('Survey: the topography exceeds the model range in depth')


class Source(object):
    ''' Source class describes the source geometry

        parameters:
        ----------
            coord: 2d array
                coordinates of sources, in meters
            wavelet: 2d array
                source wavelet with dimensions of [src_num, nt], in Pa
            f0: float
                dominant frequency of the source wavelet, in Hz
    '''

    def __init__(self, coord, wavelet, f0):
        ''' initialize source class
        '''

        # basic parameters
        self.coord = coord
        self.wavelet = wavelet
        self.f0 = f0

        # number of sources
        self.num = len(self.coord)

        # check the parameters
        self.__check__()


    def __check__(self):
        ''' check the source parameters
        '''

        # check coordinates of sources
        if self.coord.shape != (self.num, 2):
            raise ValueError('Survey: the source coordinates must be 2D array with dimensions of [src_num, 2]')

        # check the wavelet
        if self.wavelet.shape[0] != self.num:
            raise ValueError('Survey: the number of source wavelets must be the same as the number of sources')


class Receiver(object):
    ''' Receiver class describes the receiver geometry

        parameters:
        ----------
            coord: 2d array
                coordinates of receivers, in meters
            comp: str
                components of receivers, 'vx', 'vz', 'p'
    '''

    def __init__(self, coord, comp):
        ''' initialize receiver class
        '''

        # basic parameters
        self.coord = coord
        self.comp = comp

        # check the parameters
        self.__check__()


    def __check__(self):
        ''' check the receiver parameters
        '''
        # check coordinates of receivers
        if not isinstance(self.coord, list):
            raise ValueError('Survey: receiver coordinates must be a list of the length of the source number, with a 2D array for each element')

        # check the components of receivers
        if self.comp.lower() not in ['vx', 'vz', 'p']:
            msg = 'Survey: receiver component must be vx, vz or p'
            err = 'Unsupport receiver component: {}'.format(self.comp) 
            raise ValueError(msg + '\n' + err)



class Survey(object):
    ''' Survey class describes the seismic acquisition geometry

    parameters:
    ----------
        system: System class
            system configuration
        model: Model class
            model parameters
        source: Source class
            source parameters
        receiver: Receiver class
            receiver parameters
    '''

    def __init__(self, system, model, source, receiver):
        ''' initialize survey class
        '''

        # basic parameters
        self.system = system
        self.model = model
        self.source = source
        self.receiver = receiver

        # check the parameters
        self.__check__()

        # set the derived parameters
        self.__set_derived_parameters__()

        # plot the survey
        self.plot_geometry()
        self.plot_wavelet(isrc = 1, t_max = 1.0)
        self.plot_model()


    def __check__(self):
        ''' check the the overall consistency of the parameters
        '''

        # check the stability condition (4-th order FD): dt <= sqrt(3/8) * dx / vmax
        dt0 = np.sqrt(3.0/8.0) * self.model.dx / np.max(self.model.vp)
        if dt0 <= self.model.dt:
            raise ValueError('Survey Error: the stability condition of 4-th order FD method \
            is not satisfied: dt = %.4f s > dt_required = %.4f s'.format(self.model.dt, dt0))

        # check the numerical dispersion condition: dx <= vmin/(10*f0)
        dx0 = np.min(self.model.vp) / self.source.f0 / 10.
        f00 = np.min(self.model.vp) / self.model.dx  / 10.

        if dx0 < self.model.dx:
            print('Survey Warning: numerical dispersion, dx = {:6.2f} m  > dx_required = {:6.2f} m'.format(self.model.dx, dx0))
            print('Survey Warning: numerical dispersion, f0 = {:6.2f} Hz > f0_required = {:6.2f} Hz'.format(self.source.f0, f00))
  
        # check the source location
        if (self.source.coord[:,0].min() < self.model.x.min() or 
            self.source.coord[:,0].max() > self.model.x.max() or
            self.source.coord[:,1].min() < self.model.z.min() or 
            self.source.coord[:,1].max() > self.model.z.max()):
            raise ValueError('Survey Error: source location is out of model range')
        
        # check the receiver coord list
        if len(self.receiver.coord) != self.source.num:
            raise ValueError('Survey Error: receiver coord list is not consistent with the source number')

        # check the receiver location
        for rec_coord in self.receiver.coord:
            if (rec_coord[:, 0].min() < self.model.x.min() or 
                rec_coord[:, 0].max() > self.model.x.max() or
                rec_coord[:, 1].min() < self.model.z.min() or 
                rec_coord[:, 1].max() > self.model.z.max()):
                raise ValueError('Survey Error: receiver location is out of model range')
        
        ## For now, we do not have the topography
        # # check topography
        # for isrc in range(self.source.num):
        #     # check the source location is below the topography
        #     ix = np.argmin(np.abs(self.model.x - self.source.coord[isrc, 0]))
        #     if self.source.coord[isrc, 1] < self.model.topo[ix]:
        #         raise ValueError('Survey: source {} is above the topography'.format(isrc))

        #     # check the receiver location is below the topography
        #     for irec in range(self.receiver.coord[isrc].shape[0]):
        #         ix = np.argmin(np.abs(self.model.x - self.receiver.coord[isrc][irec, 0]))
        #         if self.receiver.coord[isrc][irec, 1] < self.model.topo[ix]:
        #             raise ValueError('Survey: receiver {} from source {} is above the topography'.format(irec, isrc))

        # check the mpi number and the source number
        if self.system.mpi_num > self.source.num:
            print('Survey Warning: number of mpi processes is larger than source number, reset to source number')
            self.system.mpi_num = self.source.num

        # check the configfile directory and create it if not exist
        if not os.path.exists(self.system.path + 'system'):
            os.system('mkdir -p %s' % (self.system.path + 'config'))
            os.system('mkdir -p %s' % (self.system.path + 'config/wavelet'))


    def __set_derived_parameters__(self):
        ''' Set the derived parameters
        '''

        # set the offset
        self.model.offset = []
        for i in range(self.source.num):
            self.model.offset.append(self.receiver.coord[i][:,0] - self.source.coord[i,0])

    ## For now, we do not have the topography
    # def plot_topography(self, plot_source = True):
    #     ''' Plot the topography with the source and receiver location (optional)
    #     '''

    #     # plot topography
    #     fig = plt.figure()
    #     fig.add_subplot(111)

    #     x = self.model.x / 1000.
    #     z = self.model.z / 1000.
    #     topo = self.model.topo  / 1000.

    #     fig = plt.figure(figsize=(10, 5))

    #     # plot the source and receiver
    #     if plot_source:
    #         for isrc in range(self.source.num):
    #             plt.plot(self.source.coord[isrc, 0] / 1000., self.source.coord[isrc, 1] / 1000., 'ro', ms = 5)
    #             plt.plot(self.receiver.coord[isrc][:, 0] / 1000., self.receiver.coord[isrc][:, 1] / 1000., 'bo', ms = 2)

    #     # plot the topography
    #     plt.plot(x, topo, 'k', lw = 1.5)
    #     plt.plot([x[ 0], x[ 0]], [topo[ 0], z[-1]], 'k', lw = 1.5)
    #     plt.plot([x[-1], x[-1]], [topo[-1], z[-1]], 'k', lw = 1.5)
    #     plt.plot([x[ 0], x[-1]], [   z[-1], z[-1]], 'k', lw = 1.5)
    #     plt.gca().invert_yaxis()
    #     plt.xlabel('Distance (km)', fontsize=12)
    #     plt.ylabel('Depth (km)', fontsize=12)
    #     plt.title('Topography')
    #     plt.savefig(self.system.path + 'config/topography.png', dpi=300)
    #     plt.close()


    def plot_geometry(self):
        ''' Plot source and receiver geometry in the seismic survey
        '''

        # plot geometry
        fig = plt.figure()
        fig.add_subplot(111)
        for i in range(self.source.num):
            plt.scatter(self.receiver.coord[i][:,0], np.ones(len(self.receiver.coord[i])) * i + 1, c = 'gray', marker='o', s = 2) 
            plt.scatter(self.source.coord[i, 0],  i + 1, c = 'red', marker='*', s = 6)

        plt.ylim(0, self.source.num+1)
        plt.yticks(np.arange(0, self.source.num, 2))
        plt.xlabel('Distance (m)', fontsize=12)
        plt.ylabel('Shot number', fontsize=12)
        plt.title('Survey Geometry, %d sources\n' % (self.source.num), fontsize = 14)

        plt.savefig(os.path.join(self.system.path, 'config/survey_geometry.png'), dpi=300)
        plt.close()
        # plt.show()


    def plot_wavelet(self, isrc = 1, t_max = 1.0):
        ''' Plot source time function in the time and frequency domanin.
    
        Parameters:
        -----------
            isrc: int
                source index, should be in [1, src_num]
            t_max: float
                maximum time to plot in [0, t_max], in second
        '''
        
        # check the source index
        if isrc > self.source.num or isrc < 1:
            raise ValueError('Survey Error: source index is outof range, should be in [1, {}]'.format(self.source.num))

        # check the maximum time
        if t_max > self.model.t[-1] or t_max < 0:
            print('Survey Warning: t_max is out of range, reset to model.t[-1]')
            t_max = self.model.t[-1]

        # get the source time function
        t = self.model.t
        wvlt_time = self.source.wavelet[isrc-1, :]

        # perform fft to get the frequency domain wavelet
        wvlt_freq = np.abs(np.fft.fft(wvlt_time))
        freqs = np.fft.fftfreq(len(wvlt_time), self.model.dt)
        idx = np.argsort(freqs)
        idx = idx[int(len(idx) / 2):]

        # subplot1: time domian wavelet
        fig = plt.figure()
        ax1 = fig.add_subplot(2, 1, 1)
        ax1.xaxis.set_label_text('Time (s)', fontsize=12)
        ax1.yaxis.set_label_text('Norm. Amplitude', fontsize=12)
        ax1.set_title('Source wavelet of source {}'.format(isrc), fontsize=12)
        ax1.axis((0, t[-1], -1.2, 1.2))
        ax1.plot(t, wvlt_time / abs(wvlt_time).max(), 'b-')
        ax1.set_xlim(0, t_max)

        # subplot2: frequency domian wavelet
        ax2 = fig.add_subplot(2, 1, 2)
        ax2.xaxis.set_label_text('Frequency (Hz)', fontsize=12)
        ax2.yaxis.set_label_text('Norm. Amplitude', fontsize=12)
        ax2.set_title('Amplitude spectrum', fontsize=12)
        ax2.axis((0, 50, 0, 1.2))
        ax2.fill(freqs[idx], wvlt_freq[idx] / abs(wvlt_freq).max(), 'c')
        ax2.plot(freqs[idx], wvlt_freq[idx] / abs(wvlt_freq).max(), 'b')
        fig.tight_layout()

        # save figure
        plt.savefig(os.path.join(self.system.path, 'config/wavelet.png'), dpi=300)
        plt.close()
        # plt.show()


    def plot_model(self, vmin = None, vmax = None, cmap = 'jet', fig_path = None):
        ''' Plot vp model in the seismic survey

        Parameters:
        -----------
            vmin: float
                minimum velocity to plot, default is the minimum velocity in the model
            vmax: float
                maximum velocity to plot, default is the maximum velocity in the model
            cmap: str
                colormap to plot
            fig_path: str
                path to save the figure, if None, save to the working path
        '''

        # plot vp model and save figure (could plot other models as well)
        data = self.model.vp.T

        # set figure path
        fig_path = os.path.join(self.system.path, 'config/vp.png') if fig_path is None else fig_path

        # set plot options
        opts = {
        'vmin': data.min() if vmin is None else vmin,
        'vmax': data.max() if vmax is None else vmax,
        'extent': [self.model.x[0]/1000, self.model.x[-1]/1000, self.model.z[-1]/1000, self.model.z[0]/1000],
        'aspect': self.system.fig_aspect,
        'cmap': cmap
        }

        fig = plt.figure(figsize=(10, 5))
        ax = fig.add_subplot(111)
        im = ax.imshow(data, **opts)
        fig.colorbar(im, shrink=0.5, extend='both')
        ax.xaxis.set_label_text('Depth (km)')
        ax.yaxis.set_label_text('Distance (km)')
        ax.set_title('vp', fontsize=14)
        plt.savefig(fig_path, dpi=300)
        plt.close()
