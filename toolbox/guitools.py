###############################################################################
#
# SWIT: Seismic Waveform Inversion Toolbox
#
# by Haipeng Li at USTC, haipengl@mail.ustc.edu.cn
# 
# June, 2021  
#
# guitools module, some of codes are from: https://github.com/PySimpleGUI/PySimpleGUI
#
########################################################################################

import base64
import io
import os
from threading import Thread

import matplotlib.pyplot as plt
import numpy as np
import PIL.Image
import PySimpleGUI as sg

import base
from inversion import inversion
from plot import plot_geometry, plot_model2D, plot_stf, plot_trace
from preprocess import process_workflow
from solver import forward, source_wavelet
from tools import loadfile_gui, saveparjson, smooth2d


def switasync(f):
    ''' async
    '''
    def wrapper(*args, **kwargs):
        thr = Thread(target=f, args=args, kwargs=kwargs)
        thr.start()
    return wrapper
 
 
@switasync
def forward_workflow(swit):
    ''' waveform modeling
    '''
    homepath = swit['homepath']              # home path
    mpiproc  = swit['mpiproc']               # mpi process

    # model setup
    nx   = swit['nx']               # Grids number along x axis
    nz   = swit['nz']               # Grids number along z axis
    dx   = swit['dx']               # Grid size
    dt   = swit['dt']               # Time interval
    nt   = swit['nt']               # Time step
    pml  = swit['pml']              # PML layers grid number
    fs   = swit['fs']               # Free surface or not
    vp_true  = loadfile_gui(swit['vp_true_path'], nx, nz)
    rho_true = np.power(vp_true, 0.25) * 310     # Gardner's equation (1974, Geophysics)

    ### sources setup
    f0    = swit['f0']                           # domiant frequency
    srcxz = np.loadtxt(swit['src_coor'])
    srcn  = srcxz.shape[0]                       # Source number along x axis

    wavelet  = np.zeros((srcn, nt))              # source wavelet
    for isrc in range(srcn):
        if  swit['wavelet'] == 'Ricker':         # wavelet type: Ricker or read from file
            wavelet[isrc,:] = source_wavelet(nt, dt, f0, 'ricker')
        else:
            wavelet[isrc,:] = np.loadtxt(swit['wavelet_file'])


    ### receivers setup
    temp = np.loadtxt(swit['rec_coor'])
    recn  = temp.shape[0]                         # receiver number
    recxz = np.zeros((srcn, recn, 2))             # receiver positions
    for isrc in range(srcn):
        recxz[isrc,:,0] = temp[:,0]               # receiver x position (m)
        recxz[isrc,:,1] = temp[:,1]               # receiver z position (m)

    ### simulate setup 
    sys  = base.system(homepath, mpiproc, figaspect=swit['figaspect'])
    mod  = base.model(nx, nz, dx, dt, nt, fs, pml, vp_true, rho_true)
    src  = base.source(f0, srcn, srcxz, wavelet)
    rec  = base.receiver(recn, recxz)
    simu = base.simulate(mod, src, rec, sys)

    ### plots
    plot_geometry(simu)
    plot_stf(simu, isrc=1,  stf_type='forward', t_end = 2.0)
    plot_model2D(simu, vp_true.T, vp_true.min(), vp_true.max(), 'vp-forward', colormap = 'jet')

    ### forward modeling
    forward(simu, simu_type='obs', savesnap=0)
    plot_trace(simu, 'obs', simu_type='obs', suffix='', src_space=1, trace_space=5, scale = 0.8, color='k')

    print('\n-----------  Forward modeling end  -----------\n')


@switasync
def inversion_workflow(swit):
    ''' waveform inversion 
    ''' 
    homepath = swit['homepath']              # home path
    mpiproc  = swit['mpiproc']               # mpi process

    # model setup
    nx   = swit['nx']               # Grids number along x axis
    nz   = swit['nz']               # Grids number along z axis
    dx   = swit['dx']               # Grid size
    dt   = swit['dt']               # Time interval
    nt   = swit['nt']               # Time step
    pml  = swit['pml']              # PML layers grid number
    fs   = swit['fs']               # Free surface or not
    vp_init  = loadfile_gui(swit['vp_init_path'], nx, nz)
    rho_init = np.power(vp_init, 0.25) * 310     # Gardner's equation (1974, Geophysics)

    ### sources setup
    f0    = swit['f0']                           # domiant frequency
    srcxz = np.loadtxt(swit['src_coor'])
    srcn  = srcxz.shape[0]                       # Source number along x axis
    wavelet  = np.zeros((srcn, nt))              # source wavelet
    for isrc in range(srcn):
        if  swit['wavelet'] == 'Ricker':         # wavelet type: Ricker or read from file
            wavelet[isrc,:] = source_wavelet(nt, dt, f0, 'ricker')
        else:
            wavelet[isrc,:] = np.loadtxt(swit['wavelet_file'])

    ### receivers setup
    temp = np.loadtxt(swit['rec_coor'])
    recn  = temp.shape[0]                         # receiver number
    recxz = np.zeros((srcn, recn, 2))             # receiver positions
    for isrc in range(srcn):
        recxz[isrc,:,0] = temp[:,0]               # receiver x position (m)
        recxz[isrc,:,1] = temp[:,1]               # receiver z position (m)

    ### basic inversion parameter
    misfit_type = swit['misfit_type']            # 'Waveform', 'Envelope','Traveltime', 'Globalcorrelation'
    scheme = swit['scheme']                      # 'NLCG','LBFGS'
    maxiter = swit['maxiter']
    step_length = swit['step_length']
    vpmax = swit['vpmax']  
    vpmin = swit['vpmin']  
    marine_or_land = swit['marine_or_land']
    # gradient preconditioning
    grad_mute   = swit['grad_mute']
    grad_smooth = swit['grad_smooth']
    # data filter
    fre_filter = swit['fre_filter']     #'Bandpass', 'Lowpass', 'Highpass', 'None'
    fre_low    = swit['fre_low']
    fre_high   = swit['fre_high']
    # pick first break and mute
    mute_late_arrival = swit['mute_late_arrival']
    mute_late_window = swit['mute_late_window']    # (units: time)
    # mute according to provided parameter
    mute_offset_short  = swit['mute_offset_short']
    mute_offset_long   = swit['mute_offset_long']
    mute_offset_short_dis = swit['mute_offset_short_dis']          # (units: m)
    mute_offset_long_dis  = swit['mute_offset_long_dis']           # (units: m)
    # data normalize
    normalize = swit['normalize']       # 'L1-Trace', 'L2-Trace', 'L1-Event', 'L2-Event', 'None'

    ### simulate setup 
    sys  = base.system(homepath, mpiproc, figaspect=swit['figaspect'])
    mod  = base.model(nx, nz, dx, dt, nt, fs, pml, vp_init, rho_init)
    src  = base.source(f0, srcn, srcxz, wavelet)
    rec  = base.receiver(recn, recxz)
    simu = base.simulate(mod, src, rec, sys)

    ### optimize setup 
    optim = base.optimize(misfit_type, scheme, maxiter, step_length, vpmax, vpmin, marine_or_land,
                 grad_mute, grad_smooth,
                 fre_filter, fre_low, fre_high, 
                 mute_late_arrival, mute_late_window, normalize,
                 mute_offset_short, mute_offset_long, 
                 mute_offset_short_dis, mute_offset_long_dis)

    ### Save parameter as json
    saveparjson(simu, optim)

    plot_geometry(simu)
    plot_stf(simu, isrc=1,  stf_type='in-use', t_end = 2.0)
    plot_model2D(simu, vp_init.T, vpmin, vpmax, 'vp-init', colormap = 'jet')

    ### process obs data
    if swit['field_data_path'] not in ['']:
        print('copying the obs data to the working folder...')
        print('from: %s'%(swit['field_data_path']))
        if swit['field_data_path'] != swit['homepath'] + 'data/obs/':
            os.system('cp %s %s'%(swit['field_data_path']+ '/*.su', swit['homepath'] + 'data/obs/'))

    process_workflow(simu, optim, simu_type='obs')
    plot_trace(simu, 'obs',      simu_type='obs', suffix='',      src_space=1, trace_space=5, scale = 0.8, color='k')
    plot_trace(simu, 'obs-proc', simu_type='obs', suffix='_proc', src_space=1, trace_space=5, scale = 0.8, color='k')

    ### begin inversion
    inversion(simu, optim, {'vp':vp_init,'rho':rho_init})


# push button do something
def prepare_forward_parameter(GUI_Values):
    ''' prepare parameters for forward 
    ''' 
    try:
        # make necessary adjustment
        if GUI_Values['homepath'][-1] not in ['/']:
            GUI_Values['homepath'] += '/'
        
        swit = {}
        # system setup
        swit['homepath'] = GUI_Values['homepath']
        swit['mpiproc']  = int(GUI_Values['mpiproc'])
        swit['figaspect'] = GUI_Values['figaspect']
        # model setup
        swit['nx']  = int(GUI_Values['nx'])
        swit['nz']  = int(GUI_Values['nz'])
        swit['dx']  = float(GUI_Values['dx'])
        swit['dt']  = float(GUI_Values['dt'])
        swit['nt']  = int(GUI_Values['nt'])
        swit['pml'] = int(GUI_Values['pml'])
        if GUI_Values['fs'] in ['Yes']:
            swit['fs'] = True
        else:
            swit['fs'] = False

        swit['vp_true_path'] = GUI_Values['vp_true_path']

        # source setup
        swit['f0']    = float(GUI_Values['f0'])
        swit['src_coor'] = GUI_Values['src_coor']
        swit['wavelet'] = GUI_Values['wavelet']
        swit['wavelet_file'] = GUI_Values['wavelet_file']

        # receiver setup
        swit['rec_coor'] = GUI_Values['rec_coor']

        return swit

    except: 
        print('Wrong: check the input parameters')



# push button and do something
def prepare_inversion_parameter(GUI_Values):
    ''' prepare parameters for inversion 
    ''' 
    try:
        # make necessary adjustment
        if GUI_Values['homepath'][-1] not in ['/']:
            GUI_Values['homepath'] += '/'
        
        swit = {}
        # system setup
        swit['homepath'] = GUI_Values['homepath']
        swit['mpiproc']  = int(GUI_Values['mpiproc'])
        swit['figaspect'] = GUI_Values['figaspect']
        # model setup
        swit['nx']  = int(GUI_Values['nx'])
        swit['nz']  = int(GUI_Values['nz'])
        swit['dx']  = float(GUI_Values['dx'])
        swit['dt']  = float(GUI_Values['dt'])
        swit['nt']  = int(GUI_Values['nt'])
        swit['pml'] = int(GUI_Values['pml'])
        if GUI_Values['fs'] in ['Yes']:
            swit['fs'] = True
        else:
            swit['fs'] = False

        swit['vp_init_path'] = GUI_Values['vp_init_path']
        swit['field_data_path'] = GUI_Values['field_data_path']

        # source setup
        swit['f0']    = float(GUI_Values['f0'])
        swit['src_coor'] = GUI_Values['src_coor']
        swit['wavelet'] = GUI_Values['wavelet']
        swit['wavelet_file'] = GUI_Values['wavelet_file']

        # receiver setup
        swit['rec_coor'] = GUI_Values['rec_coor']
    
        # basic inversion parameter
        swit['misfit_type'] = GUI_Values['misfit_type']
        swit['scheme']      = GUI_Values['scheme']
        swit['maxiter']     = int(GUI_Values['maxiter'])
        swit['step_length'] = float(GUI_Values['step_length'])
        swit['vpmax'] = float(GUI_Values['vpmax'])
        swit['vpmin'] = float(GUI_Values['vpmin'])
        swit['marine_or_land'] = GUI_Values['marine_or_land']
        
        # gradient preconditioning
        swit['grad_mute']   = int(GUI_Values['grad_mute'])
        swit['grad_smooth'] = int(GUI_Values['grad_smooth'])

        # data filter
        swit['fre_filter']  = GUI_Values['fre_filter']
        if swit['fre_filter']  in ['None']:
            swit['fre_low']  = 0.0
            swit['fre_high'] = 200.0
        elif swit['fre_filter']  in ['Bandpass']:
            swit['fre_low']  = float(GUI_Values['fre_low'])
            swit['fre_high'] = float(GUI_Values['fre_high'])
        elif swit['fre_filter']  in ['Lowpass']:
            swit['fre_low']  = float(GUI_Values['fre_low'])
            swit['fre_high'] = 200.0
        elif swit['fre_filter']  in ['Highpass']:
            swit['fre_low']  = 0.0
            swit['fre_high'] = float(GUI_Values['fre_low'])

        # pick first break and mute
        if GUI_Values['mute_late_arrival'] in ['Yes']:
            swit['mute_late_arrival'] = True
            swit['mute_late_window'] = float(GUI_Values['mute_late_window'])
        else:
            swit['mute_late_arrival'] = False
            swit['mute_late_window'] = 1000.

        # mute according to provided parameter
        if GUI_Values['mute_offset_short'] in ['Yes']:
            swit['mute_offset_short'] = True
            swit['mute_offset_short_dis'] = float(GUI_Values['mute_offset_short_dis'])
        else:
            swit['mute_offset_short'] = False
            swit['mute_offset_short_dis'] = 0.

        if GUI_Values['mute_offset_long'] in ['Yes']:
            swit['mute_offset_long'] = True
            swit['mute_offset_long_dis'] = float(GUI_Values['mute_offset_long_dis'])
        else:
            swit['mute_offset_long'] = False
            swit['mute_offset_long_dis'] = 1.e8

        swit['normalize'] = GUI_Values['normalize']

        return swit

    except: 
        print('Wrong: check the input parameters')



def View2D(filename, nx, nz):
    ''' View 2D
    '''
    try:
        nx = int(nx)
        nz = int(nz)

        if filename.endswith('.bin'):
            fp = open(filename, 'rb')
            data = np.fromfile(fp, dtype=np.float32)
            fp.close()
            data = data.reshape((nx,nz))
        elif filename.endswith(('.dat', '.txt')):
            data = np.loadtxt(filename)
        else:
            return 

        print('Viewing file:' + filename)
        print('Data max: %.2f, Data min: %.2f' %(data.max(), data.min()))
        print('Data size: %d %d' %(np.size(data,0), np.size(data,1)))

        fig = plt.figure(figsize=(10, 5))
        ax = fig.add_subplot(111)
        plotopts = {
            'vmin': data.min(),
            'vmax': data.max(),
            'cmap': plt.cm.jet,
            'aspect': 1.0
        }
        im = ax.imshow(data.T, **plotopts)
        fig.colorbar(im, shrink=0.5, extend='both')
        plt.xlabel('X grid')
        plt.ylabel('Z grid')
        try:
            plt.savefig('./View2D.png', dpi=300)
        except:
            print('wrong figure path')
        plt.close()
    except:
        print('Error: view2D goes wrong!')


def Smooth2D(smooth_filename, smooth_span, smooth_top_mute):
    ''' smooth 2D
    '''
    try:
        smooth_span = int(smooth_span)
        smooth_top_mute = int(smooth_top_mute)

        data = np.loadtxt(smooth_filename)
        data[:, smooth_top_mute:] = smooth2d(data[:,smooth_top_mute:], span=smooth_span)

        print('Smooth file: %s'%(smooth_filename))
        print('Data max: %.2f, Data min: %.2f' %(data.max(), data.min()))
        print('Data size: %d %d' %(np.size(data,0), np.size(data,1)))
        print('Span: %d  Top mute: %d '%(smooth_span, smooth_top_mute))

        smooth_filename = smooth_filename[0:-4] + '_smooth_%d'%smooth_span + smooth_filename[-4:]
        np.savetxt(smooth_filename, data)

        print('Save as: %s  '%(smooth_filename)) 
    except:
        print('Error: smooth2D goes wrong!')

# show the figure
def convert_to_bytes(file_or_bytes, resize=None):
    '''
    Will convert into bytes and optionally resize an image that is a file or a base64 bytes object.
    Turns into PNG format in the process so that can be displayed by tkinter
    '''
    if isinstance(file_or_bytes, str):
        img = PIL.Image.open(file_or_bytes)
    else:
        try:
            img = PIL.Image.open(io.BytesIO(base64.b64decode(file_or_bytes)))
        except Exception as e:
            dataBytesIO = io.BytesIO(file_or_bytes)
            img = PIL.Image.open(dataBytesIO)

    cur_width, cur_height = img.size
    if resize:
        new_width, new_height = resize
        scale = min(new_height/cur_height, new_width/cur_width)
        img = img.resize((int(cur_width*scale), int(cur_height*scale)), PIL.Image.ANTIALIAS)
    with io.BytesIO() as bio:
        img.save(bio, format="PNG")
        del img
        return bio.getvalue()

# system status
GRAPH_WIDTH, GRAPH_HEIGHT = 180, 100       # each individual graph size in pixels
class DashGraph(object):
    ''' Dash graph
    '''
    def __init__(self, graph_elem, starting_count, color):
        self.graph_current_item = 0
        self.graph_elem = graph_elem            # type:sg.Graph
        self.prev_value = starting_count
        self.max_sent = 1
        self.color = color
        self.graph_lines = []

    def graph_value(self, current_value):
        delta = current_value - self.prev_value
        self.prev_value = current_value
        self.max_sent = max(self.max_sent, delta)
        percent_sent = 100 * delta / self.max_sent
        
        line_id = self.graph_elem.draw_line((self.graph_current_item, 0), (self.graph_current_item, percent_sent), color=self.color)
        self.graph_lines.append(line_id)
        if self.graph_current_item >= GRAPH_WIDTH:
            self.graph_elem.delete_figure(self.graph_lines.pop(0))
            self.graph_elem.move(-1, 0)
        else:
            self.graph_current_item += 1
        return delta

    def graph_percentage_abs(self, value):
        self.graph_elem.draw_line((self.graph_current_item, 0), (self.graph_current_item, value), color=self.color)
        if self.graph_current_item >= GRAPH_WIDTH:
            self.graph_elem.move(-1, 0)
        else:
            self.graph_current_item += 1


def human_size(bytes, units=(' bytes', 'KB', 'MB', 'GB', 'TB', 'PB', 'EB')):
    '''
        Returns a human readable string reprentation of bytes
    '''
    return str(bytes) + units[0] if bytes < 1024 else human_size(bytes >> 10, units[1:])


def GraphColumn(name, key):
    ''' Graph column
    '''
    layout = [
        [sg.Text(name, size=(18,1), font=('newspaper 12'), key=key+'TXT_')],
        [sg.Graph((GRAPH_WIDTH, GRAPH_HEIGHT),
                    (0, 0),
                    (GRAPH_WIDTH, 20),
                    background_color='black',
                    key=key+'GRAPH_')]]
    return sg.Col(layout, pad=(2, 2))


