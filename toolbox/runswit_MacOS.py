###############################################################################
#
# SWIT: Seismic Waveform Inversion Toolbox
#
# by Haipeng Li at USTC, haipengl@mail.ustc.edu.cn
# 
# June, 1821  
#
# GUI
#
###############################################################################


import os
from tools import loadbinfloat32
import numpy as np

import psutil
import PySimpleGUI as sg

from guitools import (DashGraph, GraphColumn, Smooth2D, View2D,
                      convert_to_bytes, forward_workflow, human_size,
                      inversion_workflow, prepare_forward_parameter,
                      prepare_inversion_parameter)

########################################
#
# GUI begins here
#
########################################


def main():

    # signs and colors

    BORDER_COLOR = 'grey99'
    BPAD_LEFT = ((10,10), (10, 10))
    BPAD_LEFT_INSIDE = (2, 2)

    # theme
    sg.theme('Reddit')  # 'SystemDefault', 'LightGrey4'

    # set_options
    sg.set_options(border_width=1, margins=(2, 2), element_padding=(0, 0))

    font0 = ("courier 10 pitch", 14)
    font1 = ("courier 10 pitch", 14)
    font2 = ("courier 10 pitch", 16) # 'bold'
    font3 = ("courier 10 pitch", 18)

    size2 = (18, 1)
    size2 = (28, 1)
    size3 = (18, 1)
    size4 = (10, 1)
    size_txt = (18, 1)

    # Menu definition
    menu_def = [['&Information', ['&By Haipeng Li at USTC', '&Contact: haipengl@mail.ustc.edu.cn', '&Copyright © 2021 All Right Reserved', ]] ]

    # Tab layouts
    tab_sys_layout = [
                    [sg.Text('        ',      font=font1)],
                    [sg.Text('Homepath',      font=font1)],
                    [sg.Input(size=size2,     font=font1, enable_events=True, key='homepath'), sg.FolderBrowse(initial_folder='~/', font=font1)],
                    [sg.Text('CPU number',    font=font1)], 
                    [sg.Input(size=size2,     font=font1, key='mpiproc')],
                    [sg.Text('Figure aspect', font=font1)], 
                    [sg.Input(size=size2,     font=font1, key='figaspect')],
                     ]

    tab_mod_layout = [
                    [sg.Text('  ',  size=size_txt, font=font1)],
                    [sg.Text('Nx',  size=size_txt, font=font1), sg.Input(size=size3, key='nx', font=font1)],
                    [sg.Text('Nz',  size=size_txt, font=font1), sg.Input(size=size3, key='nz', font=font1)],
                    [sg.Text('Dx',  size=size_txt, font=font1), sg.Input(size=size3, key='dx', font=font1)],
                    [sg.Text('Dt',  size=size_txt, font=font1), sg.Input(size=size3, key='dt', font=font1)],
                    [sg.Text('Nt',  size=size_txt, font=font1), sg.Input(size=size3, key='nt', font=font1)],
                    [sg.Text('PML', size=size_txt, font=font1), sg.Input(size=size3, key='pml',font=font1)],
                    [sg.Text('Free surface', size=size_txt, font=font1), sg.Combo(('Yes', 'No'), default_value='Yes', size=size3, font=font1, key='fs')],
                    
                    [sg.Text(' ',                 size=size2, font=font1)],
                    [sg.Text('Vp true (forward)', size=size2, font=font1 )], 
                    [sg.Input(size=size2, enable_events=True, key='vp_true_path',    font=font1), sg.FileBrowse(initial_folder='~/', font=font1, file_types=(("dat Files", "*.dat"), ("Bin Files", "*.bin")))],
                    [sg.Text('Vp init (inversion)', size=size2, font=font1 )], 
                    [sg.Input(size=size2, enable_events=True, key='vp_init_path',    font=font1), sg.FileBrowse(initial_folder='~/', font=font1, file_types=(("dat Files", "*.dat"), ("Bin Files", "*.bin")))],
                    [sg.Text('SU data (inversion)', size=size2, font=font1 )], 
                    [sg.Input(size=size2, enable_events=True ,key='field_data_path', font=font1), sg.FolderBrowse(initial_folder='~/', font=font1)],
                     ]

    tab_acq_layout = [
                    [sg.Text('              ',    font=font1)],
                    [sg.Text('Receiver coordinate',     size=size2, font=font1)], 
                    [sg.In(size=size2, enable_events=True,key='rec_coor',  font=font1), sg.FileBrowse(initial_folder='~/', font=font1, file_types=(("dat Files", "*.dat"), ("Txt Files", "*.Txt")))],
                    [sg.Text('Source coordinate',       size=size2, font=font1)], 
                    [sg.In(size=size2, enable_events=True, key='src_coor', font=font1), sg.FileBrowse(initial_folder='~/', font=font1, file_types=(("dat Files", "*.dat"), ("Txt Files", "*.Txt")))],
                    [sg.Text('              ', size=size_txt, font=font1)],
                    [sg.Text('Land or Marine', size=size_txt, font=font1), sg.Combo(('Land', 'Marine'),                   default_value='Land',   key='marine_or_land', font=font1, size=size3)],
                    [sg.Text('Source wavelet', size=size_txt, font=font1), sg.Combo(('Ricker (Opt 1)', 'File   (Opt 2)'), default_value='Ricker', key='wavelet',        font=font1, size=size3)],
                    [sg.Text('Opt 1: F0 (Hz)', size=size_txt, font=font1), sg.Input(size=size3, key='f0', font=font1)],
                    [sg.Text('Opt 2: File   ', size=size_txt, font=font1)], 
                    [sg.In(size=size2, enable_events=True,key='wavelet_file', font=font1), sg.FileBrowse(initial_folder='~/', font=font1, file_types=(("dat Files", "*.dat")))],
                ]

    tab_inv_layout =[
                    [sg.Text('              ',    font=font1)],
                    [sg.Text('Misfit',            size=size_txt, font=font1), sg.Combo(('Waveform', 'Traveltime', 'Envelope', 'Globalcorrelation'), default_value='Waveform', key='misfit_type', font=font1, size=size3)],
                    [sg.Text('Scheme',            size=size_txt, font=font1), sg.Combo(('NLCG','LBFGS'), default_value='NLCG', key='scheme', font=font1, size=size3)],
                    [sg.Text('Step length',       size=size_txt, font=font1), sg.Input(size=size3, key='step_length', font=font1)],
                    [sg.Text('Iteration',         size=size_txt, font=font1), sg.Input(size=size3, key='maxiter',     font=font1)],
                    [sg.Text('Vpmax',             size=size_txt, font=font1), sg.Input(size=size3, key='vpmax',       font=font1)],
                    [sg.Text('Vpmin',             size=size_txt, font=font1), sg.Input(size=size3, key='vpmin',       font=font1)],
                    [sg.Text('Gradient mute',     size=size_txt, font=font1), sg.Input(size=size3, key='grad_mute',   font=font1)],
                    [sg.Text('Gradient smooth',   size=size_txt, font=font1), sg.Input(size=size3, key='grad_smooth', font=font1)],
                    [sg.Text('Normalization',     size=size_txt, font=font1), sg.Combo(('L1-Trace', 'L2-Trace', 'L1-Event', 'L2-Event', 'None'), default_value='L1-Trace', key='normalize', font=font1, size=size3)],
                    [sg.Text('Frequency filter',  size=size_txt, font=font1), sg.Combo(('None', 'Bandpass', 'Lowpass', 'Highpass'), default_value='None', key='fre_filter', font=font1, size=size3)],
                    [sg.Text('Frequency low',     size=size_txt, font=font1), sg.Input(size=size3, key='fre_low',  font=font1)],
                    [sg.Text('Frequency high',    size=size_txt, font=font1), sg.Input(size=size3, key='fre_high', font=font1)],
                    [sg.Text('Mute late arrival', size=size_txt, font=font1), sg.Combo(('Yes', 'No'), default_value='No', size=size3, font=font1, key='mute_late_arrival')],
                    [sg.Text('Mute time window',  size=size_txt, font=font1), sg.Input(size=size3, key='mute_late_window', font=font1)],
                    [sg.Text('Mute near offset',  size=size_txt, font=font1), sg.Combo(('Yes', 'No'), default_value='No', size=size3, font=font1, key='mute_offset_short')],
                    [sg.Text('Mute near distance',size=size_txt, font=font1), sg.Input(size=size3, key='mute_offset_short_dis', font=font1)],
                    [sg.Text('Mute far offset',   size=size_txt, font=font1), sg.Combo(('Yes', 'No'), default_value='No', size=size3, font=font1, key='mute_offset_long')],
                    [sg.Text('Mute far distance', size=size_txt, font=font1), sg.Input(size=size3, key='mute_offset_long_dis', font=font1)],
                ]

    # The TabgGroup layout - it must contain only Tabs
    tab_group_layout = [[sg.Tab('System',      tab_sys_layout, font=font3, key='sys'),
                         sg.Tab('Model',       tab_mod_layout, font=font3, key='mod'),
                         sg.Tab('Acquisition', tab_acq_layout, font=font3, key='acq'),
                         sg.Tab('Inversion',   tab_inv_layout, font=font3, key='inv'),
                        ]
                    ]

    # left Columns
    block_1 = [
                [sg.Text('Input Parameters', justification='left', font=font3)],
                [sg.Text('_' * 60 + '\n')],
                [sg.TabGroup(tab_group_layout, enable_events=True, key='-TABGROUP-', font=font2)],
                [sg.Text(' ' * 60 + '\n')],
                [sg.Button('Save parameters', font=font1), sg.Text(' ' * 5), sg.Button('Load parameters', font=font1)],
            ]

    block_2 = [ [sg.Text('Functions', justification='left', font=font3)],
                [sg.Text('_' * 60 + '\n')],
                [sg.In(size=(6,1), enable_events=True, key='view2D_file', font=font1),             
                sg.FileBrowse(initial_folder='~/', font=font1, file_types=(("dat Files", "*.dat"), ("Bin Files", "*.bin"))),
                sg.Text('nx', size=(2, 1), justification='left', font=font1),  sg.In(size=(6,1), enable_events=True, key='view2D_nx', font=font1),
                sg.Text('nz', size=(2, 1), justification='left', font=font1),  sg.In(size=(6,1), enable_events=True, key='view2D_nz', font=font1),
                sg.Button('View',  size=(6, 1), button_color=('white', 'black'), border_width=1, font=font1)],
                
                [sg.In(size=(6,1), enable_events=True, key='bin_file', font=font1),             
                sg.FileBrowse(initial_folder='~/', font=font1, file_types=(("Bin Files", "*.bin"), ("dat Files", "*.dat"))),
                sg.Text('nx', size=(2, 1), justification='left', font=font1),  sg.In(size=(6,1), enable_events=True, key='bin2dat_nx', font=font1),
                sg.Text('nz', size=(2, 1), justification='left', font=font1),  sg.In(size=(6,1), enable_events=True, key='bin2dat_nz', font=font1),
                sg.Button('Bin2dat',  size=(6, 1), button_color=('white', 'black'), border_width=1, font=font1)],
                

                [sg.In(size=(6,1), enable_events=True, key='smooth_file', font=font1),             
                sg.FileBrowse(initial_folder='~/', font=font1, file_types=(("dat Files", "*.dat"), ("Bin Files", "*.bin"))),
                sg.Text('sp', size=(2, 1), justification='left', font=font1),  sg.In(size=(6,1), enable_events=True, key='smooth_span',     font=font1),
                sg.Text('mt', size=(2, 1), justification='left', font=font1),  sg.In(size=(6,1), enable_events=True, key='smooth_top_mute', font=font1),
                sg.Button('Smooth',  size=(6, 1), button_color=('white', 'black'), border_width=1, font=font1)],
                
                [sg.Text(' ' * 60 + '\n')],
                [sg.Button('  Forward  ',  size=size4, button_color=('white', 'blue'),  border_width=4, font=font2),
                 sg.Button('    FWI    ',  size=size4, button_color=('white', 'red'),   border_width=4, font=font2),
                 sg.Button('   Clear   ',  size=size4, button_color=('white', 'green'), border_width=4, font=font2)],
               ]

    # middle Columns
    block_3 = [[sg.Text('Viewing Window', justification='left', font=font3)],
               [sg.Text('_' * 300 + '\n')],
               [sg.Image(key='figure')]
            ]

    block_4 = [ [sg.Text('Output Histories', justification='left', font=font3)],
                [sg.Text('_' * 300 + '\n')],
                [sg.Output(size=(85, 13), key = 'output', font=font0)],
            ]

    # right Columns
    block_5 = [[sg.Text('Viewing Options', justification='left', font=font3)],
            [sg.Text('_' * 300 + '\n')],
            [sg.Combo(('Acquisition', 'Wavelet', 'Waveform', 'Velocity', 'Gradient', 'Direction', 'Illumination'), default_value='Acquisition', enable_events=True, key='fig_type', font=font1, size=(36, 8))],
            [sg.Text(' ' * 300 + '\n')],
            [sg.Listbox(values=[], enable_events=True, size=(36,26),key='fig_list', font=font0)],
            ]

    block_6 = [[sg.Text('System Status', justification='left', font=font3)],
            [sg.Text('_' * 300 + '\n')],
            [GraphColumn('Disk Read',    '_DISK_READ_'),
             GraphColumn('Disk Write',   '_DISK_WRITE_')],
            [GraphColumn('CPU Usage',    '_CPU_'),
             GraphColumn('Memory Usage', '_MEM_')], 
            ]


    # overall layout
    layout = [[sg.Menu(menu_def, tearoff=True)],
            [
                sg.Column([ [sg.Column(block_1, size=(350,600), pad=BPAD_LEFT_INSIDE)],
                            [sg.Column(block_2, size=(350,300), pad=BPAD_LEFT_INSIDE)]], 
                            pad=BPAD_LEFT, background_color=BORDER_COLOR, vertical_alignment='center', justification='center',  k='C1'),
                sg.VSeperator(),
                sg.Column([ [sg.Column(block_3, size=(800,600), pad=BPAD_LEFT_INSIDE)],
                            [sg.Column(block_4, size=(800,300), pad=BPAD_LEFT_INSIDE)]], 
                            pad=BPAD_LEFT, background_color=BORDER_COLOR, vertical_alignment='center', justification='center',  k='C2'),
                sg.VSeperator(),
                sg.Column([ [sg.Column(block_5, size=(350,600), pad=BPAD_LEFT_INSIDE)],
                            [sg.Column(block_6, size=(350,300), pad=BPAD_LEFT_INSIDE)]], 
                            pad=BPAD_LEFT, background_color=BORDER_COLOR, vertical_alignment='center', justification='center',  k='C3'),
            ]
            ]


    # show window
    window = sg.Window('Seismic Waveform Inversion Toolbox-1.0', layout, resizable=True, finalize=True) #  icon='./SWIT.jpg'
    # keep in the center
    window['C1'].expand(True, True, True)
    window['C2'].expand(True, True, True)
    window['C3'].expand(True, True, True)

    # setup graphs & initial values
    diskio = psutil.disk_io_counters()
    disk_graph_write = DashGraph(window['_DISK_WRITE_GRAPH_'], diskio.write_bytes, '#be45be')
    disk_graph_read  = DashGraph(window['_DISK_READ_GRAPH_'], diskio.read_bytes,   '#5681d8')
    cpu_usage_graph  = DashGraph(window['_CPU_GRAPH_'], 0, '#d34545')
    mem_usage_graph  = DashGraph(window['_MEM_GRAPH_'], 0, '#BE7C29')

    view_file_is_ok = 0
    view_nx_is_ok = 0
    view_nz_is_ok = 0
    bin_file_is_ok = 0
    bin2dat_nx_is_ok = 0
    bin2dat_nz_is_ok = 0
    smooth_file_is_ok = 0
    smooth_span_is_ok = 0
    smooth_top_mute_is_ok = 0
    
    new_size = (1000,470)

    # Event Loop
    while True:
        event, values = window.read(timeout=500)
        if event == sg.WIN_CLOSED or event == 'Exit':
            break

        elif event == 'Save parameters':
            filename = sg.popup_get_file('Save parameters', save_as=True, no_window=True)
            window.SaveToDisk(filename)

        elif event == 'Load parameters':
            filename = sg.popup_get_file('Load parameters', initial_folder='~/', no_window=True)
            window.LoadFromDisk(filename)
            window.FindElement('output').Update('')

        elif event == 'view2D_file':
            view_file_is_ok = 1

        elif event == 'view2D_nx':
            view_nx_is_ok = 1

        elif event == 'view2D_nz':
            view_nz_is_ok = 1

        elif event == 'View' and view_file_is_ok:
            view_filename = values['view2D_file']
            if view_nx_is_ok and view_nz_is_ok:
                view_nx = values['view2D_nx']
                view_nz = values['view2D_nz']
            else:
                view_nx = 1
                view_nz = 1
            View2D(view_filename, view_nx, view_nz)
            window['figure'].update(data=convert_to_bytes('./View2D.png', resize=new_size))


        elif event == 'bin_file':
            bin_file_is_ok = 1

        elif event == 'bin2dat_nx':
            bin2dat_nx_is_ok = 1

        elif event == 'bin2dat_nz':
            bin2dat_nz_is_ok = 1

        elif event == 'Bin2dat' and bin_file_is_ok and bin2dat_nx_is_ok and bin2dat_nz_is_ok:
            try:
                file = values['bin_file']
                nx = int(values['bin2dat_nx'])
                nz = int(values['bin2dat_nz'])
                data = loadbinfloat32(file)
                filename = file[:-4] + '.dat'
                np.savetxt(filename, data.reshape(nx , nz))
                print('Save as: %s'%(filename))
            except:
                file = values['bin_file']
                data = loadbinfloat32(file)
                size = data.shape[0]
                print('nx or nz is inconsistent with the bin file: %s' % (values['bin_file']))
                print('the length of the bin file is: %d of float32'%size)


        elif event == 'smooth_file':
            smooth_file_is_ok = 1

        elif event == 'smooth_span':
            smooth_span_is_ok = 1

        elif event == 'smooth_top_mute':
            smooth_top_mute_is_ok = 1

        elif event == 'Smooth' and smooth_file_is_ok and smooth_span_is_ok and smooth_top_mute_is_ok:
            smooth_filename = values['smooth_file']
            smooth_span     = values['smooth_span']
            smooth_top_mute = values['smooth_top_mute']
            Smooth2D(smooth_filename, smooth_span, smooth_top_mute)

        elif event == '  Forward  ':
            swit = prepare_forward_parameter(values)
            forward_workflow(swit)

        elif event == '    FWI    ':
            swit = prepare_inversion_parameter(values)
            inversion_workflow(swit)

        elif event == '   Clear   ':
            window.FindElement('output').Update('')

        elif event == 'fig_type':
            folder = values['homepath']  + '/figures/'
            fig_type = values['fig_type']
            print('Viewing %s'%(fig_type))

            if fig_type in ['Acquisition','Wavelet']:
                suffix = ''
            elif fig_type in ['Waveform']:
                suffix = 'waveform/'
            elif fig_type in ['Velocity', 'Gradient', 'Direction', 'Illumination']:
                suffix = 'model/'
            else:
                suffix = ''
                print('No support figure type\n')
            try:
                file_list = os.listdir(folder+suffix)
            except:
                file_list = []

            fnames = [f for f in file_list if os.path.isfile(os.path.join(folder+suffix, f)) and f.lower().endswith((".png", ".jpg", "jpeg", ".tiff", ".bmp"))]
            fnames = sorted(fnames)
            # select figures for different view
            if   fig_type in ['Velocity']:
                fnames = [f for f in fnames if f[0:2] == 'vp']
            elif fig_type in ['Gradient']:
                fnames = [f for f in fnames if f[0:4] == 'grad']
            elif fig_type in ['Direction']:
                fnames = [f for f in fnames if f[0:4] == 'dire']
            elif fig_type in ['Illumination']:
                fnames = [f for f in fnames if f[5:10] == 'illum']
            else:
                pass

            window['fig_list'].update(fnames)

        elif event == 'fig_list':    # A file was chosen from the listbox
            try:
                filename = os.path.join(values['homepath']  + '/figures/' + suffix, values['fig_list'][0])
                window['figure'].update(data=convert_to_bytes(filename, resize=new_size))
            except Exception as E:
                print(f'** Error {E} **')
                pass                 # something weird happened making the full filename
        
        # update
        # Disk Graphs 
        diskio = psutil.disk_io_counters()
        write_bytes = disk_graph_write.graph_value(diskio.write_bytes)
        read_bytes = disk_graph_read.graph_value(diskio.read_bytes)
        window['_DISK_WRITE_TXT_'].update('Disk Write {}'.format(human_size(write_bytes)))
        window['_DISK_READ_TXT_'].update('Disk Read {}'.format(human_size(read_bytes)))
        # CPU Graph
        cpu = psutil.cpu_percent(0)
        cpu_usage_graph.graph_percentage_abs(cpu)
        window['_CPU_TXT_'].update('{0:2.0f}% CPU Used'.format(cpu))
        # Memory Graph
        mem_used = psutil.virtual_memory().percent
        mem_usage_graph.graph_percentage_abs(mem_used)
        window['_MEM_TXT_'].update('{}% Memory Used'.format(mem_used))

        
    window.close()


if __name__ == '__main__':
    
    main()

