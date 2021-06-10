###############################################################################
#
# SWIT: Seismic Waveform Inversion Toolbox
#
# by Haipeng Li at USTC, haipengl@mail.ustc.edu.cn
# 
# June, 2021  
#
# GUI
#
###############################################################################


import os
import psutil
import PySimpleGUI as sg


from guitools import DashGraph, GraphColumn, View2D, convert_to_bytes, human_size
from guitools import Smooth2D, inversion_workflow, forward_workflow
from guitools import prepare_forward_parameter, prepare_inversion_parameter


########################################
#
# GUI begins here
#
########################################


def main():

    # signs and colors

    BORDER_COLOR = '#C7D5E0'
    BPAD_LEFT = ((0,0), (0, 0))
    BPAD_LEFT_INSIDE = (0, 0)

    # theme
    sg.theme('Reddit')  # 'SystemDefault', 'LightGrey4'

    # set_options
    sg.set_options(border_width=1, margins=(0, 0), element_padding=(0, 0))

    # Menu definition
    menu_def = [['&Information', ['&By Haipeng Li at USTC', '&Contact: haipengl@mail.ustc.edu.cn', '&Copyright © 2021 All Right Reserved', ]] ]

    # Tab layouts
    tab_sys_layout = [
                    [sg.Text('              ',   font=("Any", 16))],
                    [sg.Text('Homepath',         font=("Any", 16))],
                    [sg.Input(size=(20,1), enable_events=True, key='homepath', font=("Any", 16)), sg.FolderBrowse(initial_folder='~/Desktop/', font=("Any", 16))],
                    [sg.Text('CPU number', font=("Any", 16))], 
                    [sg.Input(size=(20,1), key='mpiproc', font=("Any", 16))],
                    [sg.Text('Figure aspect', font=("Any", 16))], 
                    [sg.Input(size=(20,1), key='figaspect', font=("Any", 16))],
                ]

    tab_mod_layout =[
                    [sg.Text('              ',    font=("Any", 16))],
                    [sg.Text('Nx',  size=(16, 1), font=("Any", 16)), sg.Input(size=(16,1), key='nx', font=("Any", 16))],
                    [sg.Text('Nz',  size=(16, 1), font=("Any", 16)), sg.Input(size=(16,1), key='nz', font=("Any", 16))],
                    [sg.Text('Dx',  size=(16, 1), font=("Any", 16)), sg.Input(size=(16,1), key='dx', font=("Any", 16))],
                    [sg.Text('Dt',  size=(16, 1), font=("Any", 16)), sg.Input(size=(16,1), key='dt', font=("Any", 16))],
                    [sg.Text('Nt',  size=(16, 1), font=("Any", 16)), sg.Input(size=(16,1), key='nt', font=("Any", 16))],
                    [sg.Text('PML', size=(16, 1), font=("Any", 16)), sg.Input(size=(16,1), key='pml', font=("Any", 16))],
                    [sg.Text('Free surface',      size=(16, 1), font=("Any", 16)), sg.Combo(('Yes', 'No'), default_value='Yes', size=(15, 1), font=("Any", 16), key='fs')],
                    [sg.Text(' ',                 size=(16, 1), font=("Any", 16))],
                    [sg.Text('Vp true (forward)', size=(24, 1), font=("Any", 16) )], 
                    [sg.Input(size=(26,1), enable_events=True, key='vp_true_path', font=("Any", 16)), sg.FileBrowse(initial_folder='~/Desktop/', font=("Any", 16), file_types=(("dat Files", "*.dat"), ("Bin Files", "*.bin")))],
                    [sg.Text('Vp init (inversion)', size=(24, 1), font=("Any", 16) )], 
                    [sg.Input(size=(26,1), enable_events=True, key='vp_init_path', font=("Any", 16)), sg.FileBrowse(initial_folder='~/Desktop/', font=("Any", 16), file_types=(("dat Files", "*.dat"), ("Bin Files", "*.bin")))],
                    [sg.Text('SU data (inversion)', size=(24, 1), font=("Any", 16) )], 
                    [sg.Input(size=(26,1), enable_events=True ,key='field_data_path', font=("Any", 16)), sg.FolderBrowse(initial_folder='~/Desktop/', font=("Any", 16))],
                ]

    tab_acq_layout = [
                    [sg.Text('              ',    font=("Any", 16))],
                    [sg.Text('Receiver coordinate',     size=(16, 1), font=("Any", 16))], 
                    [sg.In(size=(26,1), enable_events=True,key='rec_coor', font=("Any", 16)), sg.FileBrowse(initial_folder='~/Desktop/', font=("Any", 16), file_types=(("dat Files", "*.dat")))],
                    [sg.Text('Source coordinate',       size=(16, 1), font=("Any", 16))], 
                    [sg.In(size=(26,1), enable_events=True, key='src_coor', font=("Any", 16)), sg.FileBrowse(initial_folder='~/Desktop/', font=("Any", 16), file_types=(("dat Files", "*.dat")))],
                    [sg.Text('              ',          size=(16, 1), font=("Any", 16))],
                    [sg.Text('Land or Marine',          size=(16, 1), font=("Any", 16)), sg.Combo(('Land', 'Marine'), default_value='Land', key='marine_or_land', font=("Any", 16), size=(15, 1))],
                    [sg.Text('Source wavelet',          size=(16, 1), font=("Any", 16)), sg.Combo(('Ricker (Opt 1)', 'File   (Opt 2)'), default_value='Ricker', key='wavelet', font=("Any", 16), size=(15, 1))],
                    [sg.Text('Opt 1: F0 (Hz)',          size=(16, 1), font=("Any", 16)), sg.Input(size=(16,1), key='f0', font=("Any", 16))],
                    [sg.Text('Opt 2: File   ',          size=(16, 1), font=("Any", 16))], 
                    [sg.In(size=(26,1), enable_events=True,key='wavelet_file', font=("Any", 16)), sg.FileBrowse(initial_folder='~/Desktop/', font=("Any", 16), file_types=(("dat Files", "*.dat")))],
                ]

    tab_inv_layout =[
                    [sg.Text('              ',    font=("Any", 16))],
                    [sg.Text('Misfit',            size=(20, 1), font=("Any", 16)), sg.Combo(('Waveform', 'Traveltime', 'Envelope', 'Globalcorrelation'), default_value='Waveform', key='misfit_type', font=("Any", 16), size=(12, 1))],
                    [sg.Text('Scheme',            size=(20, 1), font=("Any", 16)), sg.Combo(('NLCG','LBFGS'), default_value='NLCG', key='scheme', font=("Any", 16), size=(12, 1))],
                    [sg.Text('Step length',       size=(20, 1), font=("Any", 16)), sg.Input(size=(13,1), key='step_length', font=("Any", 16))],
                    [sg.Text('Iteration',         size=(20, 1), font=("Any", 16)), sg.Input(size=(13,1), key='maxiter', font=("Any", 16))],
                    [sg.Text('Vpmax',             size=(20, 1), font=("Any", 16)), sg.Input(size=(13,1), key='vpmax', font=("Any", 16))],
                    [sg.Text('Vpmin',             size=(20, 1), font=("Any", 16)), sg.Input(size=(13,1), key='vpmin', font=("Any", 16))],
                    [sg.Text('Gradient mute',     size=(20, 1), font=("Any", 16)), sg.Input(size=(13,1), key='grad_mute', font=("Any", 16))],
                    [sg.Text('Gradient smooth',   size=(20, 1), font=("Any", 16)), sg.Input(size=(13,1), key='grad_smooth', font=("Any", 16))],
                    [sg.Text('Normalization',     size=(20, 1), font=("Any", 16)), sg.Combo(('L1-Trace', 'L2-Trace', 'L1-Event', 'L2-Event', 'None'), default_value='L1-Trace', key='normalize', font=("Any", 16), size=(12, 1))],
                    [sg.Text('Frequency filter',  size=(20, 1), font=("Any", 16)), sg.Combo(('None', 'Bandpass', 'Lowpass', 'Highpass'), default_value='None', key='fre_filter', font=("Any", 16), size=(12, 1))],
                    [sg.Text('Frequency low',     size=(20, 1), font=("Any", 16)), sg.Input(size=(13,1), key='fre_low', font=("Any", 16))],
                    [sg.Text('Frequency high',    size=(20, 1), font=("Any", 16)), sg.Input(size=(13,1), key='fre_high', font=("Any", 16))],
                    [sg.Text('Mute late arrival', size=(20, 1), font=("Any", 16)), sg.Combo(('Yes', 'No'), default_value='No', size=(12, 1), font=("Any", 16), key='mute_late_arrival')],
                    [sg.Text('Mute time window',  size=(20, 1), font=("Any", 16)), sg.Input(size=(13,1), key='mute_late_window', font=("Any", 16))],
                    [sg.Text('Mute near offset',  size=(20, 1), font=("Any", 16)), sg.Combo(('Yes', 'No'), default_value='No', size=(12, 1), font=("Any", 16), key='mute_offset_short')],
                    [sg.Text('Mute near distance',size=(20, 1), font=("Any", 16)), sg.Input(size=(13,1), key='mute_offset_short_dis', font=("Any", 16))],
                    [sg.Text('Mute far offset',   size=(20, 1), font=("Any", 16)), sg.Combo(('Yes', 'No'), default_value='No', size=(12, 1), font=("Any", 16), key='mute_offset_long')],
                    [sg.Text('Mute far distance', size=(20, 1), font=("Any", 16)), sg.Input(size=(13,1), key='mute_offset_long_dis', font=("Any", 16))],
                ]

    # The TabgGroup layout - it must contain only Tabs
    tab_group_layout = [[sg.Tab('System',      tab_sys_layout, font=("Any", 20), key='sys'),
                         sg.Tab('Model',       tab_mod_layout, font=("Any", 20), key='mod'),
                         sg.Tab('Acquisition', tab_acq_layout, font=("Any", 20), key='acq'),
                         sg.Tab('Inversion',   tab_inv_layout, font=("Any", 20), key='inv'),
                        ]
                    ]

    # left Columns
    block_1 = [
                [sg.Text('Input Parameters', justification='left', font=("Any", 20))],
                [sg.Text('_' * 60 + '\n')],
                [sg.TabGroup(tab_group_layout, enable_events=True, key='-TABGROUP-', font=("Any", 18))],
                [sg.Text(' ' * 60 + '\n')],
                [sg.Button('Save parameters', font=("Any", 16)), sg.Text(' ' * 5), sg.Button('Load parameters', font=("Any", 16))],
            ]

    block_2 = [ [sg.Text('Functions', justification='left', font=("Any", 20))],
                [sg.Text('_' * 60 + '\n')],
                [sg.In(size=(10,1), enable_events=True, key='view2D_file'),             
                sg.FileBrowse(initial_folder='~/Desktop/', font=("Any", 16), file_types=(("dat Files", "*.dat"), ("Bin Files", "*.bin"))),
                sg.Text('nx', size=(2, 1), justification='left', font=("Any", 16)),
                sg.In(size=(6,1), enable_events=True, key='view2D_nx'),
                sg.Text('nz', size=(2, 1), justification='left', font=("Any", 16)),
                sg.In(size=(6,1), enable_events=True, key='view2D_nz'),
                sg.Button('View',  size=(6, 1), button_color=('white', 'black'), border_width=1, font='Any 16')],
                
                [sg.In(size=(10,1), enable_events=True, key='smooth_file'),             
                sg.FileBrowse(initial_folder='~/Desktop/', font=("Any", 16), file_types=(("dat Files", "*.dat"), ("Bin Files", "*.bin"))),
                sg.Text('sp', size=(2, 1), justification='left', font=("Any", 16)),
                sg.In(size=(6,1), enable_events=True, key='smooth_span'),
                sg.Text('mt', size=(2, 1), justification='left', font=("Any", 16)),
                sg.In(size=(6,1), enable_events=True, key='smooth_top_mute'),
                sg.Button('Smooth',  size=(6, 1), button_color=('white', 'black'), border_width=1, font='Any 16')],
                
                [sg.Text(' ' * 60 + '\n')],
                [sg.Button('  Forward  ',  size=(10, 1), button_color=('white', 'blue'), border_width=4, font='Any 18'),
                 sg.Button('    FWI    ',  size=(10, 1), button_color=('white', 'red'), border_width=4,  font='Any 18'),
                 sg.Button('   Clear   ',  size=(10, 1), button_color=('white', 'green'), border_width=4, font='Any 18')],
            ]

    # middle Columns
    block_3 = [[sg.Text('Viewing Window', justification='left', font=("Any", 20))],
               [sg.Text('_' * 300 + '\n')],
               [sg.Image(key='figure')]
            ]

    block_4 = [ [sg.Text('Output Histories', justification='left', font=("Any", 20))],
                [sg.Text('_' * 300 + '\n')],
                [sg.Output(size=(160, 20), key = 'output', font=('Any 16'))],
            ]

    # right Columns
    block_5 = [[sg.Text('Viewing Options', justification='left', font=("Any", 20))],
            [sg.Text('_' * 300 + '\n')],
            [sg.Combo(('Acquisition', 'Source-Wavelet', 'Velocity', 'Waveform', 'Gradient', 'Direction'), default_value='Acquisition', enable_events=True, key='fig_type', font=("Any", 16), size=(26, 8))],
            [sg.Text(' ' * 300 + '\n')],
            [sg.Listbox(values=[], enable_events=True, size=(42,38),key='fig_list', font=("Any", 16))],
            ]

    block_6 = [[sg.Text('System Status', justification='left', font=("Any", 20))],
            [sg.Text('_' * 300 + '\n')],
            [GraphColumn('Disk Read',    '_DISK_READ_'),
             GraphColumn('Disk Write',   '_DISK_WRITE_')],
            [GraphColumn('CPU Usage',    '_CPU_'),
             GraphColumn('Memory Usage', '_MEM_')], 
            ]


    # overall layout
    layout = [[sg.Menu(menu_def, tearoff=True)],
            [
                sg.Column([ [sg.Column(block_1, size=(400,600), pad=BPAD_LEFT_INSIDE)],
                            [sg.Column(block_2, size=(400,300), pad=BPAD_LEFT_INSIDE)]], 
                            pad=BPAD_LEFT, background_color=BORDER_COLOR, vertical_alignment='center', justification='center',  k='C1'),
                sg.VSeperator(),
                sg.Column([ [sg.Column(block_3, size=(900,600), pad=BPAD_LEFT_INSIDE)],
                            [sg.Column(block_4, size=(900,300), pad=BPAD_LEFT_INSIDE)]], 
                            pad=BPAD_LEFT, background_color=BORDER_COLOR, vertical_alignment='center', justification='center',  k='C2'),
                sg.VSeperator(),
                sg.Column([ [sg.Column(block_5, size=(300,600), pad=BPAD_LEFT_INSIDE)],
                            [sg.Column(block_6, size=(300,300), pad=BPAD_LEFT_INSIDE)]], 
                            pad=BPAD_LEFT, background_color=BORDER_COLOR, vertical_alignment='center', justification='center',  k='C3'),
            ]
                #    [sg.HorizontalSeparator()],
            ]



    # show window

    window = sg.Window('Seismic Waveform Inversion Toolbox-1.0', layout, resizable=True, finalize=True) # , location=(0, 0), size=(1600, 800)
    # keep in the center
    #window['C1'].expand(True, True, True)
    #window['C2'].expand(True, True, True)
    #window['C3'].expand(True, True, True)

    # setup graphs & initial values
    diskio = psutil.disk_io_counters()
    disk_graph_write = DashGraph(window['_DISK_WRITE_GRAPH_'], diskio.write_bytes, '#be45be')
    disk_graph_read = DashGraph(window['_DISK_READ_GRAPH_'], diskio.read_bytes, '#5681d8')
    cpu_usage_graph = DashGraph(window['_CPU_GRAPH_'], 0, '#d34545')
    mem_usage_graph = DashGraph(window['_MEM_GRAPH_'], 0, '#BE7C29')

    view_file_is_ok = 0
    view_nx_is_ok = 0
    view_nz_is_ok = 0
    smooth_file_is_ok = 0
    smooth_span_is_ok = 0
    smooth_top_mute_is_ok = 0

    # Event Loop
    while True:
        event, values = window.read(timeout=500)
        if event == sg.WIN_CLOSED or event == 'Exit':
            break

        elif event == 'Save parameters':
            filename = sg.popup_get_file('Save parameters', save_as=True, no_window=True)
            window.SaveToDisk(filename)

        elif event == 'Load parameters':
            filename = sg.popup_get_file('Load parameters', initial_folder='~/Desktop/', no_window=True)
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
            new_size = (900,500)
            window['figure'].update(data=convert_to_bytes('./View2D.png', resize=new_size))

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

            if fig_type in ['Acquisition','Source-Wavelet']:
                suffix = ''
            elif fig_type in ['Waveform']:
                suffix = 'waveform/'
            elif fig_type in ['Velocity']:
                suffix = 'model/'
            elif fig_type in ['Gradient']:
                suffix = 'model/'
            elif fig_type in ['Direction']:
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
            if fig_type in ['Velocity']:
                fnames = [f for f in fnames if f[0:2] == 'vp']
            elif fig_type in ['Gradient']:
                fnames = [f for f in fnames if f[0:4] == 'grad']
            elif fig_type in ['Direction']:
                fnames = [f for f in fnames if f[0:4] == 'dire']
            else:
                pass

            window['fig_list'].update(fnames)


        elif event == 'fig_list':    # A file was chosen from the listbox
            try:
                filename = os.path.join(values['homepath']  + '/figures/' + suffix, values['fig_list'][0])
                new_size = (900,800)
                window['figure'].update(data=convert_to_bytes(filename, resize=new_size))
            except Exception as E:
                print(f'** Error {E} **')
                pass        # something weird happened making the full filename
        
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

