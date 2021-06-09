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

    # ------ signs and colors -----#
    oo1 = b'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'

    BORDER_COLOR = '#C7D5E0'
    
    BPAD_LEFT = ((20,10), (0, 10))
    BPAD_LEFT_INSIDE = (0, 10)

    # ------ theme -----#
    sg.theme('SystemDefault')
    #sg.theme('LightGrey4')

    # ------ set_options -----#
    sg.set_options(border_width=0, margins=(0, 0), element_padding=(0, 0))

    # ------ Menu Definition ------ #
    menu_def = [['&Information', ['&Copyright © 2021 All Right Reserved', '&By Haipeng Li at USTC', '&Report bugs to: haipengl@mail.ustc.edu.cn']] ]

    # ------ The tab layouts ------ #
    tab_sys_layout = [
                    [sg.Text('Homepath', size=(12, 1), font=("Any", 16))],
                    [sg.Input(size=(35,1), enable_events=True ,key='homepath'), sg.FolderBrowse(initial_folder='~/Desktop/')],
                    [sg.Text('CPU number',  font=("Any", 16))], 
                    [sg.Input(size=(35,1), key='mpiproc')],
                    [sg.Text('Figure aspect',  font=("Any", 16))], 
                    [sg.Input(size=(35,1), key='figaspect')],
                ]

    tab_mod_layout = [[sg.Text('Nx',  size=(14, 1), font=("Any", 16)), sg.Input(size=(16,1), key='nx')],
                    [sg.Text('Nz',  size=(14, 1), font=("Any", 16)), sg.Input(size=(16,1), key='nz')],
                    [sg.Text('Dx',  size=(14, 1), font=("Any", 16)), sg.Input(size=(16,1), key='dx')],
                    [sg.Text('Dt',  size=(14, 1), font=("Any", 16)), sg.Input(size=(16,1), key='dt')],
                    [sg.Text('Nt',  size=(14, 1), font=("Any", 16)), sg.Input(size=(16,1), key='nt')],
                    [sg.Text('PML', size=(14, 1), font=("Any", 16)), sg.Input(size=(16,1), key='pml')],
                    [sg.Text('Free surface', size=(14, 1), font=("Any", 16)), sg.Combo(('Yes', 'No'), default_value='Yes', size=(14, 1), key='fs')],
                    [sg.Text('Vp true (forward)', size=(24, 1), font=("Any", 16) )], 
                    [sg.In(size=(36,1), enable_events=True, key='vp_true_path'), sg.FileBrowse(initial_folder='~/Desktop/', file_types=(("dat Files", "*.dat"), ("Bin Files", "*.bin")))],
                    [sg.Text('              ',          size=(18, 1), font=("Any", 16))],
                    [sg.Text('Vpmax',  size=(14, 1), font=("Any", 16)), sg.Input(size=(16,1), key='vpmax')],
                    [sg.Text('Vpmin',  size=(14, 1), font=("Any", 16)), sg.Input(size=(16,1), key='vpmin')],
                    [sg.Text('SU data (inversion)', size=(24, 1), font=("Any", 16) )], 
                    [sg.Input(size=(36,1), enable_events=True ,key='field_data_path'), sg.FolderBrowse(initial_folder='~/Desktop/')],
                    [sg.Text('Vp init (inversion)', size=(24, 1), font=("Any", 16) )], 
                    [sg.In(size=(36,1), enable_events=True, key='vp_init_path'),    sg.FileBrowse(initial_folder='~/Desktop/', file_types=(("dat Files", "*.dat"), ("Bin Files", "*.bin")))],
                ]

    tab_acq_layout = [
                    [sg.Text('Marine or Land',      size=(20, 1), font=("Any", 16)), sg.Combo(('Marine','Land'), default_value='Land', key='marine_or_land', size=(14, 1))],
                    [sg.Text('Receiver coordinate file', size=(24, 1), font=("Any", 16))], 
                    [sg.In(size=(36,1), enable_events=True,key='rec_coor'), sg.FileBrowse(initial_folder='~/Desktop/', file_types=(("dat Files", "*.dat")))],
                    [sg.Text('Source coordinate file',   size=(24, 1), font=("Any", 16))], 
                    [sg.In(size=(36,1), enable_events=True, key='src_coor'), sg.FileBrowse(initial_folder='~/Desktop/', file_types=(("dat Files", "*.dat")))],
                    [sg.Text('              ',          size=(18, 1), font=("Any", 16))],
                    [sg.Text('Source wavelet',          size=(18, 1), font=("Any", 16)), sg.Combo(('Ricker (Opt 1)', 'File   (Opt 2)'), default_value='Ricker', key='wavelet', size=(14, 1))],
                    [sg.Text('Opt 1: f0 (Hz)',          size=(18, 1), font=("Any", 16)), sg.Input(size=(16,1), key='f0')],
                    [sg.Text('Opt 2: file',             size=(30, 1), font=("Any", 16))], 
                    [sg.In(size=(36,1), enable_events=True,key='wavelet_file'), sg.FileBrowse(initial_folder='~/Desktop/', file_types=(("dat Files", "*.dat")))],
                ]

    tab_inv_layout = [[sg.Text('Misfit',    size=(20, 1), font=("Any", 16)), sg.Combo(('Waveform', 'Traveltime', 'Envelope', 'Globalcorrelation'), default_value='Waveform', key='misfit_type', size=(14, 1))],
                    [sg.Text('Scheme',      size=(20, 1), font=("Any", 16)), sg.Combo(('NLCG','LBFGS'), default_value='NLCG', key='scheme', size=(14, 1))],
                    [sg.Text('Step length', size=(20, 1), font=("Any", 16)), sg.Input(size=(16,1), key='step_length')],
                    [sg.Text('Iteration',   size=(20, 1), font=("Any", 16)), sg.Input(size=(16,1), key='maxiter')],
                    [sg.Text('Gradient mute',     size=(20, 1), font=("Any", 16)), sg.Input(size=(16,1), key='grad_mute')],
                    [sg.Text('Gradient smooth',   size=(20, 1), font=("Any", 16)), sg.Input(size=(16,1), key='grad_smooth')],
                    [sg.Text('Normalization',     size=(20, 1), font=("Any", 16)), sg.Combo(('L1-Trace', 'L2-Trace', 'L1-Event', 'L2-Event', 'None'), default_value='L1-Trace', key='normalize', size=(14, 1))],
                    [sg.Text('Frequency filter',  size=(20, 1), font=("Any", 16)), sg.Combo(('None', 'Bandpass', 'Lowpass', 'Highpass'), default_value='None', key='fre_filter', size=(14, 1))],
                    [sg.Text('Frequency low',     size=(20, 1), font=("Any", 16)), sg.Input(size=(16,1), key='fre_low')],
                    [sg.Text('Frequency high',    size=(20, 1), font=("Any", 16)), sg.Input(size=(16,1), key='fre_high')],
                    [sg.Text('Mute late arrival', size=(20, 1), font=("Any", 16)), sg.Combo(('Yes', 'No'), default_value='No', size=(14, 1), key='mute_late_arrival')],
                    [sg.Text('Mute time window',  size=(20, 1), font=("Any", 16)), sg.Input(size=(16,1), key='mute_late_window')],
                    [sg.Text('Mute near offset',  size=(20, 1), font=("Any", 16)), sg.Combo(('Yes', 'No'), default_value='No', size=(14, 1), key='mute_offset_short')],
                    [sg.Text('Mute far  offset',  size=(20, 1), font=("Any", 16)), sg.Combo(('Yes', 'No'), default_value='No', size=(14, 1), key='mute_offset_long')],
                    [sg.Text('Mute near distance',size=(20, 1), font=("Any", 16)), sg.Input(size=(16,1), key='mute_offset_short_dis')],
                    [sg.Text('Mute far  distance',size=(20, 1), font=("Any", 16)), sg.Input(size=(16,1), key='mute_offset_long_dis')],
                ]


    # The TabgGroup layout - it must contain only Tabs
    tab_group_layout = [[sg.Tab('System',      tab_sys_layout, font=("Any", 16), key='sys'),
                         sg.Tab('Model',       tab_mod_layout, font=("Any", 16), key='mod'),
                         sg.Tab('Acquisition', tab_acq_layout, font=("Any", 16), key='acq'),
                         sg.Tab('Inversion',   tab_inv_layout, font=("Any", 16), key='inv'),
                        ]
                    ]

    # ------ left Columns ------ #
    block_1 = [
                [sg.Text('Input Parameters', size=(20, 1), justification='left', font=("Any", 18), relief=sg.RELIEF_RIDGE)],
                [sg.Text('_' * 60 + '\n')],
                [sg.TabGroup(tab_group_layout, enable_events=True, key='-TABGROUP-')],
                [sg.Text(' ' * 60 + '\n')],
                [sg.Button('Save parameters'), sg.Text(' ' * 5), sg.Button('Load parameters')],
            ]

    block_2 = [ [sg.Text('Functions', size=(20, 1), justification='left', font=("Any", 18), relief=sg.RELIEF_RIDGE)],
                [sg.Text('_' * 60 + '\n')],
                [sg.In(size=(10,1), enable_events=True, key='view2D_file'),             
                sg.FileBrowse(initial_folder='~/Desktop/', file_types=(("dat Files", "*.dat"), ("Bin Files", "*.bin"))),
                sg.Text('nx', size=(2, 1), justification='left', font=("Any", 16)),
                sg.In(size=(6,1), enable_events=True, key='view2D_nx'),
                sg.Text('nz', size=(2, 1), justification='left', font=("Any", 16)),
                sg.In(size=(6,1), enable_events=True, key='view2D_nz'),
                sg.Button('View',   image_data=oo1, image_subsample=2, button_color=('black', sg.theme_background_color()), border_width=0, font='Any 15')],
                [sg.In(size=(10,1), enable_events=True, key='smooth_file'),             
                sg.FileBrowse(initial_folder='~/Desktop/', file_types=(("dat Files", "*.dat"), ("Bin Files", "*.bin"))),
                sg.Text('sp', size=(2, 1), justification='left', font=("Any", 16)),
                sg.In(size=(6,1), enable_events=True, key='smooth_span'),
                sg.Text('mt', size=(2, 1), justification='left', font=("Any", 16)),
                sg.In(size=(6,1), enable_events=True, key='smooth_top_mute'),
                sg.Button('Smooth',    image_data=oo1, image_subsample=2, button_color=('black', sg.theme_background_color()), border_width=0, font='Any 15')],
                [sg.Button('Forward',  image_data=oo1, image_subsample=2, button_color=('black', sg.theme_background_color()), border_width=0, font='Any 15'),
                 sg.Button('FWI',      image_data=oo1, image_subsample=2, button_color=('black', sg.theme_background_color()), border_width=0, font='Any 15'),
                 sg.Button('Clear',    image_data=oo1, image_subsample=2, button_color=('black', sg.theme_background_color()), border_width=0, font='Any 15')],
            ] 


    # ------ middle Columns ------ #
    block_3 = [[sg.Text('Viewing Window', size=(20, 1), justification='left', font=("Any", 18), relief=sg.RELIEF_RIDGE)],
               [sg.Text('_' * 300 + '\n')],
               [sg.Image(key='figure')]
            ]

    block_4 = [ [sg.Text('Output Historys', size=(20, 1), justification='left', font=("Any", 18), relief=sg.RELIEF_RIDGE)],
                [sg.Text('_' * 300 + '\n')],
                [sg.Output(size=(140, 20), key = 'output', font=('Any 10'))],
            ]


    # ------ right Columns ------ #
    block_5 = [[sg.Text('Viewing Options', size=(20, 1), justification='left', font=("Any", 18), relief=sg.RELIEF_RIDGE)],
            [sg.Text('_' * 300 + '\n')],
            [sg.In(size=(24,1), enable_events=True ,key='fig_folder'), sg.FolderBrowse()],           
            [sg.Combo(('Acquisition', 'Source-Wavelet', 'Velocity', 'Waveform', 'Gradient', 'Direction'), default_value='Acquisition', enable_events=True, key='fig_type', size=(32, 8))],
            [sg.Text('\nFigures for viewing\n', size=(20, 1), justification='left', font=("Any", 16))],
            [sg.Listbox(values=[], enable_events=True, size=(32,35),key='fig_list')],
            ]


    block_6 = [[sg.Text('System Status', size=(20, 1), justification='left', font=("Any", 18), relief=sg.RELIEF_RIDGE)],
            [sg.Text('_' * 300 + '\n')],
            [GraphColumn('Disk Read',    '_DISK_READ_'),
             GraphColumn('Disk Write',   '_DISK_WRITE_')],
            [GraphColumn('CPU Usage',    '_CPU_'),
             GraphColumn('Memory Usage', '_MEM_')], 
            ]


    # ------ overall layout ------ #
    layout = [[sg.Menu(menu_def, tearoff=True)],
            [
                sg.Column([ [sg.Column(block_1, size=(400,600), pad=BPAD_LEFT_INSIDE)],
                            [sg.Column(block_2, size=(400,300), pad=BPAD_LEFT_INSIDE)]], 
                            pad=BPAD_LEFT, background_color=BORDER_COLOR, vertical_alignment='center', justification='center',  k='C1'),
                #sg.VSeperator(),
                sg.Column([ [sg.Column(block_3, size=(1000,600), pad=BPAD_LEFT_INSIDE)],
                            [sg.Column(block_4, size=(1000,300), pad=BPAD_LEFT_INSIDE)]], 
                            pad=BPAD_LEFT, background_color=BORDER_COLOR, vertical_alignment='center', justification='center',  k='C2'),
                #sg.VSeperator(),
                sg.Column([ [sg.Column(block_5, size=(300,600), pad=BPAD_LEFT_INSIDE)],
                            [sg.Column(block_6, size=(300,300), pad=BPAD_LEFT_INSIDE)]], 
                            pad=BPAD_LEFT, background_color=BORDER_COLOR, vertical_alignment='center', justification='center',  k='C3'),
            ]
                #    [sg.HorizontalSeparator()],
            ]



    # ------ show window ------ #

    window = sg.Window('Seismic Waveform Inversion Toolbox-1.0', layout, resizable=True, finalize=True) # , location=(0, 0), size=(1400, 800)
    # keep in the center
    window['C1'].expand(True, True, True)
    window['C2'].expand(True, True, True)
    window['C3'].expand(True, True, True)


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
            new_size = (1000,600)
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

        elif event == 'Forward':
            print('Running forward modeling ...')
            swit = prepare_forward_parameter(values)
            forward_workflow(swit)

        elif event == 'FWI':
            print('Running full-waveform inversion ...')
            swit = prepare_inversion_parameter(values)
            inversion_workflow(swit)

        elif event == 'Clear':
            window.FindElement('output').Update('')

        elif event == 'fig_folder':
            print('Opening: ' + values['fig_folder'])

        elif event == 'fig_type':
            folder = values['fig_folder']
            fig_type = values['fig_type']
            print('Viewing %s'%fig_type)

            if fig_type in ['Acquisition','Source-Wavelet']:
                suffix = ''
            elif fig_type in ['Waveform']:
                suffix = '/waveform/'
            elif fig_type in ['Velocity']:
                suffix = '/model/'
            elif fig_type in ['Gradient']:
                suffix = '/model/'
            elif fig_type in ['Direction']:
                suffix = '/model/'
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
                filename = os.path.join(values['fig_folder']+suffix, values['fig_list'][0])
                new_size = (1000,650)
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

