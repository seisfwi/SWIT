###############################################################################
#
# SWIT v1.1: Seismic Waveform Inversion Toolbox
#
#   A Python package for seismic waveform inversion
#   By Haipeng Li at USTC & Stanford
#   Email: haipengl@mail.ustc.edu.cn, haipeng@stanford.edu 
#
#   Plot module defines the plotting functions
#
###############################################################################

import os

import matplotlib
import matplotlib.pyplot as plt
import numpy as np


def plot_waveform_comparison(path, method, iter, niter_max, src_num):
    ''' Plot waveform comparison
    '''
    pass



def plot_model(x, z, data, vmin, vmax, filename, title, figaspect = 1, colormap = 'bwr'):
    ''' Plot model (2D), e.g. vp, grad.

    Parameters
    ----------
        x : 1D array (float32)
            x axis
        z : 1D array (float32)
            z axis
        data : 2D array (float32)   
            model data
        vmin : float
            minimum value to plot
        vmax : float
            maximum value to plot
        filename : str  
            path to save the figure
        title : str
            title of the figure
        figaspect : float
            aspect ratio of the figure
        colormap : str
            colormap
    '''
    
    plotopts = {
        'vmin': vmin,
        'vmax': vmax,
    }

    plotopts['extent'] = (x[0]/1000, x[-1]/1000, z[-1]/1000, z[0]/1000)
    xlabel = 'Distance (km)'
    ylabel = 'Depth (km)'

    # colormap
    if colormap in ['my_seismic_cmap']:
        plotopts['cmap'] = my_seismic_cmap()
    else:
        plotopts['cmap'] = colormap

    fig = plt.figure(figsize=(10, 5))
    ax = fig.add_subplot(111)
    im = ax.imshow(data, **plotopts)
    fig.colorbar(im, shrink=0.5, extend='both')
    ax.xaxis.set_label_text(xlabel)
    ax.yaxis.set_label_text(ylabel)
    ax.set_title(title)
    ax.set_aspect(figaspect)
    plt.savefig(filename, dpi=300)
    plt.close()
    # plt.show()


def plot_misfit(path, method, iter, niter_max, src_num):
    ''' plot misfit

    Parameters
    ----------
        path : str
            working path
        method : str
            method name in optimization
        iter : int  
            number of iterations
        niter_max : int
            maximum number of iterations
        src_num : int
            number of sources
    '''
    
    fig_path = os.path.join(path, 'fwi/figures/')

    # plot misfit in the stacked bar chart for each iteration
    misfit = np.zeros((niter_max, src_num))
    for it in range(iter):
        misfit_file = os.path.join(path, 'fwi/misfit/fcost_all_it_%04d.npy'%(it+1))
        misfit[it] = np.load(misfit_file)
        
    x = [i+1 for i in range(niter_max)]
    misfit_sum = np.zeros(niter_max)

    plt.figure(figsize=(12, 8))
    for isrc in range(src_num):
        plt.bar(x, misfit[:, isrc], bottom = misfit_sum)
        misfit_sum += misfit[:, isrc]
        
    plt.xlabel("Iterations")
    plt.ylabel("Misfit")
    plt.xticks(x)
    plt.legend([str(isrc+1) for isrc in range(src_num)], loc = 'center left', bbox_to_anchor=(1, 0.5))
    plt.title("Misfit history of all sources")
    plt.savefig(fig_path + 'misfit_history_all_source.png', dpi=300)
    # plt.show()

    # plot misfit in curve
    misfit_file = 'iterate_{}.log'.format(method)
    misfit = []
    with open(misfit_file) as f:
        for line in f:
            values = line.split()
            try:            
                misfit.append(float(values[1]))
            except:
                pass
    misfit = np.array(misfit)

    plt.figure(figsize=(10, 8))
    plt.plot(misfit / misfit[0], marker='o')
    plt.xlim(0, niter_max + 1)
    plt.ylim(0, 1.1)
    plt.xticks(np.arange(0, niter_max, 1))
    plt.yticks(np.arange(0, 1.1, 0.2))
    plt.xlabel('iteration', fontsize=14)
    plt.ylabel('Normalized misfit', fontsize=14)
    plt.savefig(fig_path + 'misfit_history.png', dpi=300)
    plt.close()


def my_seismic_cmap():
    ''' my seismic cmap
    '''
    cdict = {'red': ((0.0, 0.0, 0.0),
                    (0.1, 0.5, 0.5),
                    (0.2, 0.0, 0.0),
                    (0.4, 0.2, 0.2),
                    (0.6, 0.0, 0.0),
                    (0.8, 1.0, 1.0),
                    (1.0, 1.0, 1.0)),
            'green':((0.0, 0.0, 0.0),
                    (0.1, 0.0, 0.0),
                    (0.2, 0.0, 0.0),
                    (0.4, 1.0, 1.0),
                    (0.6, 1.0, 1.0),
                    (0.8, 1.0, 1.0),
                    (1.0, 0.0, 0.0)),
            'blue': ((0.0, 0.0, 0.0),
                    (0.1, 0.5, 0.5),
                    (0.2, 1.0, 1.0),
                    (0.4, 1.0, 1.0),
                    (0.6, 0.0, 0.0),
                    (0.8, 0.0, 0.0),
                    (1.0, 0.0, 0.0))}

    return matplotlib.colors.LinearSegmentedColormap('my_colormap',cdict,256)
    