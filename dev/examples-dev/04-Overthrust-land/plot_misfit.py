###############################################################################
# SWIT v1.1: Seismic Waveform Inversion Toolbox
#
# Description: This script is used to plot the data misfit and model misfit
# 
###############################################################################

import os
import numpy as np
import matplotlib.pyplot as plt

## load data misfit
misfit_file = 'iterate_CG.log'
misfit_data = []
with open(misfit_file) as f:
    for line in f:
        values = line.split()
        try:            
            misfit_data.append(float(values[1]))
        except:
            pass
misfit_data = np.array(misfit_data)


## load model misfit
nx,  nz     = [501,  151]                         # Grid number along x and z directions
maxiter     = 20                                  # maximum iteration number,  i.e., 20

# true model
workpath      = '/scr2/haipeng/SWIT-1.1/04-Overthrust-land' 
misfit_model  = np.zeros((maxiter, 1))
vp_true       = np.load('./acquisition/vp_true.npy')

# load FWI models
for it in range(maxiter):
    vp_FWI = np.load(os.path.join(workpath, "fwi/model/vp_it_{:04d}.npy".format(it+1))).reshape((nx, nz))
    misfit_model[it] = np.sum(np.power((vp_true - vp_FWI), 2))

# normlize
misfit_data  /= misfit_data[0]
misfit_model /= misfit_model[0] 

# plot
fontsize, ticksize = [16, 14]
fig, axes = plt.subplots(nrows=1, ncols=2, figsize=[11, 6])
ax = fig.add_subplot(axes[0])
ax.plot(misfit_data, 'ro-')
plt.xlim([0, maxiter+1])
plt.ylim([0, 1.1])
ax.set_xticks(np.arange(0, maxiter+1, 2))
ax.set_yticks(np.arange(0, 1.1, 0.2))
ax.xaxis.set_label_text('Iterations', fontsize=fontsize)
ax.yaxis.set_label_text('Data misfit', fontsize=fontsize)

ax = fig.add_subplot(axes[1])
ax.plot(misfit_model, 'g*-')
plt.xlim([0, maxiter+1])
plt.ylim([0, 1.1])
ax.set_xticks(np.arange(0, maxiter+1, 2))
ax.set_yticks(np.arange(0, 1.1, 0.2))
ax.xaxis.set_label_text('Iterations', fontsize=fontsize)
ax.yaxis.set_label_text('model misfit', fontsize=fontsize)

# save figure to results folder
if not os.path.exists('./results'):
    os.makedirs('./results')

plt.savefig('./results/misfit-CG.png', dpi=300)
