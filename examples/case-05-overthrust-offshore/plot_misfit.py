import numpy as np
import matplotlib.pyplot as plt
from tools import loadbinfloat32

# size
nx,  nz     = [501,  171]                         # Grid number along x and z directions
maxiter     = 50                                   # maximum iteration number,  i.e., 20

# initial
misfit_model  = np.zeros((maxiter, 1))
misfit_data   = np.loadtxt("./outputs/misfit_data.dat")
vp_true       = np.loadtxt('./model/Overthrust_501_171_20m.dat').T

# load FWI models
for it in range(maxiter):
    vp_FWI = loadbinfloat32('./outputs/velocity/vp-%d.bin'%(it+1)).reshape((nx, nz)).T
    misfit_model[it] = np.sum( np.power((vp_true - vp_FWI), 2))

# normlize
misfit_data  /= misfit_data[0]
misfit_model /= misfit_model[0] 

# plot
fontsize, ticksize = [16, 14]
fig, axes = plt.subplots(nrows=1, ncols=2, figsize=[11, 6])
ax = fig.add_subplot(axes[0])
ax.plot(misfit_data, 'ro-')
plt.xlim([0, maxiter+1])
plt.ylim([0, 1.0])
ax.xaxis.set_label_text('Iterations', fontsize=fontsize)
ax.yaxis.set_label_text('Data misfit', fontsize=fontsize)

ax = fig.add_subplot(axes[1])
ax.plot(misfit_model, 'g*-')
plt.xlim([0,   maxiter+1])
plt.ylim([0., 1.0])
ax.xaxis.set_label_text('Iterations', fontsize=fontsize)
ax.yaxis.set_label_text('model misfit', fontsize=fontsize)

plt.savefig('misfit.png', dpi=300)





