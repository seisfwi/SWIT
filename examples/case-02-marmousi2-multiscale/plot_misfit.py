import numpy as np
import matplotlib.pyplot as plt
from tools import loadbinfloat32

# size
nx,  nz     = [481,  121]                           # Grid number along x and z directions
maxiter     = 20                                   # maximum iteration number,  i.e., 20
nband       = 5 

# initial
misfit_model  = np.zeros((nband, maxiter))
misfit_data   = np.zeros((nband, maxiter+1))
vp_true       = np.loadtxt('./model/Marmousi_481_121_25m.dat').T

# load FWI models
for iband in range(nband):
    for it in range(maxiter):
        temp = np.loadtxt("./freq-band-%d/outputs/misfit_data.dat"%(iband+1))
        misfit_data[iband, :] =  temp/ temp[0]
        vp_FWI = loadbinfloat32('./freq-band-%d/outputs/velocity/vp-%d.bin'%(iband+1, it+1)).reshape((nx, nz)).T
        misfit_model[iband, it] = np.sum( np.power((vp_true - vp_FWI), 2))

# normlize
misfit_data  = misfit_data.reshape((maxiter+1)*nband,1)
misfit_model = misfit_model.reshape(maxiter*nband,1)

misfit_data  /= misfit_data[0]
misfit_model /= misfit_model[0] 

# plot
fontsize, ticksize = [16, 14]
fig, axes = plt.subplots(nrows=1, ncols=2, figsize=[16, 5])
ax = fig.add_subplot(axes[0])
ax.plot(misfit_data, 'ro-')
plt.xlim([0, (maxiter+1)*nband+1])
plt.ylim([0, 1.0])
ax.xaxis.set_label_text('Iterations', fontsize=fontsize)
ax.yaxis.set_label_text('Data misfit', fontsize=fontsize)

ax = fig.add_subplot(axes[1])
ax.plot(misfit_model, 'g*-')
plt.xlim([0,   (maxiter+1)*nband+1])
plt.ylim([0., 1.0])
ax.xaxis.set_label_text('Iterations', fontsize=fontsize)
ax.yaxis.set_label_text('model misfit', fontsize=fontsize)

plt.savefig('misfit.png', dpi=300)





