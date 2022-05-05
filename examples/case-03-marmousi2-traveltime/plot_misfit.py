import numpy as np
import matplotlib.pyplot as plt
from tools import loadbinfloat32

# size
nx,  nz     = [481,  121]                           # Grid number along x and z directions
maxiter     = 130                                   # maximum iteration number,  i.e., 20
nstage      = 2

# initial
misfit_model  = np.zeros((maxiter, 1))
vp_true       = np.loadtxt('./model/Marmousi_481_121_25m.dat').T

# load FWI models
for it in range(maxiter):
    if it < 30:
        vp_FWI = loadbinfloat32('./traveltime/outputs/velocity/vp-%d.bin'%(it+1)).reshape((nx, nz)).T
    else:
        vp_FWI = loadbinfloat32('./globalcorrelation/outputs/velocity/vp-%d.bin'%(it-30+1)).reshape((nx, nz)).T

    misfit_model[it] = np.sum( np.power((vp_true - vp_FWI), 2))

# normlize
misfit_model /= misfit_model[0] 

# data misfit 
misfit_stage1 = np.loadtxt("./traveltime/outputs/misfit_data.dat")
misfit_stage2 = np.loadtxt("./globalcorrelation/outputs/misfit_data.dat")

misfit_stage1  /= misfit_stage1[0]
misfit_stage2  /= misfit_stage2[0]


# plot
fontsize, ticksize = [16, 14]
fig, axes = plt.subplots(nrows=1, ncols=3, figsize=[16, 5])
ax = fig.add_subplot(axes[0])
ax.plot(misfit_stage1, 'ro-')
plt.xlim([0, 31])
plt.ylim([0, 1.0])
ax.xaxis.set_label_text('Iterations', fontsize=fontsize)
ax.yaxis.set_label_text('Data misfit (stage 1)', fontsize=fontsize)

ax = fig.add_subplot(axes[1])
ax.plot(misfit_stage2, 'ro-')
plt.xlim([0, 101])
# plt.ylim([0, 1.0])
ax.xaxis.set_label_text('Iterations', fontsize=fontsize)
ax.yaxis.set_label_text('Data misfit (stage 2)', fontsize=fontsize)


ax = fig.add_subplot(axes[2])
ax.plot(misfit_model, 'g*-')
plt.xlim([0,   maxiter+1])
plt.ylim([0., 1.0])
ax.xaxis.set_label_text('Iterations', fontsize=fontsize)
ax.yaxis.set_label_text('model misfit', fontsize=fontsize)

plt.savefig('misfit.png', dpi=300)





