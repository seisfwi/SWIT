###############################################################################
#
# SWIT: Seismic Waveform Inversion Toolbox
#
# by Haipeng Li at USTC, haipengl@mail.ustc.edu.cn
# 
# May, 2022  
#
# Write source & receiver files
#
###############################################################################


import numpy as np

# model
nx = 501
dx = 20
x_beg = 0.
x_end = (nx-1) * dx

# source
src_beg = 0.
src_end = x_end
srcn  = 51
srcxz = np.zeros((srcn, 2))
srcxz[:,0] = np.linspace(src_beg, src_end, srcn)
srcxz[:,1] = np.linspace(   dx,    dx, srcn)

# receiver
rec_beg = x_beg
rec_end = x_end
recn  = nx
recxz = np.zeros((recn, 2))
recxz[:,0] = np.linspace(rec_beg, rec_end, recn)
recxz[:,1] = np.linspace(   dx,    dx, recn)

# save
np.savetxt('source_coordinate.dat',   srcxz)
np.savetxt('receiver_coordinate.dat', recxz)

# print
print('Source   number: %-4d, from distance %.2f ~ %.2f m'%(srcn, srcxz[0,0], srcxz[-1,0]))
print('                       from depth    %.2f ~ %.2f m'%(      srcxz[0,1], srcxz[-1,1]))
print('Receiver number: %-4d, from distance %.2f ~ %.2f m'%(recn, recxz[0,0], recxz[-1,0]))
print('                       from depth    %.2f ~ %.2f m'%(      recxz[0,1], recxz[-1,1]))
print('Successfully save source and receiver files!')