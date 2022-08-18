import sys
import matplotlib.pyplot as plt

sys.path.append("/home/haipeng/PowerEdge4/data/SWIT-1.0-Open/toolbox")

# import modules
import numpy as np
import base
from inversion import inversion, source_inversion
from preprocess import process_workflow
from solver import forward, source_wavelet
from tools import saveparjson, smooth2d, loadbinfloat32

vp_true = np.loadtxt('/home/haipeng/PowerEdge4/data/SWIT-1.0/thesis/Case13-Foothill/model/Foothill_801_291_25m_smooth25.dat')
grad = loadbinfloat32('/home/haipeng/PowerEdge4/data/SWIT-1.0/thesis/Case13-Foothill/outputs/gradient/grad-50.bin').reshape(801, 291)

grad_mask = np.ones_like(vp_true)
grad_mask[vp_true==1500] = 0.0

grad_mask_new = np.zeros_like(grad_mask)
window_len = 9
for i in range(801):
    grad_mask_new[i, :] = smooth1d(grad_mask[i,:], window_len = window_len,window='hanning')[:-window_len+1]

grad = grad * grad_mask * grad_mask_new


np.savetxt('/home/haipeng/Desktop/grad_new.dat', grad)
np.savetxt('/home/haipeng/Desktop/grad_mask.dat', grad_mask)
np.savetxt('/home/haipeng/Desktop/grad_mask_new.dat', grad_mask_new)