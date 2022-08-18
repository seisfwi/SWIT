###############################################################################
#
# SWIT: Seismic Waveform Inversion Toolbox
#
# by Haipeng Li at USTC, haipengl@mail.ustc.edu.cn
# 
# June, 2021  
#
# inversion module
#
###############################################################################

import copy
import numpy as np

from plot import plot_inv_scheme, plot_wavelet, plot_trace, plot_rtm
from misfit import misfit
from preprocess import process_workflow
from solver import adjoint, forward
from tools import array2vector, save_inv_scheme, smooth2d, vector2array
from tools import loadsu, su2array, get_offset, convert_wavelet_su


'''
         d_old            d_now    
vp_old  ------>  vp_now  ------>  vp_try  ------> make final decision
         g_old            g_now    
'''

def inversion(simu, optim, inv_model):
    ''' inversion workflow
    '''

    # set the initial model
    simu.model.vp = inv_model['vp'] 
    simu.model.rho = inv_model['rho'] 

    # initialize the inversion
    inv_scheme = optimize_init(simu, optim)
    inv_scheme['v_now'] = array2vector(simu.model.vp)
    inv_scheme['v_old'] = array2vector(simu.model.vp)

    while optim.iter < optim.maxiter:
        optim.iter += 1
        print('\n-----------  iteration %d  -----------\n'%optim.iter)

        # synthetic data from the current model
        forward(simu, simu_type='syn', savesnap=1)

        # plot_trace(simu, 'syn', simu_type='syn', suffix='', src_space=1, trace_space=5, scale = 0.8, color='b')
        # plot_trace(simu, 'syn-proc', simu_type='syn', suffix='_proc', src_space=1, trace_space=5, scale = 0.8, color='b')

        # process the synthetic data
        process_workflow(simu, optim, simu_type='syn')

        # evaluate the misfit
        inv_scheme['f_now'] = misfit(simu, optim.misfit_type)   

        # evaluate the gradient (with preconditioning and scaling)
        inv_scheme['g_now'] = adjoint(simu, optim) 
        
        # Non-linear Conjugate Gradient algorithm
        if optim.scheme in ['NLCG']:
            inv_scheme = NLCG(inv_scheme)
            
        # Limited-memory Broyden–Fletcher–Goldfarb–Shanno (quasi-Newton algorithm) 
        elif optim.scheme in ['LBFGS']:
            inv_scheme = LBFGS(inv_scheme)

        # search for the step length
        inv_scheme = backtrack(simu, optim, inv_scheme)

        # make the final decision (accept result or restart this iteration)
        inv_scheme = desicion_make(simu, optim, inv_scheme)
    
        # save and plot current outputs
        save_inv_scheme(simu, optim, inv_scheme)
        plot_inv_scheme(simu, optim, inv_scheme)

    print('\n-----------  iteration end  -----------\n')



def rtm(simu, optim, inv_model):
    ''' inversion workflow
    '''

    # set the initial model
    simu.model.vp = inv_model['vp'] 
    simu.model.rho = inv_model['rho'] 

    # initialize the inversion
    inv_scheme = optimize_init(simu, optim)
    inv_scheme['v_now'] = array2vector(simu.model.vp)
    inv_scheme['v_old'] = array2vector(simu.model.vp)

    # synthetic data from the current model
    forward(simu, simu_type='syn', savesnap=1)

    # compute the RTM image
    inv_scheme['g_now'] = adjoint(simu, optim) 
    
    # plot and save
    plot_rtm(simu, optim, inv_scheme)

    print('\n-----------  RTM end  -----------\n')




def NLCG(inv_scheme):
    ''' Non-linear Conjugate Gradient (NLCG) method
    '''

    # parameters of NLCG
    cg_thresh = 1.0  # thresh for checking conjugacy
    cg_itmax  = 6    # max allowed NLCG iteration and then set back to the steepest descent direction

    # get status
    g_now = inv_scheme['g_now']
    g_old = inv_scheme['g_old']
    d_old = inv_scheme['d_old']
    cg_it = inv_scheme['cg_it']
    
    # every search +1
    cg_it += 1

    # compute search direction
    if cg_it == 1:
        d_now = - g_now
        status = 0
    elif cg_it > cg_itmax:
        d_now = - g_now
        status = 1
        print('restarting NLCG... [periodic restart]')
    else:
        beta = pollak_ribere(g_now, g_old)
        if beta <= 0.:
            beta = 0.
            status = 1
            print('restarting NLCG... [negative beta]')
        else:
            status = 0
        # new direction
        d_now = - g_now + beta * d_old

    # check restart conditions
    '''
    # check conjugacy
    if check_conjugacy(g_now, g_old) > cg_thresh:
        d_now = -g_now
        cg_status = 1
        print('restarting NLCG... [loss of conjugacy]')
    '''
    # check descent direction
    if check_descent(d_now, g_now) > 0.:
        d_now = - g_now
        status = 1
        print('restarting NLCG... [not a descent direction]')
    else:
        print('find NLCG direction')

    # restart NLCG
    if status:
        cg_it = 0

    # update direction
    inv_scheme['d_now'] = d_now
    inv_scheme['cg_it'] = cg_it

    return inv_scheme


def LBFGS(inv_scheme):
    '''  Limited-memory Broyden–Fletcher–Goldfarb–Shanno (quasi-Newton algorithm) 
    '''

    # get status
    memory = inv_scheme['memory']
    memory_used = inv_scheme['memory_used']
    lbfgs_it = inv_scheme['lbfgs_it']

    g_now = inv_scheme['g_now']
    g_old = inv_scheme['g_old']
    v_now = inv_scheme['v_now']
    v_old = inv_scheme['v_old']
    
    nx = inv_scheme['nx']
    nz = inv_scheme['nz']

    m = len(g_now)
    n = memory
    path = './outputs/LBFGS_memory/'

    # every search +1
    lbfgs_it += 1

    # compute search direction
    if  lbfgs_it == 1:
        d_now = - g_now
        status = 0
    else:
        # update vectors s and y
        s = v_now - v_old
        y = g_now - g_old

        # implenment safeguard smooth in case of very small update
        # s = array2vector(smooth2d(vector2array(s, nx, nz), span=3))
        # y = array2vector(smooth2d(vector2array(y, nx, nz), span=3))

        if memory_used == 0:
            S = np.memmap(path+'S', mode='w+', dtype='float32', shape=(m, n))
            Y = np.memmap(path+'Y', mode='w+', dtype='float32', shape=(m, n))
            S[:, 0] = s
            Y[:, 0] = y
            memory_used = 1
        else:
            S = np.memmap(path+'S', mode='r+', dtype='float32', shape=(m, n))
            Y = np.memmap(path+'Y', mode='r+', dtype='float32', shape=(m, n))
            S[:, 1:] = S[:, :-1]
            Y[:, 1:] = Y[:, :-1]
            S[:, 0] = s
            Y[:, 0] = y

            if memory_used < memory:
                memory_used += 1

        k = memory_used
        rho = np.zeros(k)
        alpha = np.zeros(k)
        q = g_now

        # the first loop
        for i in range(k):
            rho[i] = 1. / np.dot(Y[:, i], S[:, i])
            alpha[i] = rho[i] * np.dot(S[:, i], q)
            q = q - alpha[i] * Y[:, i]

        # by Liu and Nocedal 1989
        invH = np.dot(Y[:, 0], S[:, 0]) / np.dot(Y[:, 0], Y[:, 0])
        z = invH * q
    
        # the second loop
        for i in range(k-1, -1, -1):
            beta = rho[i] * np.dot(Y[:, i], z)
            z = z + (alpha[i] - beta) * S[:, i]

        if check_descent(-z, g_now) > 0:
            d_now = - g_now
            status = 1
            print('restarting LBFGS... [not a descent direction]')
        else:
            d_now = -z
            status = 0
    
    # update direction and L-BFGS parameter
    if status:
        lbfgs_it = 1
        memory_used = 0
        S = np.memmap(path+'S', mode='r+')
        Y = np.memmap(path+'Y', mode='r+')
        S[:] = 0.
        Y[:] = 0.
    else:
        print('find LBFGS direction')

    inv_scheme['memory_used'] = memory_used
    inv_scheme['lbfgs_it'] = lbfgs_it
    inv_scheme['d_now'] = d_now

    return inv_scheme


def backtrack(simu, optim, inv_scheme):
    ''' backtrack linesearch
    '''
    # define parameters for line search
    dec = 0.5          # half the step
    inc = 2.1          # double the step
    ftol = 1e-4        # control the accuracy of the line search routine
    wolfe = 0.9        # coefficient for the Wolfe condition
    search_max = 6     # line search max iteration
    vmax_thresh = 300  # when direction is too large
    vmin_thresh = 20   # when direction is too small

    # get parameters
    misfit_type = optim.misfit_type
    vpmax = optim.vpmax
    vpmin = optim.vpmin
    nx = simu.model.nx
    nz = simu.model.nz

    # current gradient, direction and velocity
    # the gradient has been rescaled with respective to the model
    f_now = inv_scheme['f_now']
    v_now = inv_scheme['v_now']
    d_now = inv_scheme['d_now']
    g_now = inv_scheme['g_now']

    # initial step setup
    if optim.scheme in ['LBFGS']:
        if inv_scheme['lbfgs_it'] == 1:    # first interation, steepest descent direction
            step = optim.step_length
        else:                              # use the safeguard step length
            if max(abs(d_now)) > vmax_thresh:
                step = 1.0 * vmax_thresh / max(abs(d_now))
            elif max(abs(d_now)) < vmin_thresh:
                step = 1.0 * vmin_thresh / max(abs(d_now))
            else:
                step = 1.0
    elif optim.scheme in ['NLCG']:
        step = optim.step_length
    
    # start line search
    ls_iter = 0

    while True:
        ls_iter += 1
        f_try = 0.

        # try model
        v_try = model_update(v_now, d_now * step, vpmax, vpmin)
        simu.model.vp = vector2array(v_try, nx, nz)
        # forward simulation
        forward(simu, simu_type='syn', savesnap=0)
        # process the syntetic data
        process_workflow(simu, optim, simu_type='syn')
        # misfit from try model
        f_try = misfit(simu, misfit_type)

        # check the sufficient decrease condition (Armijo condition).
        if f_try >= f_now: # f_now + (step * dgtest)
            print("line search: step = %.4f, f_now = %.4e, f_try = %.4e. Not decrease." % (step, f_now, f_try))
            width = dec
        # Not implemented yet: check the Wolfe condition. g is the gradient of f(xk + step * d)  
        
        else:
            ls_status = 0
            inv_scheme = updata_inv_scheme(inv_scheme, ls_status, ls_iter, step, f_try, v_try)
            print("line search: step = %.4f, f_now = %.4e, f_try = %.4e. Succeed." % (step, f_now, f_try))
            return inv_scheme

        if ls_iter >= search_max:
            ls_status = -1
            inv_scheme = updata_inv_scheme(inv_scheme, ls_status, ls_iter, step, f_try, v_try)
            print("line search fails: exceed the iteration of linesearch.")
            return inv_scheme

        # update the step length
        step = step * width


def optimize_init(simu, optim):
    ''' initilize the optimize parameter
    '''
    nx = simu.model.nx
    nz = simu.model.nz
    inv_scheme = {}

    if optim.scheme in ['NLCG']:
        temp = np.zeros(nx*nz, dtype=np.float32)
        inv_scheme = {'g_now': temp,
                      'g_old': temp,
                      'd_now': temp,
                      'd_old': temp,
                      'v_old': temp,
                      'v_now': temp,
                      'v_try': temp,
                      'f_old':    0.,
                      'f_now':    0.,
                      'f_try':    0.,
                      'cg_it':    0 ,
                      'ls_iter':  0 ,
                      'ls_status':0 ,
                      'ls_step':  0.,
                      'nx': nx,
                      'nz': nz,
                      }

    elif optim.scheme in ['LBFGS']:
        # use the default memory of 5 in L-BFGS
        memory = 5
        temp = np.zeros(nx*nz, dtype=np.float32)
        memo = np.zeros((nx*nz, memory), dtype=np.float32)

        inv_scheme = {'memory': memory,
                      'memory_used': 0,
                      'lbfgs_it':    0,
                      's':     memo,
                      'y':     memo,
                      'g_now': temp,
                      'g_old': temp,
                      'd_now': temp,
                      'v_old': temp,
                      'v_now': temp,
                      'v_try': temp,
                      'f_old':     0.,
                      'f_now':     0.,
                      'f_try':     0.,
                      'ls_iter':   0 ,
                      'ls_status': 0 ,
                      'ls_step':   0.,
                      'nx':        nx,
                      'nz':        nz,
                      }

    return inv_scheme


def desicion_make(simu, optim, inv_scheme):
    ''' make final update desicion or restart the inversion from the last iteration
    '''
    nx = simu.model.nx
    nz = simu.model.nz

    if optim.scheme in ['NLCG']:
        # when the line search fails and CG iteration > 1
        if inv_scheme['ls_status'] == -1 and inv_scheme['cg_it'] > 1:
            inv_scheme['cg_it'] = 0
            inv_scheme['g_old'] = inv_scheme['g_old'] * 0. 
            inv_scheme['d_old'] = inv_scheme['d_old'] * 0. 
            print('line search fails: restarting the line search ...')

        # when the line search succeeds or CG iteration = 1
        else:
            inv_scheme['g_old'] = inv_scheme['g_now']
            inv_scheme['d_old'] = inv_scheme['d_now']
            inv_scheme['f_old'] = inv_scheme['f_now']
            inv_scheme['f_now'] = inv_scheme['f_try']
            inv_scheme['v_now'] = inv_scheme['v_try']
            print('accept line search result')

    elif optim.scheme in ['LBFGS']:
        # when the line search fails
        if inv_scheme['ls_status'] == -1:
            # initilize the inv_scheme again
            inv_scheme_old = copy.deepcopy(inv_scheme)
            inv_scheme = optimize_init(simu, optim)

            inv_scheme['v_old'] = inv_scheme_old['v_now']
            inv_scheme['v_now'] = inv_scheme_old['v_now']
            inv_scheme['f_now'] = inv_scheme_old['f_now']
            inv_scheme['g_now'] = inv_scheme_old['g_now']
            inv_scheme['d_now'] = inv_scheme_old['d_now']
            print('line search fails: restarting the line search ...')

        # when the line search succeeds or CG iteration = 1
        else:
            inv_scheme['g_old'] = inv_scheme['g_now']
            inv_scheme['v_old'] = inv_scheme['v_now']
            inv_scheme['v_now'] = inv_scheme['v_try']
            inv_scheme['f_old'] = inv_scheme['f_now']
            inv_scheme['f_now'] = inv_scheme['f_try']
            print('accept line search result')

    # update the velocity model
    simu.model.vp = vector2array(inv_scheme['v_now'], nx, nz)

    return inv_scheme


def updata_inv_scheme(inv_scheme, ls_status, ls_iter, ls_step, f_try, v_try):
    ''' update line search results 
    '''
    inv_scheme['ls_status'] = ls_status
    inv_scheme['ls_iter'] = ls_iter
    inv_scheme['ls_step'] = ls_step
    inv_scheme['f_try'] = f_try
    inv_scheme['v_try'] = v_try

    return inv_scheme


def model_update(vp, update, vpmax, vpmin):
    ''' update model and set max/min bound
    '''
    vp_new = np.zeros_like(vp)
    vp_new = vp + update
    vp_new[vp_new > vpmax] = vpmax
    vp_new[vp_new < vpmin] = vpmin

    return vp_new


def pollak_ribere(g_now, g_old):
    ''' pollak_ribere CG scheme
    '''
    num = np.dot(g_now, g_now-g_old)
    den = np.dot(g_old, g_old)
    beta = num/den

    return beta


def check_conjugacy(g_now, g_old):
    ''' check conjugacy
    '''
    return abs(np.dot(g_now, g_old) / np.dot(g_now, g_now))


def check_descent(d_now, g_now):
    ''' check descent direction
    '''
    return np.dot(d_now, g_now) / np.dot(g_now, g_now)


def source_inversion(simu, inv_offset=10000):
    ''' source signature inversion in the manner of Parrat (1999)
    '''
    # get prepared
    homepath = simu.system.homepath
    stf_now = simu.source.wavelet
    srcx = simu.source.xz[:,0]
    srcn = simu.source.n
    dt = simu.model.dt

    stf_inv = np.zeros_like(stf_now)

    for isrc in range(srcn):

        # load data and set offset
        obs = loadsu(homepath + 'data/obs/src%d_sg_proc.su'%(isrc+1))
        syn = loadsu(homepath + 'data/syn/src%d_sg_proc.su'%(isrc+1))
        if not np.array_equal(get_offset(obs), get_offset(syn)):
            raise ValueError("offset not consistant.")
        offset = get_offset(syn)
        obs = su2array(obs)
        syn = su2array(syn)

        # select data for source inversion
        rec_used = np.argwhere(abs(offset) < inv_offset)
        n1 = int(rec_used[0])
        n2 = int(rec_used[-1])
        data_obs = np.squeeze(obs[n1:n2, :]).T
        data_syn = np.squeeze(syn[n1:n2, :]).T
        Do = np.fft.fft(data_obs, axis=0)  # frequency domain
        Dm = np.fft.fft(data_syn, axis=0)  # frequency domain
        
        # current source wavelet
        src = np.squeeze(stf_now[isrc, :])
        S = np.fft.fft(np.squeeze(src), axis=0) # frequency domain

        # check
        if abs(np.sum(Dm * np.conj(Dm))) == 0:
            raise ValueError("No trace for source inversion, check for the reason.")
        
        # source inversion
        A = np.sum(np.conj(Dm)*Do, axis=1) / np.sum(Dm * np.conj(Dm), axis=1)
        temp = np.real(np.fft.ifft(A*S[:]))
        temp = temp / np.max(abs(temp))
        stf_inv[isrc, :] = temp

    # import matplotlib.pyplot as plt
    # plt.imshow(stf_inv[:,0:500].T, aspect=0.1)
    # plt.savefig('stf_inv.png')

    # save source wavelet plots
    plot_wavelet(simu, stf_now, 'stf_now', scale=1.0, color='k', plot_dx=5000, t_end = 1.0)
    plot_wavelet(simu, stf_inv, 'stf_inv', scale=1.0, color='r', plot_dx=5000, t_end = 1.0)
   
    # save source wavelet data
    np.savetxt(homepath+'outputs/stf_now.dat', stf_now)
    np.savetxt(homepath+'outputs/stf_inv.dat', stf_inv)

    print('Source inversion finished\n')


