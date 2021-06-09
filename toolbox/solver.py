###############################################################################
#
# SWIT: Seismic Waveform Inversion Toolbox
#
# by Haipeng Li at USTC, haipengl@mail.ustc.edu.cn
# 
# June, 2021  
#
# Forward and adjoint Solver module
#
###############################################################################


import os
import subprocess
import numpy as np
from scipy import integrate

from misfit import adjoint_source
from postprocess import grad_precond
from tools import cleandata, loadbinfloat32, savebinfloat32


def forward(simu, simu_type='obs', savesnap=0):
    ''' forward solver: fd2dmpi (2D acoustic)
    '''
    if simu_type not in ['obs', 'syn']:
        raise ValueError('Forward: unsupport simulation tyep.')

    # get prepared
    mpiproc = simu.system.mpiproc
    homepath = simu.system.homepath
    srcn = simu.source.n
    parfile = homepath + 'parfile/forward_parfile/parfile'

    # change path and clean previous data (always, even empty)
    if simu_type in ['syn']: 
        cleandata(homepath + 'data/syn/')

    # create working directory
    for isrc in range(srcn):
        ifolder = homepath + 'data/%s/src%d_snapshot'%(simu_type, isrc+1)
        if not os.path.exists(ifolder):
            os.system('mkdir %s' % ifolder)    

    # prepare the forward source
    src = integrate.cumtrapz(simu.source.wavelet, axis=-1, initial=0)

    # write parameters and model files
    write_parfile(simu, simu_type, src, savesnap=savesnap)

    # submit job
    os.chdir(homepath)
    solver_cmd = 'mpirun -np %d  fd2dmpi par=%s' % (mpiproc, parfile)
    status = subprocess.getstatusoutput(solver_cmd)
    if status[0]:
        print(status[1])
        raise ValueError('Forward solver crash')


def adjoint(simu, optim):
    ''' adjoint solver: fd2dmpi (2D acoustic)
    '''

    # get prepared
    nx = simu.model.nx
    nz = simu.model.nz
    srcn = simu.source.n
    mpiproc = simu.system.mpiproc
    homepath = simu.system.homepath
    parfile = homepath + 'parfile/adjoint_parfile/parfile'
    misfit_type = optim.misfit_type

    # always clean previous data
    cleandata(homepath + 'data/adj/')

    # create working directory
    for isrc in range(srcn):
        ifolder = homepath + 'data/adj/src%d_snapshot'%(isrc+1)
        if not os.path.exists(ifolder):
            os.system('mkdir %s' % ifolder)

    # prapare the adjoint source
    src = adjoint_source(simu, misfit_type)

    # write parameters and model files
    write_parfile(simu, 'adj', src, savesnap=0)

    # summit job
    os.chdir(homepath)
    solver_cmd = 'mpirun -np %d fd2dmpi par=%s' % (mpiproc, parfile)
    status = subprocess.getstatusoutput(solver_cmd)
    if status[0]:
        raise ValueError('Adjoint solver crash')

    # load and merge gradients and illuminations
    grad = np.zeros(nx*nz, dtype=np.float32)
    forw = np.zeros(nx*nz, dtype=np.float32)
    back = np.zeros(nx*nz, dtype=np.float32)
    for isrc in range(srcn):
        grad = grad + loadbinfloat32(homepath+'data/adj/src%d_kernel_vp.bin'  % (isrc+1))
        forw = forw + loadbinfloat32(homepath+'data/adj/src%d_illum_forw.bin' % (isrc+1))
        back = back + loadbinfloat32(homepath+'data/adj/src%d_illum_back.bin' % (isrc+1))
    
    # gradient precondtioning
    return grad_precond(simu, optim, grad, forw, back)



def write_parfile(simu, simu_type, src, savesnap = 0):
    ''' write parameter file for the fd2dmpi solver
    '''

    # get prepared
    savestep = simu.model.savestep

    # parameter file path
    path = simu.system.homepath
    if simu_type in ['obs', 'syn']:
        parpath = path + 'parfile/forward_parfile/parfile'
        srcpath = path + 'parfile/forward_source/'
    elif simu_type in ['adj']:
        parpath = path + 'parfile/adjoint_parfile/parfile'
        srcpath = path + 'parfile/adjoint_source/'

    # save source time function
    for isrc in range(simu.source.n):
        savebinfloat32(srcpath + 'src%d.bin' % (isrc+1), src[isrc, :])

    # save P-wave velocity and density files
    vp  = simu.model.vp
    rho = simu.model.rho
    savebinfloat32(path + 'parfile/model/vel.bin', vp)
    savebinfloat32(path + 'parfile/model/rho.bin', rho)

    # Write Parameter file for acoustic solver
    fp = open(parpath, "w")
    fp.write('######################################### \n')
    fp.write('#                                         \n')
    fp.write('#     fd2dmpi input parameter file        \n')
    fp.write('#                                         \n')
    fp.write('######################################### \n')
    fp.write('                                          \n')

    if simu_type in ['obs', 'syn']:
        fp.write('jobtype=forward_modeling\n')
    elif simu_type in ['adj']:
        fp.write('jobtype=adjoint_modeling\n')
    fp.write('COORD_FILE=%s\n'       % (path + 'parfile/model/coord.txt'))
    fp.write('DATA_OUT=%s\n'         % (path + 'data/'+ simu_type + '/src'))
    fp.write('VEL_IN=%s\n'           % (path + 'parfile/model/vel.bin'))
    fp.write('DENSITYFILE=%s\n'      % (path + 'parfile/model/rho.bin'))
    fp.write('FILEFORMAT=su\n'                      )
    fp.write('NX=%d\n'               % simu.model.nx)
    fp.write('NZ=%d\n'               % simu.model.nz)
    fp.write('DX=%f\n'               % simu.model.dx)
    fp.write('NPML=%d\n'             % simu.model.pml)
    fp.write('NT_WORK=%d\n'          % simu.model.nt)
    fp.write('DT_WORK=%f\n'          % simu.model.dt)
    if simu.model.fs:      # Free surface
        fp.write('FREESURFACE=1\n')
    else:
        fp.write('FREESURFACE=0\n')
    fp.write('STORE_SNAP=%d\n'       % savesnap)
    fp.write('STORE_STEP=%d\n'       % savestep)
    fp.close()

    ####################################################################
    #
    # Write geometry file for acoustic solver. The format is as follows:
    #
    # column           1        2      3    4     5     6        7
    # meaning     ( S_index  R_index  Sx    Sz    Rx    Rz  S_is_alive(0/1) )
    #
    ####################################################################

    srcn = simu.source.n
    recn = simu.receiver.n
    srcxz = simu.source.xz
    recxz = simu.receiver.xz

    fp = open(path + 'parfile/model/coord.txt', "w")
    geom = np.zeros((recn*srcn, 7))
    for isrc in range(srcn):
        for irec in range(recn):
            ig = irec + recn*(isrc-1)
            geom[ig, :] = np.array([isrc+1, irec+1, srcxz[isrc, 0], srcxz[isrc, 1],
                                    recxz[isrc, irec, 0], recxz[isrc, irec, 1], 1])
            fp.write('%6i %6i %10.1f %10.1f %10.1f %10.1f %10.1f\n' % (
                geom[ig, 0], geom[ig, 1], geom[ig, 2], geom[ig, 3], geom[ig, 4], geom[ig, 5], geom[ig, 6]))
    fp.close()

 
def source_wavelet(nt, dt, f0, srctype):
    ''' source time function
    '''
    # time and wavelet arrays
    t       = np.linspace(0, dt*nt, num=nt, endpoint=False)
    wavelet = np.zeros_like(t)

    if srctype.lower() in ['ricker']:
        t0 = 1.2/f0
        temp = (np.pi*f0) ** 2
        wavelet = (1 - 2 * temp * (t - t0) ** 2) * np.exp(- temp * (t - t0) ** 2)
    else:
        raise ValueError('Other wavelets can be implemented here.')

    return wavelet
