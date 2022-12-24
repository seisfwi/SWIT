import numpy as np


def normL2(x):
    '''The routine normL2 returns the Euclidian norm of a vector x of size n
    
    Args:

    Returns:
    '''

    return np.linalg.norm(x)


def scalL2(x, y):
    '''
     The routine scalL2 returns the scalar product between two vectors x and y
    '''

    return np.inner(x, y)



def print_info(optimize, fcost):
    ''' print  information
    '''

    ng = normL2(optimize.grad)

    if(optimize.FLAG == 'INIT'):
        f = open('iterate_%s.dat'%(optimize.method), 'w')
        if(optimize.method == 'SD'):
            f.write('**********************************************************************\n')
            f.write('         STEEEPEST DESCENT ALGORITHM         \n')
            f.write('**********************************************************************\n')
        elif(optimize.method == 'CG'):
            f.write('**********************************************************************\n')
            f.write('         NONLINEAR CONJUGATE GRADIENT ALGORITHM         \n')
            f.write('**********************************************************************\n')
        elif(optimize.method == 'LBFGS'):
            f.write('**********************************************************************\n')
            f.write('             l-BFGS ALGORITHM                \n')
            f.write('**********************************************************************\n')
        elif(optimize.method == 'PLBFGS'):
            f.write('**********************************************************************\n')
            f.write('             PRECONDITIONED l-BFGS ALGORITHM                \n')
            f.write('**********************************************************************\n')

        f.write('     Convergence criterion  : %10.2e\n' %optimize.conv)
        f.write('     Niter_max              : %10d\n'   %optimize.niter_max)
        f.write('     Initial cost is        : %10.2e\n' %optimize.f0)
        f.write('     Initial norm_grad is   : %10.2e\n' %ng)
        f.write('**********************************************************************\n')
        f.write('   Niter      fk         ||gk||       fk/f0        alpha        nls      ngrad    \n')
        f.write('%6d %12.2e %12.2e %12.2e %12.2e %8d %8d\n'
        %(optimize.cpt_iter, fcost, ng, fcost/optimize.f0, optimize.alpha, optimize.cpt_ls, optimize.nfwd_pb))
        f.close()

    elif(optimize.FLAG == 'CONV'):
        f = open('iterate_%s.dat'%(optimize.method), 'a')
        f.write('%6d %12.2e %12.2e %12.2e %12.2e %8d %8d\n'
        %(optimize.cpt_iter, fcost, ng, fcost/optimize.f0, optimize.alpha, optimize.cpt_ls, optimize.nfwd_pb))
        f.write('**********************************************************************\n')

        if(optimize.cpt_iter >= optimize.niter_max):
            f.write('  STOP: MAXIMUM NUMBER OF ITERATION REACHED    \n')
        else:
            f.write('  STOP: CONVERGENCE CRITERION SATISFIED        \n')
            f.write('**********************************************************************\n')
            f.close()

    elif(optimize.FLAG == 'FAIL'):
        f = open('iterate_%s.dat'%(optimize.method), 'a')
        f.write('**********************************************************************\n')
        f.write('  STOP: LINESEARCH FAILURE    \n')
        f.write('**********************************************************************\n')
        f.close()
    else:
        f = open('iterate_%s.dat'%(optimize.method), 'a')
        f.write('%6d %12.2e %12.2e %12.2e %12.2e %8d %8d\n'
        %(optimize.cpt_iter, fcost, ng, fcost/optimize.f0, optimize.alpha, optimize.cpt_ls, optimize.nfwd_pb))
        f.close()



def print_info_TRN(optimize, fcost):
    
    if(optimize.FLAG=='INIT') :
        f = open('iterate_TRN.dat', 'w')
        f.write('******************************************************************************************\n')
        f.write('                                 TRUNCATED NEWTON ALGORITHM                               \n')
        f.write('******************************************************************************************\n')
        f.write('     Convergence criterion  : %10.2e \n' %optimize.conv)
        f.write('     Niter_max              : %7d \n'    %optimize.niter_max)
        f.write('     Initial cost is        : %10.2e \n' %optimize.f0)
        f.write('     Initial norm_grad is   : %10.2e \n' %optimize.norm_grad)
        f.write('     Maximum CG iter        : %7d \n'    %optimize.niter_max_CG)
        f.write('******************************************************************************************\n')
        f.write('   Niter      fk         ||gk||       fk/f0          alpha      nls     nit_CG      eta       ngrad     nhess \n')
        f.write('%6d %12.2e %12.2e %12.2e %12.2e %8d %8d %12.2e %8d %8d \n'%(
        optimize.cpt_iter, fcost, optimize.norm_grad, fcost/optimize.f0, optimize.alpha, optimize.cpt_ls, optimize.cpt_iter_CG, optimize.eta, optimize.nfwd_pb, optimize.nhess))
        f.close()

        f = open('iterate_TRN_CG.dat', 'w')
        f.write('******************************************************************************************\n')
        f.write('                                 TRUNCATED NEWTON ALGORITHM                               \n')
        f.write('                                      INNER CG HISTORY                                    \n')
        f.write('******************************************************************************************\n')
        f.write('     Convergence criterion  : %10.2e \n' %optimize.conv)
        f.write('     Niter_max              : %7d \n'    %optimize.niter_max)
        f.write('     Initial cost is        : %10.2e \n' %optimize.f0)
        f.write('     Initial norm_grad is   : %10.2e \n' %optimize.norm_grad)
        f.write('     Maximum CG iter        : %7d \n'    %optimize.niter_max_CG)
        f.write('******************************************************************************************\n')
        f.close()
    elif(optimize.FLAG=='CONV'):
        f = open('iterate_TRN.dat', 'a')
        f.write('**********************************************************************\n')
        if(optimize.cpt_iter == optimize.niter_max):
            f.write('  STOP: MAXIMUM NUMBER OF ITERATION REACHED    \n')
        else:
            f.write('  STOP: CONVERGENCE CRITERION SATISFIED        \n')
        
        f.write('**********************************************************************\n')
        f.close()

    elif(optimize.FLAG=='FAIL'):
        f = open('iterate_TRN.dat', 'a')
        f.write('**********************************************************************\n')
        f.write('  STOP: LINESEARCH FAILURE    \n')
        f.write('**********************************************************************\n')
        f.close()
    elif(optimize.comm=='DESC'):
        f = open('iterate_TRN_CG.dat', 'a')
        if(optimize.CG_phase=='INIT'):
            f.write('-------------------------------------------------------------------------------------------------\n')
            f.write(' NONLINEAR ITERATION %4d  ETA IS : %12.2e \n'%(optimize.cpt_iter, optimize.eta))
            f.write('-------------------------------------------------------------------------------------------------\n')
            f.write('  Iter_CG        qk       norm_res      norm_res/||gk||  \n')
            f.write('%8d %12.2e %12.2e %12.2e \n'%(optimize.cpt_iter_CG, optimize.qk_CG, optimize.norm_residual, optimize.norm_residual/optimize.norm_grad))
        else:
            f.write('%8d %12.2e %12.2e %12.2e \n'%(optimize.cpt_iter_CG, optimize.qk_CG, optimize.norm_residual, optimize.norm_residual/optimize.norm_grad))
        f.close()            
    else:
        f = open('iterate_TRN.dat', 'a')
        f.write('%6d %12.2e %12.2e %12.2e %12.2e %8d %8d %12.2e %8d %8d \n'%(
        optimize.cpt_iter, fcost, optimize.norm_grad, fcost/optimize.f0, optimize.alpha, optimize.cpt_ls, optimize.cpt_iter_CG, optimize.eta, optimize.nfwd_pb, optimize.nhess))
        f.close()


def print_info_PTRN(optim, fcost, FLAG):

    if(FLAG=='INIT'):
        f = open('iterate_PTRN.dat', 'w')
        f.write('******************************************************************************************\n')
        f.write('                      PRECONDITIONED TRUNCATED NEWTON ALGORITHM                           \n')
        f.write('******************************************************************************************\n')
        f.write('     Convergence criterion  : %10.2e \n'%optim.conv)
        f.write('     Niter_max              : %7d \n'%optim.niter_max)
        f.write('     Initial cost is        : %10.2e \n'%optim.f0)
        f.write('     Initial norm_grad is   : %10.2e \n'%optim.norm_grad)
        f.write('     Maximum CG iter        : %7d \n'%optim.niter_max_CG)
        f.write('******************************************************************************************\n')
        f.write('   Niter      fk         ||gk||       fk/f0          alpha      nls     nit_CG      eta       ngrad     nhess \n')
        f.write('%6d %12.2e %12.2e %12.2e %12.2e %8d %8d %12.2e %8d %8d \n'%(
        optim.cpt_iter, fcost, optim.norm_grad, fcost/optim.f0, optim.alpha, optim.cpt_ls, optim.cpt_iter_CG, optim.eta, optim.nfwd_pb, optim.nhess))
        f.close()

        f = open('iterate_PTRN_CG.dat', 'w')
        f.write('******************************************************************************************\n')
        f.write('                        PRECONDITIONED TRUNCATED NEWTON ALGORITHM                         \n')
        f.write('                                      INNER CG HISTORY                                    \n')
        f.write('******************************************************************************************\n')
        f.write('     Convergence criterion  : %10.2e \n'%optim.conv)
        f.write('     Niter_max              : %7d \n'%optim.niter_max)
        f.write('     Initial cost is        : %10.2e \n'%optim.f0)
        f.write('     Initial norm_grad is   : %10.2e \n'%optim.norm_grad)
        f.write('     Maximum CG iter        : %7d \n'%optim.niter_max_CG)
        f.write('******************************************************************************************\n')
        f.close()

    elif(FLAG=='CONV'):
        f = open('iterate_PTRN.dat', 'a')
        f.write('**********************************************************************\n')
        if(optim.cpt_iter==optim.niter_max):
            f.write('  STOP: MAXIMUM NUMBER OF ITERATION REACHED    \n')
        else:
            f.write('  STOP: CONVERGENCE CRITERION SATISFIED        \n')
        f.write('**********************************************************************\n')
        f.close()

    elif(FLAG=='FAIL'):
        f = open('iterate_PTRN.dat', 'a')
        f.write('**********************************************************************\n')
        f.write('  STOP: LINESEARCH FAILURE    \n')
        f.write('**********************************************************************\n')
        f.close()     
    elif(optim.comm=='DES1'):
        f = open('iterate_PTRN_CG.dat', 'a')
        f.write('-------------------------------------------------------------------------------------------------\n')
        f.write(' NONLINEAR ITERATION %4d  ETA IS : %12.2e \n'%(optim.cpt_iter, optim.eta))
        f.write('-------------------------------------------------------------------------------------------------\n')
        f.write('  Iter_CG        qk       norm_res      norm_res/||gk||  \n')
        f.write('%8d %12.2e %12.2e %12.2e \n'%(optim.cpt_iter_CG, optim.qk_CG, optim.norm_residual, optim.norm_residual/optim.norm_grad))
        f.close()     

    elif(optim.comm=='DES2'):
        f = open('iterate_PTRN_CG.dat', 'a')
        f.write('%8d %12.2e %12.2e %12.2e \n'%(optim.cpt_iter_CG, optim.qk_CG, optim.norm_residual, optim.norm_residual/optim.norm_grad))
        f.close()     
    else:
        f = open('iterate_PTRN.dat', 'a')
        f.write('%6d %12.2e %12.2e %12.2e %12.2e %8d %8d %12.2e %8d %8d \n'%(
        optim.cpt_iter, fcost, optim.norm_grad, fcost/optim.f0, optim.alpha, optim.cpt_ls, optim.cpt_iter_CG, optim.eta, optim.nfwd_pb, optim.nhess))
        f.close()


def print_info_PTRN(optimize, fcost):
    
    if(optimize.FLAG=='INIT'):
        f = open('iterate_PTRN.dat', 'w')
        f.write('******************************************************************************************\n')
        f.write('                      PRECONDITIONED TRUNCATED NEWTON ALGORITHM                           \n')
        f.write('******************************************************************************************\n')
        f.write('     Convergence criterion  : %10.2e \n'%optimize.conv)
        f.write('     Niter_max              : %7d \n'   %optimize.niter_max)
        f.write('     Initial cost is        : %10.2e \n'%optimize.f0)
        f.write('     Initial norm_grad is   : %10.2e \n'%optimize.norm_grad)
        f.write('     Maximum CG iter        : %7d \n'   %optimize.niter_max_CG)
        f.write('******************************************************************************************\n')
        f.write('   Niter      fk         ||gk||       fk/f0          alpha      nls     nit_CG      eta       ngrad     nhess \n')
        f.write('%6d %12.2e %12.2e %12.2e %12.2e %8d %8d %12.2e %8d %8d \n'%(
        optimize.cpt_iter, fcost, optimize.norm_grad, fcost/optimize.f0, optimize.alpha, optimize.cpt_ls, optimize.cpt_iter_CG, optimize.eta, optimize.nfwd_pb, optimize.nhess))
        f.close()

        f = open('iterate_PTRN_CG.dat', 'w')
        f.write('******************************************************************************************\n')
        f.write('                        PRECONDITIONED TRUNCATED NEWTON ALGORITHM                         \n')
        f.write('                                      INNER CG HISTORY                                    \n')
        f.write('******************************************************************************************\n')
        f.write('     Convergence criterion  : %10.2e \n'%optimize.conv)
        f.write('     Niter_max              : %7d \n'   %optimize.niter_max)
        f.write('     Initial cost is        : %10.2e \n'%optimize.f0)
        f.write('     Initial norm_grad is   : %10.2e \n'%optimize.norm_grad)
        f.write('     Maximum CG iter        : %7d \n'   %optimize.niter_max_CG)
        f.write('******************************************************************************************\n')
        f.close()

    elif(optimize.FLAG=='CONV'):
        f = open('iterate_PTRN.dat', 'a')
        f.write('**********************************************************************\n')
        if(optimize.cpt_iter == optimize.niter_max):
            f.write('  STOP: MAXIMUM NUMBER OF ITERATION REACHED    \n')
        else:
            f.write('  STOP: CONVERGENCE CRITERION SATISFIED        \n')
        f.write('**********************************************************************\n')
        f.close()

    elif(optimize.FLAG=='FAIL'):
        f = open('iterate_PTRN.dat', 'a')
        f.write('**********************************************************************\n')
        f.write('  STOP: LINESEARCH FAILURE    \n')
        f.write('**********************************************************************\n')
        f.close()     
    elif(optimize.comm=='DES1'):
        f = open('iterate_PTRN_CG.dat', 'a')
        f.write('-------------------------------------------------------------------------------------------------\n')
        f.write(' NONLINEAR ITERATION %4d  ETA IS : %12.2e \n'%(optimize.cpt_iter, optimize.eta))
        f.write('-------------------------------------------------------------------------------------------------\n')
        f.write('  Iter_CG        qk       norm_res      norm_res/||gk||  \n')
        f.write('%8d %12.2e %12.2e %12.2e \n'%(optimize.cpt_iter_CG, optimize.qk_CG, optimize.norm_residual, optimize.norm_residual/optimize.norm_grad))
        f.close()     

    elif(optimize.comm=='DES2'):
        f = open('iterate_PTRN_CG.dat', 'a')
        f.write('%8d %12.2e %12.2e %12.2e \n'%(optimize.cpt_iter_CG, optimize.qk_CG, optimize.norm_residual, optimize.norm_residual/optimize.norm_grad))
        f.close()     
    else:
        f = open('iterate_PTRN.dat', 'a')
        f.write('%6d %12.2e %12.2e %12.2e %12.2e %8d %8d %12.2e %8d %8d \n'%(
        optimize.cpt_iter, fcost, optimize.norm_grad, fcost/optimize.f0, optimize.alpha, optimize.cpt_ls, optimize.cpt_iter_CG, optimize.eta, optimize.nfwd_pb, optimize.nhess))
        f.close()
