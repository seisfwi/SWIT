import numpy as np
from rosenbrock import rosenbrock, rosenbrock_hess
from optimization import Optimization

#----------------------------------------------------#
# Optimize class                                     #
#----------------------------------------------------#
optim = Optimization(niter_max = 10000, conv = 1e-8, method = 'TRN')
    
#----------------------------------------------------#
# intial guess                                       #
#----------------------------------------------------#
n = 2
x = np.zeros(n)
grad = np.zeros(n)
grad_preco = np.zeros(n)

x[0] = 1.5
x[1] = 1.5

#----------------------------------------------------#
# computation of the cost and gradient associated    #
# with the initial guess                             #
#----------------------------------------------------#
fcost, grad = rosenbrock(x)

#----------------------------------------------------#
# optimization loop: while convergence not reached or#
# linesearch not failed, iterate                     #
#----------------------------------------------------#
while ((optim.FLAG != 'CONV') and (optim.FLAG != 'FAIL')):
    x = optim.iterate(x, fcost, grad, grad_preco)

    if(optim.FLAG == 'GRAD'):
        # if FLAG is GRAD, then compute cost and gradient in  fcost, grad
        fcost, grad = rosenbrock(x)

    elif(optim.FLAG == 'HESS'):
        # if FLAG is HESS, then multiply optim.d by the Hessian operator and store the result in optim.Hd
        optim.Hd = rosenbrock_hess(x,optim.d)
    

print('END OF TEST')
print('FINAL iterate is : ', x[:])
print('See the convergence history in iterate_TRN.dat and iterate_TRN_CG.dat')



