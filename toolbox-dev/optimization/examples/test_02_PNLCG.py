import numpy as np
from rosenbrock import rosenbrock
from optimization import Optimization

#----------------------------------------------------#
# Optimize class                                     #
#----------------------------------------------------#
optim = Optimization(niter_max = 1000000, conv = 1e-8, method = 'CG')
    
#----------------------------------------------------#
# intial guess                                       #
#----------------------------------------------------#
n = 100
x = np.random.rand(n)
grad = np.zeros(n)
grad_preco = np.zeros(n)

x[:] = 1.1
print('Initial solution : ', x[:])
  
#----------------------------------------------------#
# computation of the cost and gradient associated    #
# with the initial guess                             #
#----------------------------------------------------#
fcost, grad = rosenbrock(x)
grad_preco[:] = grad[:] 

#----------------------------------------------------#
# optimization loop: while convergence not reached or#
# linesearch not failed, iterate                     #
#----------------------------------------------------#
while ((optim.FLAG != 'CONV') and (optim.FLAG != 'FAIL')):
    x = optim.iterate(x, fcost, grad, grad_preco)
    if(optim.FLAG == 'GRAD'):
        #compute cost and gradient at point x
        fcost, grad = rosenbrock(x)
        grad_preco  = np.copy(grad)
    
    
print('END OF TEST')
print('FINAL iterate is : ', x[:])
print('See the convergence history in iterate_CG.dat')




