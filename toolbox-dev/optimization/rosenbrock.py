import numpy as np

def rosenbrock(x):
    ''' n-dimension Rosenbrock function
    '''
    
    n = np.size(x)
    fcost = np.sum(100*(x[:-1]**2. - x[1:])**2. + (x[:-1] - 1.)**2.)
    grad = np.zeros(n)
    grad[1:-1] = -200.*(x[:-2]**2. - x[1:-1]) + 400.*x[1:-1]*(x[1:-1]**2. - x[2:]) + 2.*(x[1:-1]-1.)
    grad[0]    =  400.*x[0]*(x[0]**2. - x[1]) + 2.*(x[0] - 1)
    grad[-1]   = -200.*(x[-2]**2. - x[-1])

    return fcost, grad


def rosenbrock_hess(x, d):
    '''
    #---------------------------------------------!
    #  The routine Rosenbrock_Hess returns        !
    #  Hessian-vector product H(x)d in output Hd  !
    #  for input parameters x and d               !
    #  H is the Hessian matrix                    !
    #  x=(x1,x2), d=(d1,d2) are two vector of R^2 !
    #---------------------------------------------!
    '''

    if(np.size(x) > 2):
        raise ValueError('length of x <=2')

    Hd = np.zeros(2)
    Hd[0]=(1200.*x[0]**2-400.*x[1]+2.)*d[0]-400*x[0]*d[1]
    Hd[1]=-400.*x[0]*d[0] + 200.*d[1]

    return Hd