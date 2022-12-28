## Optimization class: a self-contained optimization package

This code is reimplemented by Haipeng Li in Python based on the SEISCOPE 
optimization toolbox (https://seiscope2.osug.fr/SEISCOPE-OPTIMIZATION-TOOLBOX-274). 

The following nonlinear optimization algorithms are included:
1. Preconditioned Steepest Descent: PSTD
2. Preconditioned Nonlinear Conjugate Gradient: PNLCG
3. Quasi-Newton l-BFGS method: LBFGS
4. Quasi-Newton Preconditioned l-BFGS method: PLBFGS
5. Truncated Newton method: TRN
6. Preconditioned Truncated Newton method: PTRN


The original Fortran code is written by Ludovic Métivier and Romain Modrak: 
>  Métivier, L., & Brossier, R. (2016). The SEISCOPE optimization toolbox: A large-scale 
nonlinear optimization library based on reverse communication. *Geophysics*, *81*(2), F1-F15.

