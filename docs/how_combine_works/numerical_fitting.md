# Overview

Combine provides the CascadeMinimizer, which is essentially a wrapper around the RooMinimizer, which itself is essentially a wrapper around a minimizer, typically either MINUIT or MINUIT2.

These are numerical minimizers, which provide several numerical minimization routines.


## MIGRAD

This is the default minimzation routine used by MINUIT, and it is a quasi-newtonian based minimization method. 
i.e. it relies on the functions gradient and second-derivative (hessian) matrix, but uses an approximation scheme for estimating the second derivative matrix at each point, rather than try to calculat it at each step, which is very costly.
It uses the current estimate of the covariance matrix to determine its search direction, and after each step the estimate of the covariance matrix is updated.

Essentially, from an initial point, a parabolic function can be minimized analytically. 
A newtonian method approximates the function as parabolic and finds the analytic solution.
This analytic solution to the parabolic approximation is used as the estimate.
At this estimated point, the function is again approximated as parabolic and a new estimated minimum is found. 
This processes is iterated until the true minimum is found (to within some accuracy).

This netwonian method suffers from very expensive calculations of the inverse second-derivative matrix (inverse hessian matrix) for large models (order n^3 for model with n parameters, inverting the matrix [complexity of gaussian elimination]).
For a quasi-newtonian method, the approximate inverse hessian is used -- various methods exist depending on the method used for the estimation.

Thus, as the end of the MIGRAD routine, there is an existing estimate of the covariance matrix from which the uncertainties could be derived.
However, this estimate may suffer from various precision issues, based on the non-linearity of the problem and the steps taken in finding the minimum.

Therefore, finding the uncertainty errors is best done an explicit call to an error finding routine after the minimization.

In fact, the MINUIT algorithm can be 'safer' by making more explicit calculations (Actually, still approximations, but better ones) of the second derivative matrix through calls to HESSE.
The MINUIT "strategy" can be set to values of either 0, 1, or 2. Corresponding to increasing numbers of calls to HESSE, (0 being "unsafe" with few calls, and 2 being "safe" with the most calls).


The two most often used methods are HESSE and MINOS, which also have important differences. 
In particular HESSE on makes symmetric error estimates and assumes a parabolic functional form, whereas MINOS handles asymmetric errors and non-parabolic functions.

## HESSE

HESSE is a method for calculating uncertaintes based on the full second-derivative matrix of the function.
Once the second derivatives are calculated, the matrix is inverted.

If HESSE is called at the true minimum, then the parameter errors for any real problem should positive definite.
However, sometimes non-positive definite values can be returned

## MINOS

MINOS is a method for calculating parameter errors. 
Rather than using the second derivative at the minimum, minos finds the crossing where the function being minimized is at least some value DELTA greater than the minimum.
In the case that the function is parabolic, this will correspond exactly to what is returned by HESSE. 
For the general case, however, Minos provides more complete information, because it depends on the full functional form rather than just the second derivatives at the minimum.


