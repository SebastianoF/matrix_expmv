MAtlab 2014b (.11)

=======================
Exponential map of dA times a vector
=======================

----------------------------------------------------------------------
Aim of this project is to compare some expmv methods collected in the 
folder expmv_methods, for the 4 different classes of stationary linear 
Ordinary differential equations.
----------------------------------------------------------------------

Matrices are defined by a 2x2 Stationary linear ODE system

dx/dt = a x + b y + alpha
dy/dt = a x + b y + beta

where 

dA = [a, b, alpha;
      c, d, beta;
      0, 0, 0]

is the associated matrix.
The exponential map of dA, among other things, solves the ode.

The matrices of the kind dA can be classifed according to 4 classes
according to the eigenvalues of the part

dAs = [a, b; 
       c, d]

Class 1:
real eigenvalues with the same signs (node)

Class 2:
real eigenvalues with opposite signs (saddle)

Class 3:
complex (conjugates) eigenvalues with negative real part (spiral)

Class 4:
complex (conjugates) eigenvalues with positive real part (circles)


See test_generator_class_division to see the four classes of ode

The function generate_random generates one of these, given the class, and 
a range of values where we want the rotation (theta) and the translation
(tx, ty) to belong to.

The function generate_se2, generates elements in the particular case of 
class 4 which are also elements of the Lie algebra se2.


----

Hoppenstead "Analysis and simulation of chaotic system" (ch. 2.1)
For a classification of the linear stationary ode.

Caliari et al "Comparison of various methods for computing the action of
the matrix exponential"
For the theoretical comparison of some of the selected methods.

-----

Please run the test to see that everything works!
