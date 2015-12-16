Matlab 2014b

=======================================================
Exponential map of the square matrix dA times a vector
=======================================================

----------------------------------------------------------------------
Aim of this project is to compare some expmv methods collected in the 
folder expmv_methods, for the 4 different classes of 
stationary linear ordinary differential equations.
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

The matrices of the kind dA can be classifed according to 5 types
according to the eigenvalues of the part

dAs = [a, b; 
       c, d]

Type 1:
real eigenvalues with the same signs, positive (unstable node).

Type 2:
real eigenvalues with the same signs, negative (stable node).

Type 3:
real eigenvalues with opposite signs (saddle).

Type 4:
complex (conjugates) eigenvalues with negative real part (spiral).

Type 5:
complex (conjugates) eigenvalues with positive real part (circles).


See test_generator_class_division to see the four classes of ode

The function generate_random generates one of these, given the class, and 
a range of values where we want the rotation (theta) and the translation
(tx, ty) to belong to.

The function generate_se2, generates elements in the particular case of 
class 4 which are also elements of the Lie algebra se2.


----
Bibliography:

Hoppenstead 
"Analysis and simulation of chaotic system" (ch. 2.1)
%For a classification of the linear stationary ode.

Caliari et al 
"Comparison of various methods for computing the action of
the matrix exponential"
%For the theoretical comparison of some of the selected methods.

Moler, Van Loan
"nineteen dubious ways to compute the exponential of a matrix"
% Milestone in the computation of the matrix exponential

----
Methods compared:

'exp_leja': 
http://uk.mathworks.com/matlabcentral/fileexchange/44039-matrix-exponential-times-a-vector/content/expleja.m

'expmv':
http://www.mathworks.com/matlabcentral/fileexchange/29576-matrix-exponential-times-a-vector/content/expmv.m

'expmvp': (little modifications from the original verions) 
http://www1.maths.leeds.ac.uk/~jitse/expmvp.m

'phileja': 
http://uk.mathworks.com/matlabcentral/fileexchange/40949-meshfree-exponential-integrator/content/MExpInt2D/phileja.m

'phipm': (little modifications from the original verions) 
http://www1.maths.leeds.ac.uk/~jitse/phipm.m

-----

Please run the procedure main_test to see that everything works!

-----
