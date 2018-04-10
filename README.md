Matlab 2016a


# matrix_expmv: exponential map of the square matrix dA times a vector


Research code aimed at comparing a set of numerical methods to compare the exponential of a matrix times a vector without explicitly compute the exponential of the matrix. Selected methods are the state of the art, and they are applied at 6 different classes (or tastes) of stationary linear ordinary differential equations.

## Where to start

Run the procedure `test/main_test` to run all the test in one go and reproduce the image below:

![run_example](https://github.com/SebastianoF/matrix_expmv/blob/master/screenshots/test_output.jpg)

## Not(at)ions

### 2D case

Matrices are defined by a 2x2 Stationary linear ODE system

dx/dt = a x + b y + alpha
dy/dt = a x + b y + beta

where 

dA = [a, b, alpha;
      c, d, beta;
      0, 0, 0]

is the associated matrix in homogeneous coordinates.
The exponential map of dA, among other things, solves the ode.

The matrices of the kind dA can be classifed according to 6 types (called
"tastes" here) according to the 2 eigenvalues of dA

+ **Taste 1:** real eigenvalues with the same signs, positive (unstable node).

+ **Taste 2:** real eigenvalues with the same signs, negative (stable node).

+ **Taste 3:** real eigenvalues with opposite signs (saddle).

+ **Taste 4:** complex (conjugates) eigenvalues with positive real part (unstable spiral).

+ **Taste 5:** complex (conjugates) eigenvalues with negative real part (stable spiral).

+ **Taste 6:** complex (conjugates) eigenvalues with zero real part (circles). These are elements of the Lie algebra se2


The function `generate_random` generates one of these, given the class, and 
a range of values where we want the rotation (theta) and the translation
(tx, ty) to belong to.

The function `generate_se2`, generates elements in the particular case of 
taste 6 which are elements of the Lie algebra se2.

### N-D case

For problems with more than 2-dimensions we can not classify matrices with 
the previous 6 tastes straightforwardly. We used instead a classification 
with three 'string' tastes: 

+ **Taste 'neg':** If all the eigenvalues, both real or complex, have negative or zero real 
part. 

+ **Taste 'pos':** If all the eigenvalues, both real or complex, have positive real part. 

+ **Taste 'mix':** If all the eigenvalues, both real or complex, have bot positive and 
negative real part.


## Methods compared:

+ [exp_leja](http://uk.mathworks.com/matlabcentral/fileexchange/44039-matrix-exponential-times-a-vector/content/expleja.m)

+ [expmv](http://www.mathworks.com/matlabcentral/fileexchange/29576-matrix-exponential-times-a-vector/content/expmv.m)

+ [expmvp](http://www1.maths.leeds.ac.uk/~jitse/expmvp.m) \* 

+ [phileja](http://uk.mathworks.com/matlabcentral/fileexchange/40949-meshfree-exponential-integrator/content/MExpInt2D/phileja.m)

+ [phipm](http://www1.maths.leeds.ac.uk/~jitse/phipm.m) \*

m-files have been copy-pasted in the `expmv_method` folder, in date 13-Dec-2015.

(\*) Little modifications from the originally downloaded code to make it version compatible.

## Acknowledgments

The code was developed within the [GIFT-surg research project](http://www.gift-surg.ac.uk) UCL (UK), with the aim of prototyping the exponential of large stationary velocity fields arising in diffeomorphic medical image registration algorithms.

## Bibliography:

+ For a classification of the linear stationary ode:
[Hoppenstead  "Analysis and simulation of chaotic system" (ch. 2.1)](http://www.springer.com/gb/book/9780387989433)

+ For the theoretical comparison of some of the selected methods:
[Caliari et al  "Comparison of various methods for computing the action of
the matrix exponential"](http://profs.scienze.univr.it/~caliari/pdf/preCKOR13.pdf)

+ Milestone in the computation of the matrix exponential:
[Moler, Van Loan
"Nineteen dubious ways to compute the exponential of a matrix"](http://www.cs.cornell.edu/cv/researchpdf/19ways+.pdf)


