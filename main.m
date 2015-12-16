%{
----------------------------------------------------------------------
Aim of this procedure is to compare the expmv methods collected in the 
folder, for the 4 different classes of stationary linear ODE.
----------------------------------------------------------------------

%}

addpath(genpath(pwd))

clear
close all
clc

format long

% controller
show_comparison_1 = true;
show_comparison_2 = false;

% handle exp methods
exp_method = {@exp_se2, @expm, @exp_leja, @expmv, @expmvp, @phileja, @phipm};
exp_method_name = {'exp_se2', 'exp_leja', 'expmv', 'expmvp', 'phileja', 'phipm'};

if show_comparison_1 == true

    disp('First comparison with one matrix whose exponential has a known ground truth, given by the close form:')
    dv = [0.5, 0.5, 0]';
    v = [0.5, 0.5, 1]';
    dA = generate_se2_dA([pi/8, 0.1, 0.1]);
    
    % ground. Matrix of class 4 are the only for which 
    % we have a ground truth implemented.
    
    A = exp_se2(dA); 
    
    disp('Matrix in the algebra se2: dA = ')
    disp(dA)
    disp('Corresponding matrix in the group SE2: A = ')
    disp(A)
    disp('---')
    disp('Vector v')
    disp(v)
    disp('Product A times v, ground truth')
    disp(A*v)
    disp('Products expm(dA) * v')
    disp('With expm:')
    disp(expm(dA)*v)
    disp('With expleja:')
    disp(expleja(1, dA, v))
    disp('With expmv:')
    disp(expmv(1, dA, v, 'double'))
    disp('With expmvp:')
    disp(expmvp(1, dA, v))
    disp('With phileja:')
    disp(phileja(1, dA, 1, v))
    disp('With phipm:')
    disp(phipm(1, dA, v))
    
    
    
end

if show_comparison_2 == true
    N = 500;
    disp('Second comparison with multiple random matrices of any class.')
    disp('The ground truth is given by expm(dA)*v .')
    
end



disp('Second example with :')





