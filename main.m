%{
----------------------------------------------------------------------
Aim of this procedure is to compare the expmv methods collected in the 
folder, for the 4 different classes of stationary linear ODE.
----------------------------------------------------------------------

%}

clear
close all

% controller
show_comparison_1 = true;
show_comparison_2 = true;

% handle exp methods
exp_method = {@expm, @exp_se2, @exp_leja, @expmv, @expmvp, @phileja, @phipm};
exp_method_name = {'exp_se2', 'exp_leja', 'expmv', 'expmvp', 'phileja', 'phipm'};

if show_comparison_1 == true

    disp('First comparison with one matrix whose exponential has a known ground truth, given by the close form:')
    dv = [0.5, 0.5, 0]';
    v = [0.5, 0.5, 1]';
    dA = generate_rand_dA();
    A = exp_se2(dA); % ground
    
    disp('Matrix in the algebra se2: dA = ')
    disp(dA)
    disp('Corresponding matrix in the group SE2: A = ')
    disp(A)
    disp('')
    
    A*v
    
    expm(dA)*v
    
    expleja(1, dA, v)
    

    
    
end

if show_comparison_2 == true
    N = 500;
    disp('Second comparison with multiple random matrices of any class.')
    disp('The ground truth is given by expm(dA)*v .')
    
end



disp('Second example with :')





