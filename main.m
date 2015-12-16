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

%%%%%%%%%%%%%%%%%%
%%% controller %%%
%%%%%%%%%%%%%%%%%%

show_comparison_1 = false;
show_comparison_2 = true;

% handles exp methods
exp_method = {@exp_se2, @expm, @exp_leja, @expmv, @expmvp, @phileja, @phipm};
exp_method_name = {'expm', 'exp_leja', 'expmv', 'expmvp', 'phileja', 'phipm'};

%%%%%%%%%%%%%%%%%%%
%%% model- view %%%
%%%%%%%%%%%%%%%%%%%

if show_comparison_1 == true

    disp('First comparison with one matrix whose exponential has a known ground truth, given by the close form:')
    dv = [0.5, 0.5, 0]';
    v = [0.5, 0.5, 1]';
    dA = generate_se2_dA([pi/8, 0.1, 0.1]);
    
    % ground. Matrix of class 4 are the only for which 
    % we have a ground truth implemented.
    
    A = exp_se2(dA); 
    
    A_v_gr      = A*v;
    A_v_expm    = expm(dA)*v;
    A_v_expleja = expleja(1, dA, v);
    A_v_expmv   = expmv(1, dA, v, 'double');
    A_v_expmvp  = expmvp(1, dA, v);
    A_v_phileja = phileja(1, dA, 0, v);
    A_v_phipm   = phipm(1, dA, v);
    
    format short
    
    disp('Matrix in the algebra se2: dA = ')
    disp(dA)
    disp('Corresponding matrix in the group SE2: A = ')
    disp(A)
    disp('---')
    disp('Vector v')
    disp(v)
    disp('Product A times v, ground truth')
    disp(A_v_gr)
    disp('Products expm(dA) * v')
    disp('With expm:')
    disp(A_v_expm)
    disp('With expleja:')
    disp(A_v_expleja)
    disp('With expmv:')
    disp(A_v_expmv)
    disp('With expmvp:')
    disp(A_v_expmvp)
    disp('With phileja:')
    disp(A_v_phileja)
    disp('With phipm:')
    disp(A_v_phipm)
    
    format long
    
    % errors for each method:
    errors = [norm(A_v_gr - A_v_expm    , 2),    ...
              norm(A_v_gr - A_v_expleja , 2), ...
              norm(A_v_gr - A_v_expmv   , 2),   ...
              norm(A_v_gr - A_v_expmvp  , 2),  ...
              norm(A_v_gr - A_v_phileja , 2), ...
              norm(A_v_gr - A_v_phipm, 2)];
   
    fprintf('\n Error for each method: class 4 ground truth available:')
          
    fprintf('\n\n %17s %17s %17s %17s %17s %17s\n', ...
            'expm', 'exp_leja', 'expmv', 'expmvp', 'phileja', 'phipm');
        
    fprintf(' %17.8f %17.8f %17.8f %17.8f %17.8f %17.8f \n',...
            errors');
        
    format short
    
end

if show_comparison_2 == true
    N = 500;
    disp('Second comparison with multiple random matrices of any class.')
    disp('The ground truth is given by expm(dA)*v .')
    
    
    
end


disp('END')



