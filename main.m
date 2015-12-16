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
    
    S = 100;       % Number of samples for each class.
    Methods = 5;  % expm, exp_leja, expmv, expmvp, phileja, phipm.
    Classes = 4;  % class 1, 2, 3, 4
    
    v = [10.5, 10.5, 1]';
    
    Errors = zeros(Methods, Classes, S);
    Times  = zeros(Methods, Classes, S);
    
    disp('Second comparison with multiple random matrices of any class.')
    disp('The ground truth is given by expm(dA)*v .')
    
    % 1 colour for each class
    % 1 plot for each method
    
    for s =1:S
        for c = 1:Classes
            % generate dA and compute the ground truth expm(dA)*v
            
            dA = generate_rand_dA_cl(c);  
            A_v_gr = expm(dA)*v;
            
            % - this part can be refactored using handles - 
            
            % Method 1: expleja
            tic
            A_v_expleja = expleja(1, dA, v);
            Times(1, c, s) = toc;
            Errors(1, c, s) = norm(A_v_gr - A_v_expleja);
            
            % Method 2: expmv
            tic
            A_v_expmv   = expmv(1, dA, v, 'double');
            Times(2, c, s) = toc;
            Errors(2, c, s) = norm(A_v_gr - A_v_expmv);
            
            % Method 3: expmvp
            tic
            A_v_expmvp  = expmvp(1, dA, v);
            Times(3, c, s) = toc;
            Errors(3, c, s) = norm(A_v_gr - A_v_expmvp);
            
            % Method 4: phileja
            tic
            A_v_phileja = phileja(1, dA, 0, v);
            Times(4, c, s) = toc;
            Errors(4, c, s) = norm(A_v_gr - A_v_phileja);
            
            % Method 5: phipm
            tic
            A_v_phipm   = phipm(1, dA, v);
            Times(5, c, s) = toc;
            Errors(5, c, s) = norm(A_v_gr - A_v_phipm);
       
       end
       
    end
    
    figure('units','normalized','position',[.1 .1 .8 .3]);
    
    for m = 1:Methods
        subplot(1,5,m)
        hold on
        for c = 1:Classes
            scatter(Errors(m, c, :), Times(m, c, :));
            set(gca,'xscale','log')
            %set(gca,'yscale','log')
            title(exp_method_name(m))
            xlabel('Error (norm2)')
            ylabel('Computational time (sec)') 
        end
        
        legend(gca, ...
               strcat('Class ', num2str(1)), ...
               strcat('Class ', num2str(2)), ...
               strcat('Class ', num2str(3)), ...
               strcat('Class ', num2str(4)), ...
               'Location','NorthEast');
 
        hold off
        
    end
    
end


disp('END')



