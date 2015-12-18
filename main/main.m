%{
----------------------------------------------------------------------
Aim of this procedure is to compare the expmv methods collected in the 
folder, for the 6 different tastes* of 2d stationary linear ODE.
----------------------------------------------------------------------
*) see README.dm
%}

clear
close all
clc

%%%%%%%%%%%%%%%%%%
%%% controller %%%
%%%%%%%%%%%%%%%%%%

show_comparison_1 = true;
show_comparison_2 = true;

compute = 0;  % if 1 compute and save, if 0 load the data.

% handles exp methods
exp_method = {@exp_se2, @expm, @exp_leja, @expmv, @expmvp, @phileja, @phipm};
exp_method_name = {'exp_leja', 'expmv', 'expmvp', 'phileja', 'phipm'};

%%%%%%%%%%%%%%%%%%%
%%% model- view %%%
%%%%%%%%%%%%%%%%%%%

if show_comparison_1 == true

    disp('First comparison with one matrix whose exponential has a known')
    disp('ground truth, given by the close form:')
    dv = [0.5, 0.5, 0]';
    v = [0.5, 0.5, 1]';
    dA = generate_se2_dA([pi/8, 0.1, 0.1]);
    
    % ground. Matrix of taste 6 are the only for which 
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
   
    fprintf('\n Error for each method, taste 6 (ground truth available)')
          
    fprintf('\n\n %17s %17s %17s %17s %17s %17s\n', ...
            'expm', 'exp_leja', 'expmv', 'expmvp', 'phileja', 'phipm');
        
    fprintf(' %17.8f %17.8f %17.8f %17.8f %17.8f %17.8f \n',...
            errors');
        
    format short
    
end

if show_comparison_2 == true
        
    S = 80;       % Number of samples for each taste.
    Methods = 5;  % expm, exp_leja, expmv, expmvp, phileja, phipm.
    Tastes = 6;  % taste 1, 2, 3, 4, 5, 6
    
    %%% to avoid the computaiton at each run
    if compute == 1

        v = [10.5, 10.5, 1]';

        Errors = zeros(Methods, Tastes, S);
        Times  = zeros(Methods, Tastes, S);
        ground_time_val = zeros(Tastes, S);
        
        fprintf('\n\n')
        disp('Second comparison with multiple random matrices of any taste.')
        disp('The ground truth is given by expm(dA)*v .')

        % 1 colour for each taste
        % 1 plot for each method
        
        % waitbar with cancel button creation
        h = waitbar(0, 'steps');

        for s =1:S
            for ta = 1:Tastes
                % generate dA and compute the ground truth expm(dA)*v

                dA = generate_rand_dA_by_taste(ta);  
                tic
                A_v_gr = expm(dA)*v;
                ground_time_val(ta, s) = toc;


                % - this part can be refactored using handles - 
                disp('')

                % Method 1: expleja
                tic
                A_v_expleja = expleja(1, dA, v);
                Times(1, ta, s) = toc;
                Errors(1, ta, s) = norm(A_v_gr - A_v_expleja);

                % Method 2: expmv
                tic
                A_v_expmv   = expmv(1, dA, v, 'double');
                Times(2, ta, s) = toc;
                Errors(2, ta, s) = norm(A_v_gr - A_v_expmv);

                % Method 3: expmvp
                tic
                A_v_expmvp  = expmvp(1, dA, v);
                Times(3, ta, s) = toc;
                Errors(3, ta, s) = norm(A_v_gr - A_v_expmvp);

                % Method 4: phileja
                tic
                A_v_phileja = phileja(1, dA, 0, v);
                Times(4, ta, s) = toc;
                Errors(4, ta, s) = norm(A_v_gr - A_v_phileja);

                % Method 5: phipm
                tic
                A_v_phipm   = phipm(1, dA, v);
                Times(5, ta, s) = toc;
                Errors(5, ta, s) = norm(A_v_gr - A_v_phipm);

            end
            % waitbar update
            waitbar(s/S, h, sprintf('element %d out of %d',s, S))
        end
        
        % close waitbar
        close(h)
        
        save('results/Errors.mat', 'Errors')
        save('results/Times.mat', 'Times')
        save('results/ground_time_val.mat', 'ground_time_val')
        disp('Data saved in folder results')
        
    else        
        load('results/Errors.mat', 'Errors')
        load('results/Times.mat', 'Times')
        load('results/ground_time_val.mat', 'ground_time_val')
        disp('Data loaded from folder results')
        
    end
    
    %%%%%%%%%%%%
    %%% view %%%
    %%%%%%%%%%%%
    
    figure('units','normalized','position',[.1 .1 .8 .3]);
    
    for m = 1:Methods
        subplot(1,5,m)
        hold on
        for ta = 1:Tastes
            scatter(Errors(m, ta, :), Times(m, ta, :));
            set(gca,'xscale','log')
            %set(gca,'yscale','log')
            title(exp_method_name(m))
            xlabel('Error (norm2)')
            ylabel('Computational time (sec)') 
        end
        hold off
        
    end
    
    h = legend(gca, ...
               strcat('Taste ', num2str(1)), ...
               strcat('Taste ', num2str(2)), ...
               strcat('Taste ', num2str(3)), ...
               strcat('Taste ', num2str(4)), ...
               strcat('Taste ', num2str(5)), ...
               strcat('Taste ', num2str(6)), ...
               'Location','northeastoutside');
    
    pos = get(h,'position');
    set(h, 'position',[0.9198 0.7259 pos(3:4)])
    
    
    %%%%%%%%%%%%%%%%%
    
    figure('units','normalized','position',[.1 .1 .8 .3]);
    
    for m = 1:Methods
        subplot(1,5,m)
        hold on
        for ta = 1:Tastes
            % time is given by the sum of each sampling for each taste,
            % divided by the number of sampling
            scatter(Errors(m, ta, :), (sum(Times(m, ta, :))/S) * ones(size(Times(m, ta, :))) );
            set(gca,'xscale','log')
            %set(gca,'yscale','log')
            title(exp_method_name(m))
            xlabel('Error (norm2)')
            ylabel('Mean computational time (sec)') 
        end
        hold off
 
    end
    
    h = legend(gca, ...
               strcat('Taste ', num2str(1)), ...
               strcat('Taste ', num2str(2)), ...
               strcat('Taste ', num2str(3)), ...
               strcat('Taste ', num2str(4)), ...
               strcat('Taste ', num2str(5)), ...
               strcat('Taste ', num2str(6)), ...
               'Location','northeastoutside');
    
    pos = get(h,'position');
    set(h, 'position',[0.9198 0.7259 pos(3:4)])
    
    
    %%%%%%%%%%%%%%%%
    
    figure('units','normalized','position',[.1 .1 .8 .3])
    for m = 1:Methods
        subplot(1,5,m)
        hold on
        for ta = 1:Tastes

            scatter(mean(Errors(m, ta, :)), mean(Times(m, ta, :)));
            set(gca,'xscale','log')
            %set(gca,'yscale','log')
            title(exp_method_name(m))
            xlabel('Mean error (norm2)')
            ylabel('Mean computational time (sec)') 
        end
        hold off

    end
     
     h = legend(gca, ...
               strcat('Taste ', num2str(1)), ...
               strcat('Taste ', num2str(2)), ...
               strcat('Taste ', num2str(3)), ...
               strcat('Taste ', num2str(4)), ...
               strcat('Taste ', num2str(5)), ...
               strcat('Taste ', num2str(6)), ...
               'Location','northeastoutside');

    pos = get(h,'position');
    set(h, 'position',[0.9198 0.7259 pos(3:4)])
    
    %%%%%%%%%%%%%%%%%%%%
    
    figure('units','normalized','position',[.1 .1 .8 .3]);
    
    for m = 1:Methods
        subplot(1,5,m)
        hold on
        for ta = 1:Tastes
            % time is given by the sum of each sampling for each taste,
            % divided by the number of sampling times the computational
            % time of the ground truth expm(dA)*v
            scatter(Errors(m, ta, :), (sum(Times(m, ta, :))/(S*mean(ground_time_val(ta, :)))) * ones(size(Errors(m, ta, :))) );
            set(gca,'xscale','log')
            %set(gca,'yscale','log')
            title(exp_method_name(m))
            xlabel('Error (norm2)')
            ylabel('Relative mean computational time (sec)') 
        end
        hold off
        
    end
    
    h = legend(gca, ...
               strcat('Taste ', num2str(1)), ...
               strcat('Taste ', num2str(2)), ...
               strcat('Taste ', num2str(3)), ...
               strcat('Taste ', num2str(4)), ...
               strcat('Taste ', num2str(5)), ...
               strcat('Taste ', num2str(6)), ...
               'Location','northeastoutside');
    
    pos = get(h,'position');
    set(h, 'position',[0.9198 0.7259 pos(3:4)])
    
end


disp('END')



