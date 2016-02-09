%{
----------------------------------------------------------------------
Aim of this procedure is to compare the expmv methods collected in the 
folder, for the 3 possible tastes defined for large squared matrices.
----------------------------------------------------------------------
In case of large matrices we expect to reduce the computational time
of the computation of expm(dA)*v using the expmv methods instead, 
maintaing a reasonable accuracy.

Note 1: for matrices bigger than 90x90 the method expmvp has a bug.

Note 2: if the input matrices are not big enough, it is faster to compute
directly expm(dA)*v than the same with one of the selected methods.

%}

clear
close all
clc

%%%%%%%%%%%%%%%%%%
%%% controller %%%
%%%%%%%%%%%%%%%%%%

compute = 0;  % if 1 compute and save, if 0 load the data.

%%%%%%%%%%%%%%%%%%%
%%% model- view %%%
%%%%%%%%%%%%%%%%%%%

n = 200;       % dimension of the square matrix
S = 50;       % Number of samples for each taste.
Methods = 5;  % expm, exp_leja, expmv, expmvp, phileja, phipm.
Tastes = 3;   % taste 'pos', 'neg', 'mix'
tastes_names = {'pos', 'neg', 'mix'};

if n < 90
    exp_method_name = {'exp_leja', 'expmv', 'expmvp', 'phileja', 'phipm'};
    Methods = 5;
else
    exp_method_name = {'exp_leja', 'expmv', 'phileja', 'phipm'};
    Methods = 4;
end

%%% to avoid the computaiton at each run:

if compute == 1
    
    v = 10*rand(n,1);

    Errors_n_dim = zeros(Methods, Tastes, S);
    Times_n_dim  = zeros(Methods, Tastes, S);
    ground_time_val_n_dim = zeros(Tastes, S);

    fprintf('\n\n')
    disp('Comparison with multiple random matrices of any taste.')
    disp(strcat('Matrices have dimension', num2str(n), 'times',  num2str(n)));
    disp('The ground truth is given by expm(dA)*v .')

    % 1 colour for each taste
    % 1 plot for each method
    
    % waitbar with cancel button creation
    h = waitbar(0, 'steps');
    

    for s = 1:S
        for ta = 1:Tastes
            % generate dA and compute the ground truth expm(dA)*v

            dA = generate_rand_dA_by_taste_n_dim(n, tastes_names(ta));  
            tic
            A_v_gr = expm(dA)*v;
            ground_time_val_n_dim(ta, s) = toc;


            % - this part can be refactored using handles - 
            disp('')

            % Method 1: expleja
            tic
            A_v_expleja = expleja(1, dA, v);
            Times_n_dim(1, ta, s) = toc;
            Errors_n_dim(1, ta, s) = norm(A_v_gr - A_v_expleja);

            % Method 2: expmv
            tic
            A_v_expmv   = expmv(1, dA, v, 'double');
            Times_n_dim(2, ta, s) = toc;
            Errors_n_dim(2, ta, s) = norm(A_v_gr - A_v_expmv);
            
            % Method 3: phileja
            tic
            A_v_phileja = phileja(1, dA, 0, v);
            Times_n_dim(3, ta, s) = toc;
            Errors_n_dim(3, ta, s) = norm(A_v_gr - A_v_phileja);

            % Method 4: phipm
            tic
            A_v_phipm   = phipm(1, dA, v);
            Times_n_dim(4, ta, s) = toc;
            Errors_n_dim(4, ta, s) = norm(A_v_gr - A_v_phipm);

            % Method 5: expmvp
            if n < 90  % infinite loop if expmvp is used for n > 90 or 100
                tic
                A_v_expmvp  = expmvp(1, dA, v);
                Times_n_dim(5, ta, s) = toc;
                Errors_n_dim(5, ta, s) = norm(A_v_gr - A_v_expmvp);
            end
            
        end
        % waitbar update
        waitbar(s/S, h, sprintf('element %d out of %d',s, S))
    end
    
    % waitbar close
    close(h)
   
    save('results/Errors_n_dim.mat', 'Errors_n_dim')
    save('results/Times_n_dim.mat', 'Times_n_dim')
    save('results/ground_time_val_n_dim.mat', 'ground_time_val_n_dim')
    disp('Data saved in folder results')
    
else
    load('results/Errors_n_dim.mat', 'Errors_n_dim')
    load('results/Times_n_dim.mat', 'Times_n_dim')
    load('results/ground_time_val_n_dim.mat', 'ground_time_val_n_dim')
    disp('Data loaded from folder results')

end

%%%%%%%%%%%%
%%% view %%%
%%%%%%%%%%%%


figure('units','normalized','position',[.1 .1 .8 .3]); 
for m = 1:Methods
    subplot(1,Methods,m)
    hold on
    for ta = 1:Tastes
        scatter(Errors_n_dim(m, ta, :), Times_n_dim(m, ta, :));
        set(gca,'xscale','log')
        %set(gca,'yscale','log')
        title(exp_method_name(m))
        xlabel('Error (norm2)')
        ylabel('Computational time (sec)') 
    end
    hold off

end

h = legend(gca, ...
           strcat('Taste ', ' pos'), ...
           strcat('Taste ', ' neg'), ...
           strcat('Taste ', ' mix'), ...
               'Location','northeastoutside');
    
pos = get(h,'position');
set(h, 'position',[0.9198 0.7259 pos(3:4)])

%%%%%%%%%%%%%%%%%%%

figure('units','normalized','position',[.1 .1 .8 .3]);
    
for m = 1:Methods
    subplot(1,Methods,m)
    hold on
    for ta = 1:Tastes
        % time is given by the sum of each sampling for each taste,
        % divided by the number of sampling times the computational
        % time of the ground truth expm(dA)*v
        scatter(Errors_n_dim(m, ta, :), (sum(Times_n_dim(m, ta, :))/...
            (S*mean(ground_time_val_n_dim(ta, :)))) * ones(size(Errors_n_dim(m, ta, :))) );
        set(gca,'xscale','log')
        %set(gca,'yscale','log')
        title(exp_method_name(m))
        xlabel('Error (norm2)')
        ylabel('Relative mean computational time (sec)') 
    end
    hold off

end

h = legend(gca, ...
           strcat('Taste ', ' pos'), ...
           strcat('Taste ', ' neg'), ...
           strcat('Taste ', ' mix'), ...
               'Location','northeastoutside');
    
pos = get(h,'position');
set(h, 'position',[0.9198 0.7259 pos(3:4)])

