% Some tests to assess the validity of the function exp_se2

% test 1

clear 

A = generate_se2_dA([0,1,1]);

exp_A      = exp_se2(A);
exp_ground = eye(3) + A;

risp = isequal(exp_ground, exp_A);

if risp == 0
    disp('test_se2 1 NOT passed');
elseif risp == 1
    disp('test_se2 1 passed');
end

% test 2

clear

A = generate_se2_dA();

exp_A  = exp_se2(A);
expm_A = expm(A);

decimal = 3;

risp = isequal(round(expm_A, decimal), round(exp_A, decimal));

if risp == 0
    disp('test_se2 2 not passed');
elseif risp == 1
    disp('test_se2 2 passed');
end

% test 3

clear

A = generate_se2_dA([.0, 3, 2]);

exp_A      = exp_se2(A);
exp_ground = eye(3);
exp_ground(1,3) = 3; 
exp_ground(2,3) = 2;

risp = isequal(exp_ground, exp_A);

if risp == 0
    disp('test_se2 3 not passed');
elseif risp == 1
    disp('test_se2 3 passed');
end

clear
