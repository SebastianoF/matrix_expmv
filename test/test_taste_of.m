% test 1: to see if a matrix of a taste 1 really has this taste according 
% to taste_of_dA

clear 

dA = diag([2,3,0]);
   
if taste_of_dA(dA) == 1
    disp('test_taste_of_dA 1 passed');
else
    disp('test_taste_of_dA 1 NOT passed');
end

% test 2: to see if a matrix of a taste 2 really has this taste according 
% to taste_of_dA

clear 

dA = diag([-2,-3,0]);
   
if taste_of_dA(dA) == 2
    disp('test_taste_of_dA 2 passed');
else
    disp('test_taste_of_dA 2 NOT passed');
end


% test 3: to see if a matrix of a taste 2 really has this taste according 
% to taste_of_dA

clear 

dA = diag([2,-3,0]);
   
if taste_of_dA(dA) == 3
    disp('test_taste_of_dA 3 passed');
else
    disp('test_taste_of_dA 3 NOT passed');
end

% test 3: to see if a matrix of a class 3 in input_class is really of this
% class

clear 

rng(5);

ta_input = 3;

dA = generate_rand_dA_cl(ta_input);
S = dA(1:2, 1:2);

e = eig(S); 
l1 = e(1); 
l2 = e(2);

ok = 0;

if imag(l1) > (10^4)*eps
    l1 = real(l1);
    if abs(real(l1)) > (10^4)*eps 
        ok = 1;
    end
end
   
if ok == 0
    disp('test_generator_class_division 3 NOT passed');
elseif ok == 1
    disp('test_generator_class_division 3 passed');
end

% test 4: to see if a matrix of a class 3 in input_class is really of this
% class

clear

ta_input = 4;

dA = generate_rand_dA_cl(ta_input);
S = dA(1:2, 1:2);

e = eig(S); l1 = e(1); l2 = e(2);

ok = 0;

if imag(l1) > (10^4)*eps
    l1 = real(l1); l2 = real(l2);
    if abs(real(l1)) < (10^4)*eps 
        ok = 1;
    end
end
   
if ok == 0
    disp('test_generator_class_division 4 NOT passed');
elseif ok == 1
    disp('test_generator_class_division 4 passed');
end