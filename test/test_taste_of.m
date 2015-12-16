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


% test 3: to see if a matrix of a taste 3 really has this taste according 
% to taste_of_dA

clear 

dA = diag([2,-3,0]);
   
if taste_of_dA(dA) == 3
    disp('test_taste_of_dA 3 passed');
else
    disp('test_taste_of_dA 3 NOT passed');
end

% test 4: to see if a matrix of a taste 4 really has this taste according 
% to taste_of_dA

clear 

dA = diag([2 + 3i, 2 - 3i,0]);
   
if taste_of_dA(dA) == 4
    disp('test_taste_of_dA 4 passed');
else
    disp('test_taste_of_dA 4 NOT passed');
end

% test 5: to see if a matrix of a taste 5 really has this taste according 
% to taste_of_dA

clear 

dA = diag([-2 + 3i, -2 - 3i,0]);
   
if taste_of_dA(dA) == 5
    disp('test_taste_of_dA 5 passed');
else
    disp('test_taste_of_dA 5 NOT passed');
end


% test 6: to see if a matrix of a taste 5 really has this taste according 
% to taste_of_dA

clear 

dA = diag([+3i, -3i, 0]);
   
if taste_of_dA(dA) == 6
    disp('test_taste_of_dA 6 passed');
else
    disp('test_taste_of_dA 6 NOT passed');
end
