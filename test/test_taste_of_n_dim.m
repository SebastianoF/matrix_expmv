% test 1: to see if taste of 2 dim case and n dim cases are coherent.

% taste 1 is 'pos'
% taste 2 is 'neg'
% taste 3 is 'mix'
% taste 4 is 'pos'
% taste 5 is 'neg'
% taste 6 is 'neg'

generate_rand_dA_by_taste(1)

assert(strcmp(taste_of_dA_n_dim(generate_rand_dA_by_taste(1)), 'pos'))
assert(strcmp(taste_of_dA_n_dim(generate_rand_dA_by_taste(2)), 'neg'))
assert(strcmp(taste_of_dA_n_dim(generate_rand_dA_by_taste(3)), 'mix'))
assert(strcmp(taste_of_dA_n_dim(generate_rand_dA_by_taste(4)), 'pos'))
assert(strcmp(taste_of_dA_n_dim(generate_rand_dA_by_taste(5)), 'neg'))
assert(strcmp(taste_of_dA_n_dim(generate_rand_dA_by_taste(6)), 'neg'))

disp('test_taste_of_n_dim test 1 passed')

% test 2: one closed form for each case

dA_pos = eye(10);
dA_neg = -eye(10);
dA_mix = eye(10); dA_mix(3,3) = -5;
dA_zer = zeros(10);

assert(strcmp(taste_of_dA_n_dim(dA_pos), 'pos'))
assert(strcmp(taste_of_dA_n_dim(dA_neg), 'neg'))
assert(strcmp(taste_of_dA_n_dim(dA_mix), 'mix'))
assert(strcmp(taste_of_dA_n_dim(dA_zer), 'neg'))

disp('test_taste_of_n_dim test 2 passed')
