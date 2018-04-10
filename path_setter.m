
here = mfilename('fullpath');
root = fileparts(here);

cd(root)

addpath(pwd)
addpath(fullfile(pwd, 'expmv_methods'))
addpath(fullfile(pwd, 'test'))
addpath(fullfile(pwd, 'utils'))
addpath(fullfile(pwd, 'main'))
addpath(fullfile(pwd, 'homographies'))
addpath(fullfile(pwd, 'stats'))