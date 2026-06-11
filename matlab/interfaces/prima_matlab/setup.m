base_dir = fileparts(mfilename('fullpath'));
addpath(fullfile(base_dir, 'common'));
addpath(fullfile(base_dir, 'cobyla'));
addpath(fullfile(base_dir, 'bobyqa'));
addpath(fullfile(base_dir, 'lincoa'));
addpath(fullfile(base_dir, 'newuoa'));
addpath(fullfile(base_dir, 'uobyqa'));
addpath(fullfile(base_dir, 'examples/newuoa'));

run('setup_mex.m');
