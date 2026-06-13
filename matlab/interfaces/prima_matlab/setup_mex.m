% Compile Fortran MEX files for full numerical consistency
original_dir = pwd;
base_dir = fileparts(mfilename('fullpath'));
mex_dir = fullfile(base_dir, '+fortran', 'private');
if ~isfolder(mex_dir); mkdir(mex_dir); end
cd(mex_dir);
mex -R2018a FOPTIMFLAGS='-O3' '../fortran_cos.F90'
mex -R2018a FOPTIMFLAGS='-O3' '../fortran_exp.F90'
mex -R2018a FOPTIMFLAGS='-O3' '../fortran_log.F90'
mex -R2018a FOPTIMFLAGS='-O3' '../fortran_power_integer.F90'
mex -R2018a FOPTIMFLAGS='-O3' '../fortran_power_real.F90'
mex -R2018a FOPTIMFLAGS='-O3' '../fortran_sin.F90'
mex -R2018a FOPTIMFLAGS='-O3' '../fortran_sqrt.F90'
mex -R2018a FOPTIMFLAGS='-O3' '../fortran_tan.F90'
cd(original_dir);
