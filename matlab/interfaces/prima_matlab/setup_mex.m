% Compile Fortran MEX files for full numerical consistency
original_dir = pwd;
base_dir = fileparts(mfilename('fullpath'));
mex_dir = fullfile(base_dir, '+fortran', 'private');
if ~isfolder(mex_dir); mkdir(mex_dir); end
cd(mex_dir);

% Extra compiler flags matching compile.m (gfortran on Linux)
extra_fflags = ['FFLAGS="$FFLAGS', ...
    ' -Wno-missing-include-dirs', ...
    ' -fno-stack-arrays -frecursive', ...
    ' -fbacktrace -fno-stack-protector'];
% Conditionally add -ftrampoline-impl=heap (gfortran >= 14 + libgcc >= 14)
% matching the logic in setup_tools/compile.m
try
    compiler_info = mex.getCompilerConfigurations('fortran', 'selected');
    if contains(lower(compiler_info.Manufacturer), 'gnu')
        compiler_major = sscanf(compiler_info.Version, '%d');
        if ~isempty(compiler_major) && compiler_major >= 14
            extra_fflags = [extra_fflags, ' -ftrampoline-impl=heap'];
        end
    end
catch
end
extra_fflags = [extra_fflags, '"'];
mex('-R2018a', 'FOPTIMFLAGS=-O2', extra_fflags, '../fortran_cos.F90');
mex('-R2018a', 'FOPTIMFLAGS=-O2', extra_fflags, '../fortran_exp.F90');
mex('-R2018a', 'FOPTIMFLAGS=-O2', extra_fflags, '../fortran_log.F90');
mex('-R2018a', 'FOPTIMFLAGS=-O2', extra_fflags, '../fortran_power_integer.F90');
mex('-R2018a', 'FOPTIMFLAGS=-O2', extra_fflags, '../fortran_power_real.F90');
mex('-R2018a', 'FOPTIMFLAGS=-O2', extra_fflags, '../fortran_sin.F90');
mex('-R2018a', 'FOPTIMFLAGS=-O2', extra_fflags, '../fortran_sqrt.F90');
mex('-R2018a', 'FOPTIMFLAGS=-O2', extra_fflags, '../fortran_tan.F90');
cd(original_dir);
