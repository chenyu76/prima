function varargout = prima_matlab_call(solver_name, varargin)
%PRIMA_MATLAB_CALL dispatches to the pure MATLAB implementation of PRIMA solvers.
%
%   This function should only be called by uobyqa.m, newuoa.m, bobyqa.m,
%   lincoa.m, and cobyla.m when options.fortran is false.
%
%   Usage:
%   [x, fx, nf, xhist, fhist, exitflag] = prima_matlab_call('uobyqa', ...)
%   [x, fx, nf, xhist, fhist, exitflag] = prima_matlab_call('newuoa', ...)
%   [x, fx, nf, xhist, fhist, exitflag] = prima_matlab_call('bobyqa', ...)
%   [x, fx, constrviolation, nf, xhist, fhist, chist, exitflag] = prima_matlab_call('lincoa', ...)
%   [x, fx, constrviolation, nlconstr, nf, xhist, fhist, chist, nlchist, exitflag] = prima_matlab_call('cobyla', ...)
%
%   The pure MATLAB implementation resides in the prima_matlab/ subdirectory
%   of the interfaces/ directory. This function temporarily adds that
%   directory to the MATLAB search path.

% Attribute: private (not supposed to be called by users)

mfiledir = fileparts(fileparts(mfilename('fullpath')));  % Go up from private/ to interfaces/
prima_matlab_dir = fullfile(mfiledir, 'prima_matlab');

old_path = path;
addpath(genpath(prima_matlab_dir));
try
    solver = lower(solver_name);
    if strcmp(solver, 'uobyqa')
        uo = uobyqa_mod();
        [varargout{1:nargout}] = uo.uobyqa(varargin{:});
    elseif strcmp(solver, 'newuoa')
        no = newuoa_mod();
        [varargout{1:nargout}] = no.newuoa(varargin{:});
    elseif strcmp(solver, 'bobyqa')
        bo = bobyqa_mod();
        [varargout{1:nargout}] = bo.bobyqa(varargin{:});
    elseif strcmp(solver, 'lincoa')
        lo = lincoa_mod();
        [varargout{1:nargout}] = lo.lincoa(varargin{:});
    elseif strcmp(solver, 'cobyla')
        co = cobyla_mod();
        [varargout{1:nargout}] = co.cobyla(varargin{:});
    else
        error('%s: UNEXPECTED ERROR: unknown solver ''%s''.', mfilename(), solver_name);
    end
    path(old_path);
catch exception
    path(old_path);
    rethrow(exception);
end
return
