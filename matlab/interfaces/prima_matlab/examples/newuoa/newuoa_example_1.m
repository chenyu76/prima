%---------------------------------------- THE MAIN PROGRAM ----------------------------------------%


% The following line makes the solver available.
newuoa_obj = newuoa_mod();

% The following line specifies which module provides CALFUN.
calfun_obj = calfun_mod();


n = 2;
nf = NaN; info = NaN;
f = NaN; x = NaN(n, 1); x0 = NaN(n, 1);

% Define the starting point.
x0 = zeros(size(x0));

% The following lines illustrates how to call the solver.
x(:) = x0;
[x, f] = newuoa_obj.newuoa(@calfun_obj.calfun, x); % This call will not print anything.

% In addition to the compulsory arguments, the following illustration specifies also RHOBEG and
% IPRINT, which are optional. All the unspecified optional arguments (RHOEND, MAXFUN, etc.) will
% take their default values coded in the solver.
x(:) = x0;
[x, f, nf, ~, ~, info] = newuoa_obj.newuoa(@calfun_obj.calfun, x, 'rhobeg', 1.0, 'iprint', 1, 'callback_fcn', @calfun_obj.callback_fcn);

