classdef bobyqa_mod
    %--------------------------------------------------------------------------------------------------%
    % BOBYQA_MOD is a module providing the reference implementation of Powell's BOBYQA algorithm in
    %
    % M. J. D. Powell, The BOBYQA algorithm for bound constrained optimization without derivatives,
    % Technical Report DAMTP 2009/NA06, Department of Applied Mathematics and Theoretical Physics,
    % Cambridge University, Cambridge, UK, 2009
    %
    % BOBYQA approximately solves
    %
    %   min F(X) subject to XL <= X <= XU,
    %
    % where X is a vector of variables that has N components, and F is a real-valued objective function.
    % XL and XU are a pair of N-dimensional vectors indicating the lower and upper bounds of X. The
    % algorithm assumes that XL < XU entrywise. It tackles the problem by applying a trust region method
    % that forms quadratic models by interpolation. There is usually some freedom in the interpolation
    % conditions, which is taken up by minimizing the Frobenius norm of the change to the second
    % derivative of the model, beginning with the ZERO matrix. The values of the variables are
    % constrained by upper and lower bounds. The arguments of the subroutine are as follows.
    %
    % Coded by Zaikun ZHANG (www.zhangzk.net) based on the BOBYQA paper and Powell's code, with
    % modernization, bug fixes, and improvements.
    %
    % Dedicated to the late Professor M. J. D. Powell FRS (1936--2015).
    %
    % Started: February 2022
    %
    % Last Modified: Thursday, February 22, 2024 PM03:30:31
    %--------------------------------------------------------------------------------------------------%

    methods
        function [x, f_loc, nf_loc, xhist, fhist, info] = bobyqa(~, calfun, x, varargin)
            %--------------------------------------------------------------------------------------------------%
            % Among all the arguments, only CALFUN and X are obligatory. The others are OPTIONAL and you can
            % neglect them unless you are familiar with the algorithm. Any unspecified optional input will take
            % the default value detailed below. For instance, we may invoke the solver as follows.
            %
            % % First define CALFUN and X, and then do the following.
            % call bobyqa(calfun, x, f)
            %
            % or
            %
            % % First define CALFUN, X, and XL, and then do the following.
            % call bobyqa(calfun, x, f, xl = xl, rhobeg = 1.0D0, rhoend = 1.0D-6)
            %
            % See examples/bobyqa_exmp.f90 for a concrete example.
            %
            % A detailed introduction to the arguments is as follows.
            % N.B.: RP and IK are defined in the module CONSTS_MOD. See consts.F90 under the directory named
            % "common". By default, RP = kind(0.0D0) and IK = kind(0), with REAL(RP) being the double-precision
            % real, and INTEGER(IK) being the default integer. For ADVANCED USERS, RP and IK can be defined by
            % setting PRIMA_REAL_PRECISION and PRIMA_INTEGER_KIND in common/ppf.h. Use the default if unsure.
            %
            % CALFUN
            %   Input, subroutine.
            %   CALFUN(X, F) should evaluate the objective function at the given REAL(RP) vector X and set the
            %   value to the REAL(RP) scalar F. It must be provided by the user, and its definition must conform
            %   to the following interface:
            %   %-------------------------------------------------------------------------%
            %    subroutine calfun(x, f)
            %    real(RP), intent(in) :: x(:)
            %    real(RP), intent(out) :: f
            %    end subroutine calfun
            %   %-------------------------------------------------------------------------%
            %
            % X
            %   Input and output, REAL(RP) vector.
            %   As an input, X should be an N dimensional vector that contains the starting point, N being the
            %   dimension of the problem. As an output, X will be set to an approximate minimizer.
            %
            % F
            %   Output, REAL(RP) scalar.
            %   F will be set to the objective function value of X at exit.
            %
            % XL, XU
            %   Input, REAL(RP) vectors, default: XL = [], XU = [].
            %   XL is the lower bound for X. Its size is either N or 0, the latter signifying that X has no
            %   lower bound. Any entry of XL that is NaN or below -BOUNDMAX will be taken as -BOUNDMAX, which
            %   effectively means there is no lower bound for the corresponding entry of X. The value of
            %   BOUNDMAX is 0.25*HUGE(X), which is about 8.6E37 for single precision and 4.5E307 for double
            %   precision. XU is similar.
            %   N.B.:
            %   1. It is required that XU - XL > 2*EPSILON(X), which is about 2.4E-7 for single precision and
            %   4.5E-16 for double precision. Otherwise, the solver will return after printing a warning.
            %   2. Why don't we set BOUNDMAX to REALMAX? Because we want to avoid overflow when calculating
            %   XU - XL and when defining/updating SU and SL. This is not a problem in MATLAB/Python/Julia/R.
            %
            % NF
            %   Output, INTEGER(IK) scalar.
            %   NF will be set to the number of calls of CALFUN at exit.
            %
            % RHOBEG, RHOEND
            %   Inputs, REAL(RP) scalars, default: RHOBEG = 1, RHOEND = 10^-6. RHOBEG and RHOEND must be set to
            %   the initial and final values of a trust-region radius, both being positive and RHOEND <= RHOBEG.
            %   Typically RHOBEG should be about one tenth of the greatest expected change to a variable, and
            %   RHOEND should indicate the accuracy that is required in the final values of the variables.
            %
            % FTARGET
            %   Input, REAL(RP) scalar, default: -Inf.
            %   FTARGET is the target function value. The algorithm will terminate when a point with a function
            %   value <= FTARGET is found.
            %
            % MAXFUN
            %   Input, INTEGER(IK) scalar, default: MAXFUN_DIM_DFT*N with MAXFUN_DIM_DFT defined in the module
            %   CONSTS_MOD (see common/consts.F90). MAXFUN is the maximal number of calls of CALFUN.
            %
            % NPT
            %   Input, INTEGER(IK) scalar, default: 2N + 1.
            %   NPT is the number of interpolation conditions for each trust region model. Its value must be in
            %   the interval [N+2, (N+1)(N+2)/2]. Powell commented that "the value NPT = 2*N+1 being recommended
            %   for a start ... much larger values tend to be inefficient, because the amount of routine work of
            %   each iteration is of magnitude NPT**2, and because the achievement of adequate accuracy in some
            %   matrix calculations becomes more difficult. Some excellent numerical results have been found in
            %   the case NPT=N+6 even with more than 100 variables." And "choices that exceed 2*N+1 are not
            %   recommended" by Powell.
            %
            % IPRINT
            %   Input, INTEGER(IK) scalar, default: 0.
            %   The value of IPRINT should be set to 0, 1, -1, 2, -2, 3, or -3, which controls how much
            %   information will be printed during the computation:
            %   0: there will be no printing;
            %   1: a message will be printed to the screen at the return, showing the best vector of variables
            %      found and its objective function value;
            %   2: in addition to 1, each new value of RHO is printed to the screen, with the best vector of
            %      variables so far and its objective function value;
            %   3: in addition to 2, each function evaluation with its variables will be printed to the screen;
            %   -1, -2, -3: the same information as 1, 2, 3 will be printed, not to the screen but to a file
            %      named BOBYQA_output.txt; the file will be created if it does not exist; the new output will
            %      be appended to the end of this file if it already exists.
            %   Note that IPRINT = +/-3 can be costly in terms of time and/or space.
            %
            % ETA1, ETA2, GAMMA1, GAMMA2
            %   Input, REAL(RP) scalars, default: ETA1 = 0.1, ETA2 = 0.7, GAMMA1 = 0.5, and GAMMA2 = 2.
            %   ETA1, ETA2, GAMMA1, and GAMMA2 are parameters in the updating scheme of the trust-region radius
            %   detailed in the subroutine TRRAD in trustregion.f90. Roughly speaking, the trust-region radius
            %   is contracted by a factor of GAMMA1 when the reduction ratio is below ETA1, and enlarged by a
            %   factor of GAMMA2 when the reduction ratio is above ETA2. It is required that 0 < ETA1 <= ETA2
            %   < 1 and 0 < GAMMA1 < 1 < GAMMA2. Normally, ETA1 <= 0.25. It is NOT advised to set ETA1 >= 0.5.
            %
            % XHIST, FHIST, MAXHIST
            %   XHIST: Output, ALLOCATABLE rank 2 REAL(RP) array;
            %   FHIST: Output, ALLOCATABLE rank 1 REAL(RP) array;
            %   MAXHIST: Input, INTEGER(IK) scalar, default: MAXFUN
            %   XHIST, if present, will output the history of iterates, while FHIST, if present, will output the
            %   history function values. MAXHIST should be a nonnegative integer, and XHIST/FHIST will output
            %   only the history of the last MAXHIST iterations. Therefore, MAXHIST = 0 means XHIST/FHIST will
            %   output nothing, while setting MAXHIST = MAXFUN requests XHIST/FHIST to output all the history.
            %   If XHIST is present, its size at exit will be [N, min(NF, MAXHIST)]; if FHIST is present, its
            %   size at exit will be min(NF, MAXHIST).
            %
            %   IMPORTANT NOTICE:
            %   Setting MAXHIST to a large value can be costly in terms of memory for large problems.
            %   MAXHIST will be reset to a smaller value if the memory needed exceeds MAXHISTMEM defined in
            %   CONSTS_MOD (see consts.F90 under the directory named "common").
            %   Use *HIST with caution! (N.B.: the algorithm is NOT designed for large problems).
            %
            % HONOUR_X0
            %  Input, LOGICAL scalar, default: it is .false. if RHOBEG is present and 0 < RHOBEG < Inf, and it
            %  is .true. otherwise. HONOUR_X0 indicates whether to respect the user-defined X0 or not.
            %  BOBYQA requires that the distance between X0 and the inactive bounds is at least RHOBEG. X0 or
            %  RHOBEG is revised if this requirement is not met. If HONOUR_X0 == TRUE, revise RHOBEG if needed;
            %  otherwise, revise X0 if needed. See the PREPROC subroutine for more information.
            %
            % CALLBACK_FCN
            %   Input, function to report progress and optionally request termination.
            %
            % INFO
            %   Output, INTEGER(IK) scalar.
            %   INFO is the exit flag. It will be set to one of the following values defined in the module
            %   INFOS_MOD (see common/infos.f90):
            %   SMALL_TR_RADIUS: the lower bound for the trust region radius is reached;
            %   FTARGET_ACHIEVED: the target function value is reached;
            %   MAXFUN_REACHED: the objective function has been evaluated MAXFUN times;
            %   MAXTR_REACHED: the trust region iteration has been performed MAXTR times (MAXTR = 2*MAXFUN);
            %   NAN_INF_MODEL: NaN or Inf occurs in the model;
            %   NAN_INF_X: NaN or Inf occurs in X;
            %   DAMAGING_ROUNDING: the rounding error becomes damaging;
            %   NO_SPACE_BETWEEN_BOUNDS: there is not enough space between some lower and upper bounds, namely
            %   one of the difference XU(I)-XL(I) is less than 2*RHOBEG.
            %   %--------------------------------------------------------------------------%
            %   The following case(s) should NEVER occur unless there is a bug.
            %   NAN_INF_F: the objective function returns NaN or +Inf;
            %   TRSUBP_FAILED: a trust region step failed to reduce the model.
            %   %--------------------------------------------------------------------------%
            %--------------------------------------------------------------------------------------------------%


            evaluate_obj = prima_mat.common.evaluate_mod();
            history_obj = prima_mat.common.history_mod();

            preproc_obj = prima_mat.common.preproc_mod();

            % Solver-specific modules
            bobyqb_obj = prima_mat.bobyqa.bobyqb_mod();

            % X(N)


            % XL(N)
            % XU(N)


            % FHIST(MAXFHIST)
            % XHIST(N, MAXXHIST)


            solver = "BOBYQA";

            nf_loc = NaN;

            eta1_loc = NaN;

            f_loc = NaN;

            xl_loc = NaN(size(x));
            xu_loc = NaN(size(x));
            % FHIST_LOC(MAXFHIST)
            % XHIST_LOC(N, MAXXHIST)


            n = numel(x);

            ipObj = inputParser();
            addParameter(ipObj, 'f', NaN);
            addParameter(ipObj, 'xl', NaN);
            addParameter(ipObj, 'xu', NaN);
            addParameter(ipObj, 'nf', NaN);
            addParameter(ipObj, 'rhobeg', NaN);
            addParameter(ipObj, 'rhoend', NaN);
            addParameter(ipObj, 'ftarget', -realmax);
            addParameter(ipObj, 'maxfun', NaN);
            addParameter(ipObj, 'npt', NaN);
            addParameter(ipObj, 'iprint', 0);
            addParameter(ipObj, 'eta1', NaN);
            addParameter(ipObj, 'eta2', NaN);
            addParameter(ipObj, 'gamma1', 0.5);
            addParameter(ipObj, 'gamma2', 2.0);
            addParameter(ipObj, 'xhist', NaN);
            addParameter(ipObj, 'fhist', NaN);
            addParameter(ipObj, 'maxhist', NaN);
            addParameter(ipObj, 'honour_x0', false);
            addParameter(ipObj, 'callback_fcn', struct());
            addParameter(ipObj, 'info', NaN);
            parse(ipObj, varargin{:});

            xl = ipObj.Results.xl;
            xu = ipObj.Results.xu;

            rhobeg = ipObj.Results.rhobeg;
            rhoend = ipObj.Results.rhoend;
            ftarget_loc = ipObj.Results.ftarget;
            maxfun = ipObj.Results.maxfun;
            npt = ipObj.Results.npt;
            iprint_loc = ipObj.Results.iprint;
            eta1 = ipObj.Results.eta1;
            eta2 = ipObj.Results.eta2;
            gamma1_loc = ipObj.Results.gamma1;
            gamma2_loc = ipObj.Results.gamma2;
            xhist = ipObj.Results.xhist;
            fhist = ipObj.Results.fhist;
            maxhist = ipObj.Results.maxhist;
            honour_x0 = ipObj.Results.honour_x0;
            callback_fcn = ipObj.Results.callback_fcn;
            info = ipObj.Results.info;

            % Read the inputs

            xl_loc(:) = -(0.25 * realmax);
            if ~ismember('xl', ipObj.UsingDefaults)
                if numel(xl) > 0
                    xl_loc = xl;
                end
            end
            xl_loc(isnan(xl_loc) | xl_loc < -(0.25 * realmax)) = -(0.25 * realmax);

            xu_loc(:) = 0.25 * realmax;
            if ~ismember('xu', ipObj.UsingDefaults)
                if numel(xu) > 0
                    xu_loc = xu;
                end
            end
            xu_loc(isnan(xu_loc) | xu_loc > 0.25 * realmax) = 0.25 * realmax;

            % The solver requires that MINVAL(XU-XL) >= 2*RHOBEG, and we return if MINVAL(XU-XL) < 2*EPS.
            % It would be better to fix the variables at (XU+XL)/2 wherever XU and XL almost equal, as is done
            % in the MATLAB/Python interface of the solvers. In Fortran, this is doable using internal functions,
            % but we choose not to implement it in the current version.
            if any(xu_loc - xl_loc < 2.0 * eps(1.0), 'all')
                if nargout >= 6
                    info = 6;
                end

                return
            end

            x = max(xl_loc, min(xu_loc, evaluate_obj.moderatex(x)));

            % If RHOBEG is present, then RHOBEG_LOC is a copy of RHOBEG; otherwise, RHOBEG_LOC takes the default
            % value for RHOBEG, taking the value of RHOEND into account. Note that RHOEND is considered only if
            % it is present and it is VALID (i.e., finite and positive). The other inputs are read similarly.
            if ~ismember('rhobeg', ipObj.UsingDefaults)
                rhobeg_loc = rhobeg;
            elseif ~ismember('rhoend', ipObj.UsingDefaults)
                % Fortran does not take short-circuit evaluation of logic expressions. Thus it is WRONG to
                % combine the evaluation of PRESENT(RHOEND) and the evaluation of IS_FINITE(RHOEND) as
                % "IF (PRESENT(RHOEND) .AND. IS_FINITE(RHOEND))". The compiler may choose to evaluate the
                % IS_FINITE(RHOEND) even if PRESENT(RHOEND) is false!
                if isfinite(rhoend) && rhoend > 0
                    rhobeg_loc = max(10.0 * rhoend, 1.0);
                else
                    rhobeg_loc = 1.0;
                end
            else
                rhobeg_loc = 1.0;
            end

            if ~ismember('rhoend', ipObj.UsingDefaults)
                rhoend_loc = rhoend;
            elseif rhobeg_loc > 0
                rhoend_loc = max(eps(1.0), min((1.0e-6 / 1.0) * rhobeg_loc, 1.0e-6));
            else
                rhoend_loc = 1.0e-6;
            end

            if ismember('maxfun', ipObj.UsingDefaults)
                maxfun_loc = 500 * n;
            else
                maxfun_loc = maxfun;
            end

            if ~ismember('npt', ipObj.UsingDefaults)
                npt_loc = npt;
            elseif maxfun_loc >= n + 3
                % Take MAXFUN into account if it is valid.
                npt_loc = min(maxfun_loc - 1, 2 * n + 1);
            else
                npt_loc = 2 * n + 1;
            end

            if ~ismember('eta1', ipObj.UsingDefaults)
                eta1_loc = eta1;
            elseif ~ismember('eta2', ipObj.UsingDefaults)
                if eta2 > 0 && eta2 < 1
                    eta1_loc = max(eps(1.0), eta2 / 7.0);
                end
            else
                eta1_loc = 0.1;
            end

            if ~ismember('eta2', ipObj.UsingDefaults)
                eta2_loc = eta2;
            elseif eta1_loc > 0 && eta1_loc < 1
                eta2_loc = (eta1_loc + 2.0) / 3.0;
            else
                eta2_loc = 0.7;
            end

            if ismember('maxhist', ipObj.UsingDefaults)
                maxhist_loc = max([maxfun_loc, n + 3, 500 * n], [], 'all');
            else
                maxhist_loc = maxhist;
            end

            has_rhobeg = ~ismember('rhobeg', ipObj.UsingDefaults);
            honour_x0_loc = true;
            if ~ismember('honour_x0', ipObj.UsingDefaults)
                honour_x0_loc = honour_x0;
            elseif has_rhobeg
                % HONOUR_X0 is FALSE if user provides a valid RHOBEG. Is this the best choice?
                honour_x0_loc = (~(isfinite(rhobeg) && rhobeg > 0));
            end

            % Preprocess the inputs in case some of them are invalid. It does nothing if all inputs are valid.
            [iprint_loc, maxfun_loc, maxhist_loc, ftarget_loc, rhobeg_loc, rhoend_loc, npt_loc, ~, ~, ~, eta1_loc, eta2_loc, gamma1_loc, gamma2_loc, x] = preproc_obj.preproc(solver, n, iprint_loc, maxfun_loc, maxhist_loc, ftarget_loc, rhobeg_loc, rhoend_loc, 'npt', npt_loc, 'eta1', eta1_loc, 'eta2', eta2_loc, 'gamma1', gamma1_loc, 'gamma2', gamma2_loc, 'has_rhobeg', has_rhobeg, 'honour_x0', honour_x0_loc, 'xl', xl_loc, 'xu', xu_loc, 'x0', x);

            % Further revise MAXHIST_LOC according to MAXHISTMEM, and allocate memory for the history.
            % In MATLAB/Python/Julia/R implementation, we should simply set MAXHIST = MAXFUN and initialize
            % FHIST = NaN(1, MAXFUN), XHIST = NaN(N, MAXFUN)
            % if they are requested; replace MAXFUN with 0 for the history that is not requested.
            [maxhist_loc, xhist_loc, fhist_loc] = history_obj.prehist(maxhist_loc, n, nargout >= 4, nargout >= 5);

            %-------------------- Call BOBYQB, which performs the real calculations. --------------------------%
            if ismember('callback_fcn', ipObj.UsingDefaults)
                [x, nf_loc, f_loc, fhist_loc, xhist_loc, info_loc] = bobyqb_obj.bobyqb(calfun, iprint_loc, maxfun_loc, npt_loc, eta1_loc, eta2_loc, ftarget_loc, gamma1_loc, gamma2_loc, rhobeg_loc, rhoend_loc, xl_loc, xu_loc, x, fhist_loc, xhist_loc);
            else
                [x, nf_loc, f_loc, fhist_loc, xhist_loc, info_loc] = bobyqb_obj.bobyqb(calfun, iprint_loc, maxfun_loc, npt_loc, eta1_loc, eta2_loc, ftarget_loc, gamma1_loc, gamma2_loc, rhobeg_loc, rhoend_loc, xl_loc, xu_loc, x, fhist_loc, xhist_loc, 'callback_fcn', callback_fcn);
            end
            %--------------------------------------------------------------------------------------------------%

            % Write the outputs.


            if nargout >= 6
                info = info_loc;
            end

            % Copy XHIST_LOC to XHIST if needed.
            if nargout >= 4
                nhist = min(nf_loc, size(xhist_loc, 2));
                %----------------------------------------------------%
                % Removable in F2003.
                %----------------------------------------------------%
                xhist = xhist_loc(:, 1:nhist);
                % N.B.:
                % 0. Allocate XHIST as long as it is present, even if the size is 0; otherwise, it will be
                % illegal to enquire XHIST after exit.
                % 1. Even though Fortran 2003 supports automatic (re)allocation of allocatable arrays upon
                % intrinsic assignment, we keep the line of SAFEALLOC, because some very new compilers (Absoft
                % Fortran 21.0) are still not standard-compliant in this respect.
                % 2. NF may not be present. Hence we should NOT use NF but NF_LOC.
                % 3. When SIZE(XHIST_LOC, 2) > NF_LOC, which is the normal case in practice, XHIST_LOC contains
                % GARBAGE in XHIST_LOC(:, NF_LOC + 1 : END). Therefore, we MUST cap XHIST at NF_LOC so that
                % XHIST contains only valid history. For this reason, there is no way to avoid allocating
                % two copies of memory for XHIST unless we declare it to be a POINTER instead of ALLOCATABLE.

            end
            % F2003 automatically deallocate local ALLOCATABLE variables at exit, yet we prefer to deallocate
            % them immediately when they finish their jobs.


            % Copy FHIST_LOC to FHIST if needed.
            if nargout >= 5
                nhist = min(nf_loc, numel(fhist_loc));
                %--------------------------------------------------%
                % Removable in F2003.
                %--------------------------------------------------%
                fhist = fhist_loc(1:nhist); % The same as XHIST, we must cap FHIST at NF_LOC.

            end

            % If NF_LOC > MAXHIST_LOC, warn that not all history is recorded.
            if (nargout >= 4 || nargout >= 5) && maxhist_loc < nf_loc
            end

        end

    end
end