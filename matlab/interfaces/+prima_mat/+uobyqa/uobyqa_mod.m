classdef uobyqa_mod
    %--------------------------------------------------------------------------------------------------%
    % UOBYQA_MOD is a module providing the reference implementation of Powell's UOBYQA algorithm in
    %
    % M. J. D. Powell, UOBYQA: unconstrained optimization by quadratic approximation, Math. Program.,
    % 92(B):555--582, 2002
    %
    % UOBYQA approximately solves
    %
    %   min F(X),
    %
    % where X is a vector of variables that has N components and F is a real-valued objective function.
    % It tackles the problem by a trust region method that forms quadratic models by interpolation.
    %
    % Coded by Zaikun ZHANG (www.zhangzk.net) based on the UOBYQA paper and Powell's code, with
    % modernization, bug fixes, and improvements.
    %
    % Dedicated to the late Professor M. J. D. Powell FRS (1936--2015).
    %
    % Started: February 2022
    %
    % Last Modified: Wed 10 Sep 2025 02:03:43 AM CST
    %--------------------------------------------------------------------------------------------------%

    methods
        function [x, f_loc, nf_loc, xhist, fhist, info_loc] = uobyqa(~, calfun, x, varargin)
            %--------------------------------------------------------------------------------------------------%
            % Among all the arguments, only CALFUN and X are obligatory. The others are OPTIONAL and you can
            % neglect them unless you are familiar with the algorithm. Any unspecified optional input will take
            % the default value detailed below. For instance, we may invoke the solver as follows.
            %
            % % First define CALFUN and X, and then do the following.
            % call uobyqa(calfun, x, f)
            %
            % or
            %
            % % First define CALFUN and X, and then do the following.
            % call uobyqa(calfun, x, f, rhobeg = 1.0D0, rhoend = 1.0D-6)
            %
            % See examples/uobyqa_exmp.f90 for a concrete example.
            %
            % A detailed introduction to the arguments is as follows.
            % N.B.: RP and IK are defined in the module CONSTS_MOD. See consts.F90 under the directory name
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
            %      named UOBYQA_output.txt; the file will be created if it does not exist; the new output will
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
            %   If XHIST is present, its size at exit will be (N, min(NF, MAXHIST)); if FHIST is present, its
            %   size at exit will be min(NF, MAXHIST).
            %
            %   IMPORTANT NOTICE:
            %   Setting MAXHIST to a large value can be costly in terms of memory for large problems.
            %   MAXHIST will be reset to a smaller value if the memory needed exceeds MAXHISTMEM defined in
            %   CONSTS_MOD (see consts.F90 under the directory named "common").
            %   Use *HIST with caution!!! (N.B.: the algorithm is NOT designed for large problems).
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
            %   NAN_INF_X: NaN or Inf occurs in X.
            %   %--------------------------------------------------------------------------%
            %   The following case(s) should NEVER occur unless there is a bug.
            %   NAN_INF_F: the objective function returns NaN or +Inf;
            %   TRSUBP_FAILED: a trust region step has failed to reduce the model;
            %   %--------------------------------------------------------------------------%
            %--------------------------------------------------------------------------------------------------%


            evaluate_obj = prima_mat.common.evaluate_mod();
            history_obj = prima_mat.common.history_mod();

            preproc_obj = prima_mat.common.preproc_mod();

            % Solver-specific modules
            uobyqb_obj = prima_mat.uobyqa.uobyqb_mod();

            % X(N)


            % FHIST(MAXFHIST)
            % XHIST(N, MAXXHIST)


            solver = "UOBYQA";

            info_loc = NaN;

            nf_loc = NaN;

            eta1_loc = NaN;

            f_loc = NaN;

            % FHIST_LOC(MAXFHIST)
            % XHIST_LOC(N, MAXXHIST)


            n = numel(x);
            npt = (n + 1) * (n + 2) / 2;
            if ~(npt > 0)
                error("NPT > 0");
            end % Validate that NPT does not overflow.

            % Replace any NaN in X by ZERO and Inf/-Inf in X by REALMAX/-REALMAX.
            x(:) = evaluate_obj.moderatex(x);

            % Read the inputs.

            % If RHOBEG is present, then RHOBEG_LOC is a copy of RHOBEG; otherwise, RHOBEG_LOC takes the default
            % value for RHOBEG, taking the value of RHOEND into account. Note that RHOEND is considered only if
            % it is present and it is VALID (i.e., finite and positive). The other inputs are read similarly.
            ipObj = inputParser();
            addParameter(ipObj, 'f', NaN);
            addParameter(ipObj, 'nf', NaN);
            addParameter(ipObj, 'rhobeg', NaN);
            addParameter(ipObj, 'rhoend', NaN);
            addParameter(ipObj, 'ftarget', -realmax);
            addParameter(ipObj, 'maxfun', NaN);
            addParameter(ipObj, 'iprint', 0);
            addParameter(ipObj, 'eta1', NaN);
            addParameter(ipObj, 'eta2', NaN);
            addParameter(ipObj, 'gamma1', 0.5);
            addParameter(ipObj, 'gamma2', 2.0);
            addParameter(ipObj, 'xhist', NaN);
            addParameter(ipObj, 'fhist', NaN);
            addParameter(ipObj, 'maxhist', NaN);
            addParameter(ipObj, 'callback_fcn', struct());
            addParameter(ipObj, 'info', NaN);
            parse(ipObj, varargin{:});

            rhobeg = ipObj.Results.rhobeg;
            rhoend = ipObj.Results.rhoend;
            ftarget_loc = ipObj.Results.ftarget;
            maxfun = ipObj.Results.maxfun;
            iprint_loc = ipObj.Results.iprint;
            eta1 = ipObj.Results.eta1;
            eta2 = ipObj.Results.eta2;
            gamma1_loc = ipObj.Results.gamma1;
            gamma2_loc = ipObj.Results.gamma2;
            xhist = ipObj.Results.xhist;
            fhist = ipObj.Results.fhist;
            maxhist = ipObj.Results.maxhist;
            callback_fcn = ipObj.Results.callback_fcn;

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
                maxfun_loc = max(500 * n, npt + 1);
            else
                maxfun_loc = maxfun;
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
                maxhist_loc = max([maxfun_loc, npt + 1, 500 * n], [], 'all');
            else
                maxhist_loc = maxhist;
            end

            % Preprocess the inputs in case some of them are invalid.
            [iprint_loc, maxfun_loc, maxhist_loc, ftarget_loc, rhobeg_loc, rhoend_loc, ~, ~, ~, ~, eta1_loc, eta2_loc, gamma1_loc, gamma2_loc] = preproc_obj.preproc(solver, n, iprint_loc, maxfun_loc, maxhist_loc, ftarget_loc, rhobeg_loc, rhoend_loc, 'eta1', eta1_loc, 'eta2', eta2_loc, 'gamma1', gamma1_loc, 'gamma2', gamma2_loc);

            % Further revise MAXHIST_LOC according to MAXHISTMEM, and allocate memory for the history.
            % In MATLAB/Python/Julia/R implementation, we should simply set MAXHIST = MAXFUN and initialize
            % FHIST = NaN(1, MAXFUN), XHIST = NaN(N, MAXFUN) if they are requested; replace MAXFUN with 0 for
            % the history that is not requested.
            [maxhist_loc, xhist_loc, fhist_loc] = history_obj.prehist(maxhist_loc, n, nargout >= 4, nargout >= 5);

            %-------------------- Call UOBYQB, which performs the real calculations. --------------------------%
            if ismember('callback_fcn', ipObj.UsingDefaults)
                [x, nf_loc, f_loc, fhist_loc, xhist_loc, info_loc] = uobyqb_obj.uobyqb(calfun, iprint_loc, maxfun_loc, eta1_loc, eta2_loc, ftarget_loc, gamma1_loc, gamma2_loc, rhobeg_loc, rhoend_loc, x, fhist_loc, xhist_loc);
            else
                [x, nf_loc, f_loc, fhist_loc, xhist_loc, info_loc] = uobyqb_obj.uobyqb(calfun, iprint_loc, maxfun_loc, eta1_loc, eta2_loc, ftarget_loc, gamma1_loc, gamma2_loc, rhobeg_loc, rhoend_loc, x, fhist_loc, xhist_loc, 'callback_fcn', callback_fcn);
            end
            %--------------------------------------------------------------------------------------------------%


            % Write the outputs.


            % Copy XHIST_LOC to XHIST if needed.
            if nargout >= 4
                nhist = min(nf_loc, size(xhist_loc, 2));
                %----------------------------------------------------%
                xhist = NaN(n, nhist); % Removable in F2003.
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
                fhist = NaN(nhist, 1); % Removable in F2003.
                %--------------------------------------------------%
                fhist = fhist_loc(1:nhist); % The same as XHIST, we must cap FHIST at NF_LOC.

            end

            % If MAXFHIST_IN >= NF_LOC > MAXFHIST_LOC, warn that not all history is recorded.
            if (nargout >= 4 || nargout >= 5) && maxhist_loc < nf_loc
            end

        end

    end
end