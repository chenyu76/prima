classdef lincoa_mod
    %--------------------------------------------------------------------------------------------------%
    % LINCOA_MOD is a module providing the reference implementation of Powell's LINCOA algorithm.
    %
    % The algorithm approximately solves
    %
    %   min F(X) subject to Aineq*X <= Bineq, Aeq*x = Beq, XL <= X <= XU,
    %
    % where X is a vector of variables that has N components, F is a real-valued objective function,
    % Aineq is an Mineq-by-N matrix, Bineq is an Mineq-dimensional real vector, Aeq is an Meq-by-N
    % matrix, Beq is an Meq-dimensional real vector, XL is an N-dimensional real vector, and XU is
    % an N-dimensional real vector.
    %
    % It tackles the problem by a trust region method that forms quadratic models by interpolation.
    % Usually there is much freedom in each new model after satisfying the interpolation conditions,
    % which is taken up by minimizing the Frobenius norm of the change to the second derivative matrix
    % of the model. One new function value is calculated on each iteration, usually at a point where
    % the current model predicts a reduction in the least value so far of the objective function subject
    % to the linear constraints. Alternatively, a new vector of variables may be chosen to replace an
    % interpolation point that may be too far away for reliability, and the new point does not have to
    % satisfy the constraints.
    %
    % Coded by Zaikun ZHANG (www.zhangzk.net) based on the paper
    %
    % M. J. D. Powell, On fast trust region methods for quadratic models with linear constraints,
    % Math. Program. Comput., 7:237--267, 2015
    %
    % and Powell's code, with modernization, bug fixes, and improvements.
    %
    % N.B.:
    % 1. Powell did not publish a paper to introduce the algorithm. The above paper does not describe
    % LINCOA but discusses how to solve linearly-constrained trust-region subproblems.
    % 2. Powell's code does not accept linear equality constraints or bound constraints.
    %
    % Dedicated to the late Professor M. J. D. Powell FRS (1936--2015).
    %
    % Started: February 2022
    %
    % Last Modified: Sunday, April 07, 2024 PM04:11:56
    %--------------------------------------------------------------------------------------------------%

    methods
        function [x, f_loc, cstrv_loc, nf_loc, xhist, fhist, chist, info_loc] = lincoa(obj, calfun, x, varargin)
            %--------------------------------------------------------------------------------------------------%
            % Among all the arguments, only CALFUN, and X are obligatory. The others are OPTIONAL and you can
            % neglect them unless you are familiar with the algorithm. Any unspecified optional input will take
            % the default value detailed below. For instance, we may invoke the solver as follows.
            %
            % % First define CALFUN and X, and then do the following.
            % call lincoa(calfun, x, f)
            %
            % or
            %
            % % First define CALFUN, X, Aineq, and Bineq, and then do the following.
            % call lincoa(calfun, x, f, cstrv, Aineq = Aineq, bineq = bineq, rhobeg = 1.0D0, rhoend = 1.0D-6)
            %
            % See examples/lincoa_exmp.f90 for a concrete example.
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
            % CSTRV
            %   Output, REAL(RP) scalar.
            %   CSTRV will be set to the L-infinity constraint violation of X at exit, namely
            %   MAXVAL([0, Aineq*X - Bineq, abs(Aeq*X - Beq), XL - X, X - XU])
            %   N.B.: We use the original constraints to evaluate CSTRV, even though they may be modified during
            %   the computation.
            %
            % Aineq, Bineq
            %   Input, REAL(RP) matrix of size [Mineq, N] and REAL vector of size Mineq unless they are both
            %   empty, default: [] and [].
            %   Aineq and Bineq represent the linear inequality constraints: Aineq*X <= Bineq.
            %
            % Aeq, Beq
            %   Input, REAL(RP) matrix of size [Meq, N] and REAL vector of size Meq unless they are both
            %   empty, default: [] and [].
            %   Aeq and Beq represent the linear equality constraints: Aeq*X = Beq.
            %
            % XL, XU
            %   Input, REAL(RP) vectors of size N unless they are both empty, default: [] and [].
            %   XL is the lower bound for X. Its size is either N or 0, the latter signifying that X has no
            %   lower bound. Any entry of XL that is NaN or below -BOUNDMAX will be taken as -BOUNDMAX, which
            %   effectively means there is no lower bound for the corresponding entry of X. The value of
            %   BOUNDMAX is 0.25*HUGE(X), which is about 8.6E37 for single precision and 4.5E307 for double
            %   precision. XU is similar.
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
            % CTOL
            %   Input, REAL(RP) scalar, default: machine epsilon.
            %   CTOL is the tolerance of constraint violation. X is considered feasible if CSTRV(X) <= CTOL.
            %   N.B.: 1. CTOL is absolute, not relative.
            %   2. CTOL is used for choosing the returned X. It does not affect the iterations of the algorithm.
            %
            % CWEIGHT
            %   Input, REAL(RP) scalar, default: CWEIGHT_DFT defined in the module CONSTS_MOD in common/consts.F90.
            %   CWEIGHT is the weight that the constraint violation takes in the selection of the returned X.
            %
            % MAXFUN
            %   Input, INTEGER(IK) scalar, default: MAXFUN_DIM_DFT*N with MAXFUN_DIM_DFT defined in the module
            %   CONSTS_MOD (see common/consts.F90). MAXFUN is the maximal number of calls of CALFUN.
            %
            % NPT
            %   Input, INTEGER(IK) scalar, default: 2N + 1.
            %   NPT is the number of interpolation conditions for each trust region model. Its value must be in
            %   the interval [N+2, (N+1)(N+2)/2]. Typical choices of Powell were NPT=N+6 and NPT=2*N+1. Powell
            %   commented that "larger values tend to be highly inefficient when the number of variables is
            %   substantial, due to the amount of work and extra difficulty of adjusting more points."
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
            %      named LINCOA_output.txt; the file will be created if it does not exist; the new output will
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
            % XHIST, FHIST, CHIST, MAXHIST
            %   XHIST: Output, ALLOCATABLE rank 2 REAL(RP) array;
            %   FHIST: Output, ALLOCATABLE rank 1 REAL(RP) array;
            %   CHIST: Output, ALLOCATABLE rank 1 REAL(RP) array;
            %   MAXHIST: Input, INTEGER(IK) scalar, default: MAXFUN
            %   XHIST, if present, will output the history of iterates, while FHIST/CHIST, if present, will output
            %   the history function values/constraint violations. MAXHIST should be a nonnegative integer, and
            %   XHIST/FHIST/CHIST will output only the history of the last MAXHIST iterations. Therefore,
            %   MAXHIST = 0 means XHIST/FHIST/CHIST will output nothing, while setting MAXHIST = MAXFUN requests
            %   XHIST/FHIST/CHIST to output all the history.
            %   If XHIST is present, its size at exit will be [N, min(NF, MAXHIST)]; if FHIST/CHIST is present,
            %   its size at exit will be min(NF, MAXHIST).
            %
            %   IMPORTANT NOTICE:
            %   Setting MAXHIST to a large value can be costly in terms of memory for large problems.
            %   MAXHIST will be reset to a smaller value if the memory needed exceeds MAXHISTMEM defined in
            %   CONSTS_MOD (see consts.F90 under the directory named "common").
            %   Use *HIST with caution! (N.B.: the algorithm is NOT designed for large problems).
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
            %   ZERO_LINEAR_CONSTRAINT: one of the linear constraints has a zero gradient
            %   %--------------------------------------------------------------------------%
            %   The following case(s) should NEVER occur unless there is a bug.
            %   NAN_INF_F: the objective function returns NaN or +Inf;
            %   TRSUBP_FAILED: a trust region step failed to reduce the model.
            %   %--------------------------------------------------------------------------%
            %--------------------------------------------------------------------------------------------------%

            % Common modules
            consts_obj = prima_mat.common.consts_mod();


            debug_obj = prima_mat.common.debug_mod();
            evaluate_obj = prima_mat.common.evaluate_mod();
            history_obj = prima_mat.common.history_mod();
            infnan_obj = prima_mat.common.infnan_mod();
            linalg_obj = prima_mat.common.linalg_mod();
            memory_obj = prima_mat.common.memory_mod();

            preproc_obj = prima_mat.common.preproc_mod();
            selectx_obj = prima_mat.common.selectx_mod();
            string_obj = prima_mat.common.string_mod();

            % Solver-specific modules
            lincob_obj = prima_mat.lincoa.lincob_mod();


            % Compulsory arguments
            % N.B.: INTENT cannot be specified if a dummy procedure is not a POINTER
            % X(N)

            % Optional inputs



            % Aeq(Meq, N)
            % Aineq(Mineq, N)
            % Beq(Meq)
            % Bineq(Mineq)



            % XL(N)
            % XU(N)

            % Optional outputs



            % CHIST(MAXCHIST)
            % FHIST(MAXFHIST)
            % XHIST(N, MAXXHIST)

            % Local variables
            solver = "LINCOA";
            srname = "LINCOA";
            info_loc = NaN;


            nf_loc = NaN;


            cstrv_loc = NaN;


            eta1_loc = NaN;

            f_loc = NaN;


            xl_loc = NaN(numel(x), 1);
            xu_loc = NaN(numel(x), 1);
            % Aeq_LOC(Meq, N)
            % Aineq_LOC(Mineq, N)
            amat = NaN; % AMAT(N, M); each column corresponds to a constraint
            % Beq_LOC(Meq)
            % Bineq_LOC(Mineq)
            bvec = NaN; % BVEC(M)
            % CHIST_LOC(MAXCHIST)
            % FHIST_LOC(MAXFHIST)
            % XHIST_LOC(N, MAXXHIST)

            % Sizes
            ipObj = inputParser();
            addParameter(ipObj, 'f', NaN);
            addParameter(ipObj, 'cstrv', NaN);
            addParameter(ipObj, 'Aineq', NaN);
            addParameter(ipObj, 'bineq', NaN);
            addParameter(ipObj, 'Aeq', NaN);
            addParameter(ipObj, 'beq', NaN);
            addParameter(ipObj, 'xl', NaN);
            addParameter(ipObj, 'xu', NaN);
            addParameter(ipObj, 'nf', NaN);
            addParameter(ipObj, 'rhobeg', NaN);
            addParameter(ipObj, 'rhoend', NaN);
            addParameter(ipObj, 'ftarget', NaN);
            addParameter(ipObj, 'ctol', NaN);
            addParameter(ipObj, 'cweight', NaN);
            addParameter(ipObj, 'maxfun', NaN);
            addParameter(ipObj, 'npt', NaN);
            addParameter(ipObj, 'iprint', NaN);
            addParameter(ipObj, 'eta1', NaN);
            addParameter(ipObj, 'eta2', NaN);
            addParameter(ipObj, 'gamma1', NaN);
            addParameter(ipObj, 'gamma2', NaN);
            addParameter(ipObj, 'xhist', NaN);
            addParameter(ipObj, 'fhist', NaN);
            addParameter(ipObj, 'chist', NaN);
            addParameter(ipObj, 'maxhist', NaN);
            addParameter(ipObj, 'maxfilt', NaN);
            addParameter(ipObj, 'callback_fcn', struct());
            addParameter(ipObj, 'info', NaN);
            parse(ipObj, varargin{:});


            Aineq = ipObj.Results.Aineq;
            bineq = ipObj.Results.bineq;
            Aeq = ipObj.Results.Aeq;
            beq = ipObj.Results.beq;
            xl = ipObj.Results.xl;
            xu = ipObj.Results.xu;

            rhobeg = ipObj.Results.rhobeg;
            rhoend = ipObj.Results.rhoend;
            ftarget = ipObj.Results.ftarget;
            ctol = ipObj.Results.ctol;
            cweight = ipObj.Results.cweight;
            maxfun = ipObj.Results.maxfun;
            npt = ipObj.Results.npt;
            iprint = ipObj.Results.iprint;
            eta1 = ipObj.Results.eta1;
            eta2 = ipObj.Results.eta2;
            gamma1 = ipObj.Results.gamma1;
            gamma2 = ipObj.Results.gamma2;
            xhist = ipObj.Results.xhist;
            fhist = ipObj.Results.fhist;
            chist = ipObj.Results.chist;
            maxhist = ipObj.Results.maxhist;
            maxfilt = ipObj.Results.maxfilt;
            callback_fcn = ipObj.Results.callback_fcn;

            if ismember('bineq', ipObj.UsingDefaults)
                mineq = 0;
            else
                mineq = numel(bineq);
            end
            if ismember('beq', ipObj.UsingDefaults)
                meq = 0;
            else
                meq = numel(beq);
            end
            n = numel(x);

            % Preconditions
            if consts_obj.DEBUGGING
                debug_obj.assert(mineq >= 0, "Mineq >= 0", srname);
                debug_obj.assert(n >= 1, "N >= 1", srname);
                debug_obj.assert(~ismember('Aineq', ipObj.UsingDefaults) == ~ismember('bineq', ipObj.UsingDefaults), "Aineq and Bineq are both present or both absent", srname);
                if ~ismember('Aineq', ipObj.UsingDefaults)
                    debug_obj.assert((size(Aineq, 1) == mineq && size(Aineq, 2) == n) || (size(Aineq, 1) == 0 && size(Aineq, 2) == 0 && mineq == 0), "SIZE(Aineq) == [Mineq, N] unless Aineq and Bineq are both empty", srname);
                end
                debug_obj.assert(~ismember('Aeq', ipObj.UsingDefaults) == ~ismember('beq', ipObj.UsingDefaults), "Aeq and Beq are both present or both absent", srname);
                if ~ismember('Aeq', ipObj.UsingDefaults)
                    debug_obj.assert((size(Aeq, 1) == meq && size(Aeq, 2) == n) || (size(Aeq, 1) == 0 && size(Aeq, 2) == 0 && meq == 0), "SIZE(Aeq) == [Meq, N] unless Aeq and Beq are both empty", srname);
                end
                if ~ismember('xl', ipObj.UsingDefaults)
                    debug_obj.assert(numel(xl) == n || numel(xl) == 0, "SIZE(XL) == N unless XL is empty", srname);
                end
                if ~ismember('xu', ipObj.UsingDefaults)
                    debug_obj.assert(numel(xu) == n || numel(xu) == 0, "SIZE(XU) == N unless XU is empty", srname);
                end
            end

            % Read the inputs

            x(:) = evaluate_obj.moderatex(x);

            Aineq_loc = memory_obj.alloc_rmatrix_sp(mineq, n); % NOT removable even in F2003, as Aineq may be absent or of size 0-by-0.
            if ~ismember('Aineq', ipObj.UsingDefaults) && mineq > 0
                % We must check Mineq > 0. Otherwise, the size of Aineq_LOC may be changed to 0-by-0 due to
                % automatic (re)allocation if that is the size of Aineq; we allow Aineq to be 0-by-0, but
                % Aineq_LOC should be n-by-0.
                Aineq_loc = Aineq;
            end

            bineq_loc = memory_obj.alloc_rvector_sp(mineq); % NOT removable even in F2003, as Bineq may be absent.
            if ~ismember('bineq', ipObj.UsingDefaults)
                bineq_loc = bineq;
            end

            Aeq_loc = memory_obj.alloc_rmatrix_sp(meq, n); % NOT removable even in F2003, as Aeq may be absent or of size 0-by-0.
            if ~ismember('Aeq', ipObj.UsingDefaults) && meq > 0
                % We must check Meq > 0. Otherwise, the size of Aeq_LOC may be changed to 0-by-0 due to
                % automatic (re)allocation if that is the size of Aeq; we allow Aeq to be 0-by-0, but
                % Aeq_LOC should be n-by-0.
                Aeq_loc = Aeq;
            end

            beq_loc = memory_obj.alloc_rvector_sp(meq); % NOT removable even in F2003, as Beq may be absent.
            if ~ismember('beq', ipObj.UsingDefaults)
                beq_loc = beq;
            end

            xl_loc(:) = -consts_obj.BOUNDMAX;
            if ~ismember('xl', ipObj.UsingDefaults)
                if numel(xl) > 0
                    xl_loc(:) = xl;
                end
            end
            xl_loc(linalg_obj.trueloc(infnan_obj.is_nan_sp(xl_loc) | xl_loc < -consts_obj.BOUNDMAX)) = -consts_obj.BOUNDMAX;

            xu_loc(:) = consts_obj.BOUNDMAX;
            if ~ismember('xu', ipObj.UsingDefaults)
                if numel(xu) > 0
                    xu_loc(:) = xu;
                end
            end
            xu_loc(linalg_obj.trueloc(infnan_obj.is_nan_sp(xu_loc) | xu_loc > consts_obj.BOUNDMAX)) = consts_obj.BOUNDMAX;

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
                if infnan_obj.is_finite(rhoend) && rhoend > 0
                    rhobeg_loc = max(consts_obj.TEN * rhoend, consts_obj.RHOBEG_DFT);
                else
                    rhobeg_loc = consts_obj.RHOBEG_DFT;
                end
            else
                rhobeg_loc = consts_obj.RHOBEG_DFT;
            end

            if ~ismember('rhoend', ipObj.UsingDefaults)
                rhoend_loc = rhoend;
            elseif rhobeg_loc > 0
                rhoend_loc = max(consts_obj.EPS, min((consts_obj.RHOEND_DFT / consts_obj.RHOBEG_DFT) * rhobeg_loc, consts_obj.RHOEND_DFT));
            else
                rhoend_loc = consts_obj.RHOEND_DFT;
            end

            if ismember('ctol', ipObj.UsingDefaults)
                ctol_loc = consts_obj.CTOL_DFT;
            else
                ctol_loc = ctol;
            end

            if ismember('cweight', ipObj.UsingDefaults)
                cweight_loc = consts_obj.CWEIGHT_DFT;
            else
                cweight_loc = cweight;
            end

            if ismember('ftarget', ipObj.UsingDefaults)
                ftarget_loc = consts_obj.FTARGET_DFT;
            else
                ftarget_loc = ftarget;
            end

            if ismember('maxfun', ipObj.UsingDefaults)
                maxfun_loc = consts_obj.MAXFUN_DIM_DFT * n;
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

            if ismember('iprint', ipObj.UsingDefaults)
                iprint_loc = consts_obj.IPRINT_DFT;
            else
                iprint_loc = iprint;
            end

            if ~ismember('eta1', ipObj.UsingDefaults)
                eta1_loc = eta1;
            elseif ~ismember('eta2', ipObj.UsingDefaults)
                if eta2 > 0 && eta2 < 1
                    eta1_loc = max(consts_obj.EPS, eta2 / 7.0);
                end
            else
                eta1_loc = consts_obj.TENTH;
            end

            if ~ismember('eta2', ipObj.UsingDefaults)
                eta2_loc = eta2;
            elseif eta1_loc > 0 && eta1_loc < 1
                eta2_loc = (eta1_loc + consts_obj.TWO) / 3.0;
            else
                eta2_loc = 0.7;
            end

            if ismember('gamma1', ipObj.UsingDefaults)
                gamma1_loc = consts_obj.HALF;
            else
                gamma1_loc = gamma1;
            end

            if ismember('gamma2', ipObj.UsingDefaults)
                gamma2_loc = consts_obj.TWO;
            else
                gamma2_loc = gamma2;
            end

            if ismember('maxhist', ipObj.UsingDefaults)
                maxhist_loc = max([maxfun_loc, n + 3, consts_obj.MAXFUN_DIM_DFT * n], [], 'all');
            else
                maxhist_loc = maxhist;
            end

            if ismember('maxfilt', ipObj.UsingDefaults)
                maxfilt_loc = consts_obj.MAXFILT_DFT;
            else
                maxfilt_loc = maxfilt;
            end

            % Preprocess the inputs in case some of them are invalid. It does nothing if all inputs are valid.
            [iprint_loc, maxfun_loc, maxhist_loc, ftarget_loc, rhobeg_loc, rhoend_loc, npt_loc, maxfilt_loc, ctol_loc, cweight_loc, eta1_loc, eta2_loc, gamma1_loc, gamma2_loc] = preproc_obj.preproc(solver, n, iprint_loc, maxfun_loc, maxhist_loc, ftarget_loc, rhobeg_loc, rhoend_loc, 'npt', npt_loc, 'ctol', ctol_loc, 'cweight', cweight_loc, 'eta1', eta1_loc, 'eta2', eta2_loc, 'gamma1', gamma1_loc, 'gamma2', gamma2_loc, 'maxfilt', maxfilt_loc);

            % Further revise MAXHIST_LOC according to MAXHISTMEM, and allocate memory for the history.
            % In MATLAB/Python/Julia/R implementation, we should simply set MAXHIST = MAXFUN and initialize
            % CHIST = NaN(1, MAXFUN), FHIST = NaN(1, MAXFUN), XHIST = NaN(N, MAXFUN)
            % if they are requested; replace MAXFUN with 0 for the history that is not requested.
            [maxhist_loc, xhist_loc, fhist_loc, chist_loc] = history_obj.prehist(maxhist_loc, n, nargout >= 5, nargout >= 6, 'output_chist', nargout >= 7);

            % Wrap the linear and bound constraints into a single constraint: AMAT^T*X <= BVEC.
            [amat, bvec] = obj.get_lincon(Aeq_loc, Aineq_loc, beq_loc, bineq_loc, rhoend_loc, xl_loc, xu_loc, x, amat, bvec);

            %-------------------- Call LINCOB, which performs the real calculations. --------------------------%
            if ismember('callback_fcn', ipObj.UsingDefaults)
                [x, nf_loc, chist_loc, cstrv_loc, f_loc, fhist_loc, xhist_loc, info_loc] = lincob_obj.lincob(calfun, iprint_loc, maxfilt_loc, maxfun_loc, npt_loc, Aeq_loc, Aineq_loc, amat, beq_loc, bineq_loc, bvec, ctol_loc, cweight_loc, eta1_loc, eta2_loc, ftarget_loc, gamma1_loc, gamma2_loc, rhobeg_loc, rhoend_loc, xl_loc, xu_loc, x, chist_loc, fhist_loc, xhist_loc);
            else
                [x, nf_loc, chist_loc, cstrv_loc, f_loc, fhist_loc, xhist_loc, info_loc] = lincob_obj.lincob(calfun, iprint_loc, maxfilt_loc, maxfun_loc, npt_loc, Aeq_loc, Aineq_loc, amat, beq_loc, bineq_loc, bvec, ctol_loc, cweight_loc, eta1_loc, eta2_loc, ftarget_loc, gamma1_loc, gamma2_loc, rhobeg_loc, rhoend_loc, xl_loc, xu_loc, x, chist_loc, fhist_loc, xhist_loc, 'callback_fcn', callback_fcn);
            end
            %--------------------------------------------------------------------------------------------------%

            % Deallocate variables not needed any more. We prefer explicit deallocation to the automatic one.



            % Write the outputs.



            % Copy XHIST_LOC to XHIST if needed.
            if nargout >= 5
                nhist = min(nf_loc, size(xhist_loc, 2));
                %----------------------------------------------------%
                xhist = memory_obj.alloc_rmatrix_sp(n, nhist); % Removable in F2003.
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
            if nargout >= 6
                nhist = min(nf_loc, numel(fhist_loc));
                %--------------------------------------------------%
                fhist = memory_obj.alloc_rvector_sp(nhist); % Removable in F2003.
                %--------------------------------------------------%
                fhist = fhist_loc(1:nhist); % The same as XHIST, we must cap FHIST at NF_LOC.

            end


            % Copy CHIST_LOC to CHIST if needed.
            if nargout >= 7
                nhist = min(nf_loc, numel(chist_loc));
                %--------------------------------------------------%
                chist = memory_obj.alloc_rvector_sp(nhist); % Removable in F2003.
                %--------------------------------------------------%
                chist = chist_loc(1:nhist); % The same as XHIST, we must cap CHIST at NF_LOC.

            end


            % If NF_LOC > MAXHIST_LOC, warn that not all history is recorded.
            if (nargout >= 5 || nargout >= 6 || nargout >= 7) && maxhist_loc < nf_loc
                debug_obj.warning(solver, "Only the history of the last " + string_obj.int2str(maxhist_loc) + " function evaluation(s) is recorded");
            end

            % Postconditions
            if consts_obj.DEBUGGING
                debug_obj.assert(nf_loc <= maxfun_loc, "NF <= MAXFUN", srname);
                debug_obj.assert(numel(x) == n && ~any(infnan_obj.is_nan_sp(x), 'all'), "SIZE(X) == N, X does not contain NaN", srname);
                nhist = min(nf_loc, maxhist_loc);
                if nargout >= 5
                    debug_obj.assert(size(xhist, 1) == n && size(xhist, 2) == nhist, "SIZE(XHIST) == [N, NHIST]", srname);
                    debug_obj.assert(~any(infnan_obj.is_nan_sp(xhist), 'all'), "XHIST does not contain NaN", srname);
                end
                if nargout >= 6
                    debug_obj.assert(numel(fhist) == nhist, "SIZE(FHIST) == NHIST", srname);
                    debug_obj.assert(~any(infnan_obj.is_nan_sp(fhist) | infnan_obj.is_posinf(fhist), 'all'), "FHIST does not contain NaN/+Inf", srname);
                end
                if nargout >= 7
                    debug_obj.assert(numel(chist) == nhist, "SIZE(CHIST) == NHIST", srname);
                    debug_obj.assert(~any(infnan_obj.is_nan_sp(chist) | infnan_obj.is_posinf(chist), 'all'), "CHIST does not contain NaN/+Inf", srname);
                end
                if nargout >= 6 && nargout >= 7
                    debug_obj.assert(~any(selectx_obj.isbetter10(fhist(1:nhist), chist(1:nhist), f_loc, cstrv_loc, ctol_loc), 'all'), "No point in the history is better than X", srname);
                end
            end

        end
        function [amat, bvec] = get_lincon(~, Aeq, Aineq, beq, bineq, rhoend, xl, xu, x0, amat, bvec)
            %--------------------------------------------------------------------------------------------------%
            % This subroutine wraps the linear and bound constraints into a single constraint: AMAT^T*X <= BVEC.
            % N.B.:
            % 1. LINCOA modifies the right hand sides of the constraints to make the starting point feasible if
            % it is not. This is not ideal, but Powell's code was implemented in this way. In the
            % MATLAB/Python/Julia/R code, we should include a preprocessing subroutine to project the starting
            % point to the feasible region if it is infeasible, so that the modification will not occur.
            % 2. The linear inequality constraints received by LINCOB is AMAT^T * X <= BVEC. Note that Each
            % column of AMAT corresponds to a constraint. This is different from Aineq and Aeq, whose rows
            % correspond to constraints. AMAT is defined in this way because it is accessed in columns during
            % the computation, and because Fortran saves arrays in the column-major order. In Python/C
            % implementations, AMAT should be transposed.
            % 3. LINCOA normalizes the linear constraints so that each constraint has a gradient of norm 1. This
            % is essential for LINCOA.
            %--------------------------------------------------------------------------------------------------%

            % Common modules
            consts_obj = prima_mat.common.consts_mod();
            debug_obj = prima_mat.common.debug_mod();
            linalg_obj = prima_mat.common.linalg_mod();
            memory_obj = prima_mat.common.memory_mod();


            % Inputs



            % Outputs



            % Local variables
            solver = "LINCOA";
            srname = "GET_LINCON";


            Aeq_norm = NaN(size(Aeq, 1), 1);
            Aeqx0 = NaN(size(Aeq, 1), 1);
            Aineq_norm = NaN(size(Aineq, 1), 1);
            Aineqx0 = NaN(size(Aineq, 1), 1);
            idmat = NaN(numel(x0));


            % Sizes
            n = numel(x0);

            % Preconditions
            if consts_obj.DEBUGGING
                debug_obj.assert(size(Aineq, 1) == numel(bineq) && size(Aineq, 2) == n, "SIZE(AINEQ) == [SIZE(BINEQ), N]", srname);
                debug_obj.assert(size(Aeq, 1) == numel(beq) && size(Aeq, 2) == n, "SIZE(AEQ) == [SIZE(BEQ), N]", srname);
                debug_obj.assert(numel(xl) == n && numel(xu) == n, "SIZE(XL) == SIZE(XU) == N", srname);
            end

            %====================%
            % Calculation starts %
            %====================%

            % Decide the number of nontrivial and valid (gradient is nonzero) constraints.
            mxl = nnz(xl > -consts_obj.BOUNDMAX);
            mxu = nnz(xu < consts_obj.BOUNDMAX);
            Aeq_norm(:) = sqrt(sum(Aeq .^ 2, 2));
            meq = nnz(Aeq_norm > 0);
            Aineq_norm(:) = sqrt(sum(Aineq .^ 2, 2));
            mineq = nnz(Aineq_norm > 0);
            m = mxl + mxu + 2 * meq + mineq; % The final number of linear inequality constraints.

            % Print a warning if some constraints are invalid. They will be ignored (Powell's code would stop).
            if meq < size(Aeq, 1) || mineq < size(Aineq, 1)
                debug_obj.warning(solver, "Some linear constraints have zero gradients; they are ignored");
            end

            % Allocate memory. Removable in F2003.
            memory_obj.alloc_ivector(mxl);
            memory_obj.alloc_ivector(mxu);
            memory_obj.alloc_ivector(meq);
            memory_obj.alloc_ivector(mineq);
            amat = memory_obj.alloc_rmatrix_sp(n, m);
            bvec = memory_obj.alloc_rvector_sp(m);
            memory_obj.alloc_rvector_sp(2 * meq + mineq);

            % Define the indices of the valid and nontrivial constraints.
            ixl = linalg_obj.trueloc(xl > -consts_obj.BOUNDMAX);
            ixu = linalg_obj.trueloc(xu < consts_obj.BOUNDMAX);
            ieq = linalg_obj.trueloc(Aeq_norm > 0);
            iineq = linalg_obj.trueloc(Aineq_norm > 0);

            % Wrap the linear constraints.
            % The bound constraint XL <= X <= XU is handled as two constraints -X <= -XL, X <= XU.
            % The equality constraint Aeq*X = Beq is handled as two constraints -Aeq*X <= -Beq, Aeq*X <= Beq.
            % N.B.:
            % 1. The treatment of the equality constraints is naive. One may choose to eliminate them instead.
            % 2. The code below is quite inefficient in terms of memory, but we prefer readability.
            idmat(:, :) = linalg_obj.eye2(n, n);
            amat = reshape([reshape(-idmat(:, ixl), 1, []), reshape(idmat(:, ixu), 1, []), reshape(-Aeq(ieq, :).', 1, []), reshape(Aeq(ieq, :).', 1, []), reshape(Aineq(iineq, :).', 1, [])], size(amat));
            bvec = [reshape(-xl(ixl), [], 1); reshape(xu(ixu), [], 1); reshape(-beq(ieq), [], 1); reshape(beq(ieq), [], 1); reshape(bineq(iineq), [], 1)];
            %%MATLAB code:
            %%amat = [-idmat(:, ixl), idmat(:, ixu), -Aeq(ieq, :)', Aeq(ieq, :)', Aineq(iineq, :)'];
            %%bvec = [-xl(ixl); xu(ixu); -beq(ieq); beq(ieq); bineq(iineq)];

            % Modify BVEC if necessary so that the initial point is feasible.
            Aeqx0(:) = linalg_obj.matprod21(Aeq, x0);
            Aineqx0(:) = linalg_obj.matprod21(Aineq, x0);
            bvec = max(bvec, [reshape(-x0(ixl), [], 1); reshape(x0(ixu), [], 1); reshape(-Aeqx0(ieq), [], 1); reshape(Aeqx0(ieq), [], 1); reshape(Aineqx0(iineq), [], 1)]);

            % Normalize the linear constraints so that each constraint has a gradient of norm 1.
            Anorm = [reshape(Aeq_norm(ieq), [], 1); reshape(Aeq_norm(ieq), [], 1); reshape(Aineq_norm(iineq), [], 1)];
            amat(:, mxl + mxu + 1:m) = amat(:, mxl + mxu + 1:m) ./ reshape(Anorm, 1, []);
            bvec(mxl + mxu + 1:m) = bvec(mxl + mxu + 1:m) ./ Anorm;

            % Deallocate memory.


            % Print a warning if the starting point is sufficiently infeasible and the constraints are modified.
            smallx = consts_obj.TEN ^ max(-6, -consts_obj.MAXPOW10) * rhoend;
            constr_modified = (any(x0 + smallx < xl, 'all') || any(x0 - smallx > xu, 'all') || any(abs(Aeqx0 - beq) > smallx * Aeq_norm, 'all') || any(Aineqx0 - bineq > smallx * Aineq_norm, 'all'));
            if constr_modified
                debug_obj.warning(solver, "The starting point is infeasible. " + solver + " modified the right-hand sides of the constraints to make it feasible");
            end

            %====================%
            %  Calculation ends  %
            %====================%

            % Postconditions
            if consts_obj.DEBUGGING
                debug_obj.assert(size(amat, 1) == numel(x0) && size(amat, 2) == numel(bvec), "SIZE(AMAT) == [SIZE(X), SIZE(BVEC)]", srname);
                debug_obj.assert(all(linalg_obj.matprod12(x0, amat) - bvec <= max(consts_obj.TEN ^ max(-12, -consts_obj.MAXPOW10), 100.0 * consts_obj.EPS) * (consts_obj.ONE + sum(abs(x0), 'all') + sum(abs(bvec), 'all')), 'all'), "The starting point is feasible", srname);
            end
        end

    end
end