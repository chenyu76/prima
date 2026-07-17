classdef preproc_mod
    %--------------------------------------------------------------------------------------------------%
    % PREPROC_MOD is a module that preprocesses the inputs.
    %
    % Coded by Zaikun ZHANG (www.zhangzk.net) based on Powell's code and papers.
    %
    % Started: July 2020
    %
    % Last Modified: Mon 06 Apr 2026 10:54:37 PM CST
    %--------------------------------------------------------------------------------------------------%

    methods
        function obj = preproc_mod()
            % N.B.:
            % 1. If all the inputs are valid, then PREPROC should do nothing.
            % 2. In PREPROC, we use VALIDATE instead of ASSERT, so that the parameters are validated even if we
            %    are not in debug mode.

        end
        function [iprint, maxfun, maxhist, ftarget, rhobeg, rhoend, npt, maxfilt, ctol, cweight, eta1, eta2, gamma1, gamma2, x0] = preproc(~, solver, n, iprint, maxfun, maxhist, ftarget, rhobeg, rhoend, varargin)
            %--------------------------------------------------------------------------------------------------%
            % This subroutine preprocesses the inputs. It does nothing to the inputs that are valid.
            %--------------------------------------------------------------------------------------------------%
            consts_obj = prima_mat.common.consts_mod();


            debug_obj = prima_mat.common.debug_mod();
            infnan_obj = prima_mat.common.infnan_mod();
            linalg_obj = prima_mat.common.linalg_mod();
            memory_obj = prima_mat.common.memory_mod();
            string_obj = prima_mat.common.string_mod();

            % Compulsory inputs



            % Optional inputs



            % Compulsory in-outputs



            % Optional in-outputs



            % Local variables
            srname = "PREPROC";

            % INTEGER(IK) may overflow if IK corresponds to the 16-bit integer.
            % INTEGER(IK) may overflow if IK corresponds to the 16-bit integer.



            lbx = false(n, 1);
            ubx = false(n, 1);


            x0_in = NaN(n, 1);

            % Preconditions
            ipObj = inputParser();
            addParameter(ipObj, 'm', NaN);
            addParameter(ipObj, 'npt', NaN);
            addParameter(ipObj, 'maxfilt', NaN);
            addParameter(ipObj, 'ctol', NaN);
            addParameter(ipObj, 'cweight', NaN);
            addParameter(ipObj, 'eta1', NaN);
            addParameter(ipObj, 'eta2', NaN);
            addParameter(ipObj, 'gamma1', NaN);
            addParameter(ipObj, 'gamma2', NaN);
            addParameter(ipObj, 'is_constrained', false);
            addParameter(ipObj, 'has_rhobeg', false);
            addParameter(ipObj, 'honour_x0', false);
            addParameter(ipObj, 'xl', NaN);
            addParameter(ipObj, 'xu', NaN);
            addParameter(ipObj, 'x0', NaN);
            parse(ipObj, varargin{:});
            m = ipObj.Results.m;
            npt = ipObj.Results.npt;
            maxfilt = ipObj.Results.maxfilt;
            ctol = ipObj.Results.ctol;
            cweight = ipObj.Results.cweight;
            eta1 = ipObj.Results.eta1;
            eta2 = ipObj.Results.eta2;
            gamma1 = ipObj.Results.gamma1;
            gamma2 = ipObj.Results.gamma2;
            is_constrained = ipObj.Results.is_constrained;
            has_rhobeg = ipObj.Results.has_rhobeg;
            honour_x0 = ipObj.Results.honour_x0;
            xl = ipObj.Results.xl;
            xu = ipObj.Results.xu;
            x0 = ipObj.Results.x0;
            if consts_obj.DEBUGGING
                debug_obj.validate(n >= 1, "N >= 1", srname);
                debug_obj.validate((~ismember('npt', ipObj.UsingDefaults) || nargout >= 7) == (string_obj.lower(solver) == "newuoa" || string_obj.lower(solver) == "bobyqa" || string_obj.lower(solver) == "lincoa"), "NPT is present if and only if SOLVER is NEWUOA, BOBYQA, or LINCOA", srname);
                if ~ismember('m', ipObj.UsingDefaults)
                    debug_obj.validate(m >= 0, "M >= 0", srname);
                    debug_obj.validate(m == 0 || string_obj.lower(solver) == "cobyla", "M == 0 unless the solver is COBYLA", srname);
                end
                if string_obj.lower(solver) == "cobyla" && ~ismember('m', ipObj.UsingDefaults) && ~ismember('is_constrained', ipObj.UsingDefaults)
                    debug_obj.validate(m == 0 || is_constrained, "For COBYLA, M == 0 unless the problem is constrained", srname);
                end
                debug_obj.validate((~ismember('maxfilt', ipObj.UsingDefaults) || nargout >= 8) == (string_obj.lower(solver) == "lincoa" || string_obj.lower(solver) == "cobyla"), "MAXFILT is present if and only if the solver is LINCOA or COBYLA", srname);
                if string_obj.lower(solver) == "bobyqa"
                    debug_obj.validate(~ismember('xl', ipObj.UsingDefaults) && ~ismember('xu', ipObj.UsingDefaults), "XL and XU are present if the solver is BOBYQA", srname);
                    debug_obj.validate(all(xu - xl >= consts_obj.TWO * consts_obj.EPS, 'all'), "MINVAL(XU-XL) > 2*EPS", srname);
                end
                debug_obj.validate((~ismember('honour_x0', ipObj.UsingDefaults) == (~ismember('x0', ipObj.UsingDefaults) || nargout >= 15)) && (~ismember('honour_x0', ipObj.UsingDefaults) == ~ismember('has_rhobeg', ipObj.UsingDefaults)), "HONOUR_X0, X0, and HAS_RHOBEG are present or absent simultaneously", srname);
                debug_obj.validate(~ismember('honour_x0', ipObj.UsingDefaults) == (string_obj.lower(solver) == "bobyqa"), "HONOUR_X0 is present if and only if the solver is BOBYQA", srname);
                % N.B.: LINCOA and COBYLA will have HONOUR_X0 as well if we intend to make them respect bounds.
                % %call validate(present(honour_x0) .eqv. &
                % %    & (lower(solver) == 'bobyqa' .or. lower(solver) == 'lincoa' .or. lower(solver) == 'cobyla'), &
                % %    & 'HONOUR_X0 is present if and only if the solver is BOBYQA, LINCOA, or COBYLA', srname)

            end

            %====================%
            % Calculation starts %
            %====================%

            % Read M, if necessary
            if string_obj.lower(solver) == "cobyla" && ~ismember('m', ipObj.UsingDefaults)
                m_loc = m;
            else
                m_loc = 0;
            end

            % Decide whether the problem is truly constrained
            if ismember('is_constrained', ipObj.UsingDefaults)
                is_constrained_loc = (m_loc > 0);
            else
                is_constrained_loc = is_constrained;
            end

            % Validate IPRINT
            if abs(iprint) > 3
                iprint_in = iprint;
                iprint = consts_obj.IPRINT_DFT;
                debug_obj.warning(solver, "Invalid IPRINT: " + string_obj.int2str(iprint_in) + "; it should be 0, 1, -1, 2, -2, 3, or -3; it is set to " + string_obj.int2str(iprint));
            end

            % Validate MAXFUN
            % N.B.: The INT(N), INT(N+1), and INT(N+2) below convert integers to the default integer kind,
            % which is the kind of MIN_MAXFUN. Fortran compilers may complain without the conversion. It is
            % not needed in Python/MATLAB/Julia/R.
            switch string_obj.lower(solver)
            case "uobyqa"
                min_maxfun = (fix(n + 1) * fix(n + 2)) / 2 + 1; % INT(*) avoids overflow when IK is 16-bit.
                min_maxfun_str = "(N+1)(N+2)/2 + 1";
            case "cobyla"
                min_maxfun = fix(n) + 2;
                min_maxfun_str = "N + 2";
            otherwise                % CASE ('NEWUOA', 'BOBYQA', 'LINCOA')
                min_maxfun = fix(n) + 3;
                min_maxfun_str = "N + 3";
            end
            if maxfun <= max(0, min_maxfun - 1)
                maxfun_in = maxfun;
                if maxfun > 0
                    maxfun = fix(min_maxfun);
                else                    % We assume that non-positive values of MAXFUN are produced by overflow.
                    maxfun = fix(max(min_maxfun, 10 ^ min(4, floor(log10(realmax(class(maxfun))))))); %%MATLAB: maxfun =  max(min_maxfun, 10^4);
                    % N.B.: Do NOT set MAXFUN to HUGE(MAXFUN), as it may cause overflow and infinite cycling
                    % when used as the upper bound of DO loops. This occurred on 20240225 with gfortran 13. See
                    % https://fortran-lang.discourse.group/t/loop-variable-reaching-integer-huge-causes-infinite-loop
                    % https://fortran-lang.discourse.group/t/loops-dont-behave-like-they-should
                end
                debug_obj.warning(solver, "Invalid MAXFUN: " + string_obj.int2str(maxfun_in) + "; it should be at least " + min_maxfun_str + " with N = " + string_obj.int2str(n) + "; it is set to " + string_obj.int2str(maxfun));
            end

            % Validate MAXHIST
            if maxhist <= 0
                maxhist_in = maxhist;
                maxhist = maxfun;
                debug_obj.warning(solver, "Invalid MAXHIST: " + string_obj.int2str(maxhist_in) + "; it should be a positive integer; it is set to " + string_obj.int2str(maxhist));
            end
            maxhist = min(maxhist, maxfun); % MAXHIST > MAXFUN is never needed.

            % Validate FTARGET
            if infnan_obj.is_nan_sp(ftarget)
                % No warning if FTARGET is NaN, which is interpreted as no target function value is provided.
                ftarget = -realmax;
            end

            % Validate NPT
            if ~ismember('npt', ipObj.UsingDefaults) || nargout >= 7
                if npt < n + 2 || npt >= maxfun || 2 * fix(npt) > fix(n + 2) * fix(n + 1)
                    %INT(*) avoids overflow when IK is 16-bit
                    npt_in = npt;
                    npt = fix(min(maxfun - 1, 2 * n + 1));
                    debug_obj.warning(solver, "Invalid NPT: " + string_obj.int2str(npt_in) + "; it should be an integer in the interval [N+2, (N+1)(N+2)/2] with N = " + string_obj.int2str(n) + " and less than MAXFUN = " + string_obj.int2str(maxfun) + "; it is set to " + string_obj.int2str(npt));
                end
            end

            % Validate MAXFILT
            if ~ismember('maxfilt', ipObj.UsingDefaults) || nargout >= 8
                maxfilt_in = maxfilt;
                if maxfilt <= 0
                    maxfilt = consts_obj.MAXFILT_DFT;
                else
                    maxfilt = max(consts_obj.MIN_MAXFILT, maxfilt); % The inputted MAXFILT is too small.
                end
                % Further revise MAXFILT according to MAXHISTMEM.
                switch string_obj.lower(solver)
                case "lincoa"
                    unit_memo = fix(n + 2) * fix(memory_obj.size_of_sp(0.0)); % INT(*) avoids overflow when IK is 16-bit.
                case "cobyla"
                    unit_memo = fix(m_loc + n + 2) * fix(memory_obj.size_of_sp(0.0)); % INT(*) avoids overflow when IK is 16-bit.
                otherwise                    % The following should not be reached unless there is a bug, but we keep it for safety.
                    unit_memo = 1;
                end
                % We cannot simply set MAXFILT = MIN(MAXFILT, MAXHISTMEM/...), as they may not have
                % the same kind, and compilers may complain. We may convert them, but overflow may occur.
                if maxfilt > consts_obj.MAXHISTMEM / unit_memo
                    maxfilt = fix(consts_obj.MAXHISTMEM / unit_memo); % Integer division.

                end
                maxfilt = min(maxfun, max(consts_obj.MIN_MAXFILT, maxfilt));
                if is_constrained_loc
                    if maxfilt_in <= 0
                        debug_obj.warning(solver, "Invalid MAXFILT: " + string_obj.int2str(maxfilt_in) + "; it should be a positive integer; it is set to " + string_obj.int2str(maxfilt));
                    elseif maxfilt_in < min(maxfun, consts_obj.MIN_MAXFILT)
                        debug_obj.warning(solver, "MAXFILT = " + string_obj.int2str(maxfilt_in) + " is too small; it is set to " + string_obj.int2str(maxfilt));
                    elseif maxfilt < min(maxfilt_in, maxfun)
                        debug_obj.warning(solver, "MAXFILT is reduced from " + string_obj.int2str(maxfilt_in) + " to " + string_obj.int2str(maxfilt) + " due to memory limit");
                    end
                end
            end

            % Validate ETA1 and ETA2
            if ~(eta1 >= 0 && eta1 < 1)
                % ETA1 = NaN falls into this case.
                eta1_in = eta1;
                % Take ETA2 into account if it has a valid value.
                if eta2 >= 0 && eta2 < 1
                    eta1 = eta2 / 7.0;
                else
                    eta1 = consts_obj.ETA1_DFT;
                end
                debug_obj.warning(solver, "Invalid ETA1: " + string_obj.real2str_scalar(eta1_in) + "; it should be in the interval [0, 1) and not more than ETA2 = " + string_obj.real2str_scalar(eta2) + "; it is set to " + string_obj.real2str_scalar(eta1));
            end

            if ~(eta2 >= eta1 && eta2 < 1)
                % ETA2 = NaN falls into this case.
                eta2_in = eta2;
                % Take ETA1 into account if it has a valid value.
                if eta1 >= 0 && eta1 < 1
                    eta2 = (eta1 + consts_obj.TWO) / 3.0;
                else
                    eta2 = consts_obj.ETA2_DFT;
                end
                debug_obj.warning(solver, "Invalid ETA2: " + string_obj.real2str_scalar(eta2_in) + "; it should be in the interval [0, 1) and not less than ETA1 = " + string_obj.real2str_scalar(eta1) + "; it is set to " + string_obj.real2str_scalar(eta2));
            end

            % The following revision may update ETA1 slightly. It prevents ETA1 > ETA2 due to rounding
            % errors, which would not be accepted by the solvers.
            eta1 = min(eta1, eta2);

            % Validate GAMMA1 and GAMMA2
            if ~(gamma1 > 0 && gamma1 < 1)
                % GAMMA1 = NaN falls into this case.
                gamma1_in = gamma1;
                gamma1 = consts_obj.GAMMA1_DFT;
                debug_obj.warning(solver, "Invalid GAMMA1: " + string_obj.real2str_scalar(gamma1_in) + "; it should in the interval (0, 1); it is set to " + string_obj.real2str_scalar(gamma1));
            end

            if ~(infnan_obj.is_finite(gamma2) && gamma2 >= 1)
                % GAMMA2 = NaN falls into this case.
                gamma2_in = gamma2;
                gamma2 = consts_obj.GAMMA2_DFT;
                debug_obj.warning(solver, "Invalid GAMMA2: " + string_obj.real2str_scalar(gamma2_in) + "; it should be a real number not less than 1; it is set to " + string_obj.real2str_scalar(gamma2));
            end

            % Validate RHOBEG and RHOEND

            rhobeg_in = rhobeg;
            rhoend_in = rhoend;

            % Revise the default values for RHOBEG/RHOEND according to the solver.
            if string_obj.lower(solver) == "bobyqa"
                rhobeg_default = max(consts_obj.EPS, min(consts_obj.RHOBEG_DFT, min(xu - xl, [], 'all') / 4.0));
                rhoend_default = max(consts_obj.EPS, min((consts_obj.RHOEND_DFT / consts_obj.RHOBEG_DFT) * rhobeg_default, consts_obj.RHOEND_DFT));
            else
                rhobeg_default = consts_obj.RHOBEG_DFT;
                rhoend_default = consts_obj.RHOEND_DFT;
            end

            if string_obj.lower(solver) == "bobyqa"
                % Do NOT merge the IF below into the ELSEIF above! Otherwise, XU and XL may be accessed even if
                % the solver is not BOBYQA, because the logical evaluation is not short-circuit.
                if rhobeg > min(xu - xl, [], 'all') / consts_obj.TWO
                    % Do NOT make this revision if RHOBEG not positive or not finite, because otherwise RHOBEG
                    % will get a huge value when XU or XL contains huge values that indicate unbounded variables.
                    rhobeg = min(xu - xl, [], 'all') / 4.0; % Here, we do not take RHOBEG_DEFAULT.
                    debug_obj.warning(solver, "Invalid RHOBEG: " + string_obj.real2str_scalar(rhobeg_in) + "; " + solver + " requires 0 < RHOBEG <= MINVAL(XU-XL)/2 = " + string_obj.real2str_scalar(min(xu - xl, [], 'all') / 2.0) + "; it is set to " + string_obj.real2str_scalar(rhobeg));
                end
            end

            if ~(infnan_obj.is_finite(rhobeg) && rhobeg > 0)
                % RHOBEG = NaN falls into this case.
                % Take RHOEND into account if it has a valid value. We do not do this if the solver is BOBYQA,
                % which requires that RHOBEG <= (XU-XL)/2.
                if infnan_obj.is_finite(rhoend) && rhoend > 0 && string_obj.lower(solver) ~= "bobyqa"
                    rhobeg = max(consts_obj.TEN * rhoend, rhobeg_default);
                else
                    rhobeg = rhobeg_default;
                end
                debug_obj.warning(solver, "Invalid RHOBEG: " + string_obj.real2str_scalar(rhobeg_in) + "; it should be a positive number; it is set to " + string_obj.real2str_scalar(rhobeg));
            end

            if ~(infnan_obj.is_finite(rhoend) && rhoend >= 0 && rhoend <= rhobeg)
                % RHOEND = NaN falls into this case.
                rhoend = max(consts_obj.EPS, min((consts_obj.RHOEND_DFT / consts_obj.RHOBEG_DFT) * rhobeg, rhoend_default));
                debug_obj.warning(solver, "Invalid RHOEND: " + string_obj.real2str_scalar(rhoend_in) + "; we should have " + string_obj.real2str_scalar(rhobeg) + " = RHOBEG >= RHOEND >= 0; it is set to " + string_obj.real2str_scalar(rhoend));
            end

            % For BOBYQA, revise X0 or RHOBEG so that the distance between X0 and the inactive bounds is at
            % least RHOBEG. If HONOUR_X0 == FALSE, revise X0 if needed; then revise RHOBEG if needed.
            % N.B.: We should do the same for LINCOA and COBYLA if we make them respect the bounds in the future.
            % %if (lower(solver) == 'bobyqa' .or. lower(solver) == 'lincoa' .or. lower(solver) == 'cobyla') then
            if string_obj.lower(solver) == "bobyqa"
                % Revise X0 if allowed and needed.
                if ~honour_x0
                    x0_in(:) = x0; % Recorded to see whether X0 is really revised.
                    % N.B.: The following revision is valid only if XL <= X0 <= XU and RHOBEG <= MINVAL(XU-XL)/2,
                    % which should hold at this point due to the revision of RHOBEG and moderation of X0.
                    % The cases below are mutually exclusive in precise arithmetic as MINVAL(XU-XL) >= 2*RHOBEG.
                    mask00 = x0 <= xl + consts_obj.HALF * rhobeg; %Unsupported statement inside WHERE block: StatementLineBreak 1
                    x0(mask00) = xl(mask00); %Unsupported statement inside WHERE block: StatementLineBreak 1
                    mask01 = ~mask00 & x0 < xl + rhobeg; %Unsupported statement inside WHERE block: StatementLineBreak 1
                    x0(mask01) = xl(mask01) + rhobeg; %Unsupported statement inside WHERE block: StatementLineBreak 1

                    mask00 = x0 >= xu - consts_obj.HALF * rhobeg; %Unsupported statement inside WHERE block: StatementLineBreak 1
                    x0(mask00) = xu(mask00); %Unsupported statement inside WHERE block: StatementLineBreak 1
                    mask01 = ~mask00 & x0 > xu - rhobeg; %Unsupported statement inside WHERE block: StatementLineBreak 1
                    x0(mask01) = xu(mask01) - rhobeg; %Unsupported statement inside WHERE block: StatementLineBreak 1

                    %%MATLAB code:
                    %%lbx = (x0 <= xl + 0.5 * rhobeg);
                    %%lbx_plus = (x0 > xl + 0.5 * rhobeg .and. x0 < xl + rhobeg);
                    %%ubx = (x0 >= xu - 0.5 * rhobeg);
                    %%ubx_minus = (x0 < xu - 0.5 * rhobeg .and. x0 > xu - rhobeg);
                    %%x0(lbx) = xl(lbx);
                    %%x0(lbx_plus) = xl(lbx_plus) + rhobeg;
                    %%x0(ubx) = xu(ubx);
                    %%x0(ubx_minus) = xu(ubx_minus) - rhobeg;

                    if any(abs(x0_in - x0) > 0, 'all')
                        debug_obj.warning(solver, "X0 is revised so that the distance between X0 and the inactive bounds is at least RHOBEG = " + string_obj.real2str_scalar(rhobeg) + "; revise RHOBEG or set HONOUR_X0 to .TRUE. if you prefer to keep X0 unchanged");
                    end
                end

                % Revise RHOBEG if needed.
                % N.B.: If X0 has been revised above (i.e., HONOUR_X0 is FALSE), then the following revision
                % is unnecessary in precise arithmetic. However, it may still be needed due to rounding errors.
                lbx(:) = (infnan_obj.is_finite(xl) & x0 - xl <= consts_obj.EPS * max(consts_obj.ONE, abs(xl))); % X0 essentially equals XL
                ubx(:) = (infnan_obj.is_finite(xu) & x0 - xu >= -consts_obj.EPS * max(consts_obj.ONE, abs(xu))); % X0 essentially equals XU
                x0(linalg_obj.trueloc(lbx)) = xl(linalg_obj.trueloc(lbx));
                x0(linalg_obj.trueloc(ubx)) = xu(linalg_obj.trueloc(ubx));
                rhobeg = max(consts_obj.EPS, min([rhobeg; reshape(x0(linalg_obj.falseloc(lbx)) - xl(linalg_obj.falseloc(lbx)), [], 1); reshape(xu(linalg_obj.falseloc(ubx)) - x0(linalg_obj.falseloc(ubx)), [], 1)], [], 'all'));
                if rhobeg_in - rhobeg > consts_obj.EPS * max(consts_obj.ONE, rhobeg_in)
                    rhoend = max(consts_obj.EPS, min((rhoend / rhobeg_in) * rhobeg, rhoend)); % We do not revise RHOEND unless RHOBEG is truly revised.
                    if has_rhobeg
                        debug_obj.warning(solver, "RHOBEG is revised from " + string_obj.real2str_scalar(rhobeg_in) + " to " + string_obj.real2str_scalar(rhobeg) + " and RHOEND from " + string_obj.real2str_scalar(rhoend_in) + " to " + string_obj.real2str_scalar(rhoend) + " so that the distance between X0 and the inactive bounds is at least RHOBEG");
                    end
                end
            end

            % The following revision may update RHOBEG and RHOEND slightly. It particularly prevents
            % RHOEND > RHOBEG due to rounding errors, which would not be accepted by the solvers.
            rhobeg = max(rhobeg, consts_obj.EPS);
            rhoend = min(max(rhoend, consts_obj.EPS), rhobeg);

            % Validate CTOL (it can be 0)
            if ~ismember('ctol', ipObj.UsingDefaults) || nargout >= 9
                if ~(ctol >= 0)
                    % CTOL = NaN falls into this case.
                    ctol_in = ctol;
                    ctol = consts_obj.CTOL_DFT;
                    if is_constrained_loc
                        debug_obj.warning(solver, "Invalid CTOL: " + string_obj.real2str_scalar(ctol_in) + "; it should be a nonnegative number; it is set to " + string_obj.real2str_scalar(ctol));
                    end
                end
            end

            % Validate CWEIGHT (it can be +Inf)
            if ~ismember('cweight', ipObj.UsingDefaults) || nargout >= 10
                if ~(cweight >= 0)
                    % CWEIGHT = NaN falls into this case.
                    cweight_in = cweight;
                    cweight = consts_obj.CWEIGHT_DFT;
                    if is_constrained_loc
                        debug_obj.warning(solver, "Invalid CWEIGHT: " + string_obj.real2str_scalar(cweight_in) + "; it should be a nonnegative number; it is set to " + string_obj.real2str_scalar(cweight));
                    end
                end
            end

            %====================%
            %  Calculation ends  %
            %====================%

            % Postconditions
            if consts_obj.DEBUGGING
                debug_obj.validate(abs(iprint) <= 3, "IPRINT is 0, 1, -1, 2, -2, 3, or -3", solver);
                debug_obj.validate(maxhist >= 0 && maxhist <= maxfun, "0 <= MAXHIST <= MAXFUN", solver);
                debug_obj.validate(maxfun >= min_maxfun, "MAXFUN >= MIN_MAXFUN", solver);
                if ~ismember('npt', ipObj.UsingDefaults) || nargout >= 7
                    debug_obj.validate(npt >= n + 2 && npt < maxfun && 2 * fix(npt) <= fix(n + 2) * fix(n + 1), "N+2 <= NPT < MAXFUN and 2*NPT <= (N+1)(N+2)", solver);
                end
                if ~ismember('maxfilt', ipObj.UsingDefaults) || nargout >= 8
                    debug_obj.validate(maxfilt >= min(consts_obj.MIN_MAXFILT, maxfun) && maxfilt <= maxfun, "MIN(MIN_MAXFILT, MAXFUN) <= MAXFILT <= MAXFUN", solver);
                end
                debug_obj.validate(eta1 >= 0 && eta1 <= eta2 && eta2 < 1, "0 <= ETA1 <= ETA2 < 1", solver);
                debug_obj.validate(gamma1 > 0 && gamma1 < 1 && gamma2 > 1, "0 < GAMMA1 < 1 < GAMMA2", solver);
                debug_obj.validate(rhobeg >= rhoend && rhoend > 0, "RHOBEG >= RHOEND > 0", solver);
                if string_obj.lower(solver) == "bobyqa"
                    debug_obj.validate(all(rhobeg <= (xu - xl) ./ consts_obj.TWO, 'all'), "RHOBEG <= MINVAL(XU-XL)/2", solver);
                    debug_obj.validate(all(infnan_obj.is_finite(x0), 'all'), "X0 is finite", solver);
                    debug_obj.validate(all(x0 >= xl & (x0 <= xl | x0 - xl >= rhobeg), 'all'), "X0 == XL or X0 - XL >= RHOBEG", solver);
                    debug_obj.validate(all(x0 <= xu & (x0 >= xu | xu - x0 >= rhobeg), 'all'), "X0 == XU or XU - X0 >= RHOBEG", solver);
                end
                if ~ismember('ctol', ipObj.UsingDefaults) || nargout >= 9
                    debug_obj.validate(ctol >= 0, "CTOL >= 0", solver);
                end
            end

        end

    end
end