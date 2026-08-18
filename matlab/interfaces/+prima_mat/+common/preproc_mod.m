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



            % Compulsory inputs



            % Optional inputs



            % Compulsory in-outputs



            % Optional in-outputs



            % Local variables


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


            %====================%
            % Calculation starts %
            %====================%

            % Read M, if necessary
            if lower(solver) == "cobyla" && ~ismember('m', ipObj.UsingDefaults)
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

                iprint = 0;
            end

            % Validate MAXFUN
            % N.B.: The INT(N), INT(N+1), and INT(N+2) below convert integers to the default integer kind,
            % which is the kind of MIN_MAXFUN. Fortran compilers may complain without the conversion. It is
            % not needed in Python/MATLAB/Julia/R.
            switch lower(solver)
            case "uobyqa"
                min_maxfun = ((n + 1) * (n + 2)) / 2 + 1; % INT(*) avoids overflow when IK is 16-bit.

            case "cobyla"
                min_maxfun = n + 2;

            otherwise                % CASE ('NEWUOA', 'BOBYQA', 'LINCOA')
                min_maxfun = n + 3;

            end
            if maxfun <= max(0, min_maxfun - 1)

                if maxfun > 0
                    maxfun = min_maxfun;
                else                    % We assume that non-positive values of MAXFUN are produced by overflow.
                    maxfun = fix(max(min_maxfun, 10 ^ min(4, floor(log10(realmax(class(maxfun))))))); %%MATLAB: maxfun =  max(min_maxfun, 10^4);
                    % N.B.: Do NOT set MAXFUN to HUGE(MAXFUN), as it may cause overflow and infinite cycling
                    % when used as the upper bound of DO loops. This occurred on 20240225 with gfortran 13. See
                    % https://fortran-lang.discourse.group/t/loop-variable-reaching-integer-huge-causes-infinite-loop
                    % https://fortran-lang.discourse.group/t/loops-dont-behave-like-they-should
                end
            end

            % Validate MAXHIST
            if maxhist <= 0

                maxhist = maxfun;
            end
            maxhist = min(maxhist, maxfun); % MAXHIST > MAXFUN is never needed.

            % Validate FTARGET
            if isnan(ftarget)
                % No warning if FTARGET is NaN, which is interpreted as no target function value is provided.
                ftarget = -realmax;
            end

            % Validate NPT
            if ~ismember('npt', ipObj.UsingDefaults) || nargout >= 7
                if npt < n + 2 || npt >= maxfun || 2 * npt > (n + 2) * (n + 1)
                    %INT(*) avoids overflow when IK is 16-bit

                    npt = min(maxfun - 1, 2 * n + 1);
                end
            end

            % Validate MAXFILT
            if ~ismember('maxfilt', ipObj.UsingDefaults) || nargout >= 8
                maxfilt_in = maxfilt;
                if maxfilt <= 0
                    maxfilt = 2000;
                else
                    maxfilt = max(200, maxfilt); % The inputted MAXFILT is too small.
                end
                % Further revise MAXFILT according to MAXHISTMEM.
                switch lower(solver)
                case "lincoa"
                    unit_memo = (n + 2) * fix(8); % INT(*) avoids overflow when IK is 16-bit.
                case "cobyla"
                    unit_memo = (m_loc + n + 2) * fix(8); % INT(*) avoids overflow when IK is 16-bit.
                otherwise                    % The following should not be reached unless there is a bug, but we keep it for safety.
                    unit_memo = 1;
                end
                % We cannot simply set MAXFILT = MIN(MAXFILT, MAXHISTMEM/...), as they may not have
                % the same kind, and compilers may complain. We may convert them, but overflow may occur.
                if maxfilt > 100000000 / unit_memo
                    maxfilt = fix(100000000 / unit_memo); % Integer division.

                end
                maxfilt = min(maxfun, max(200, maxfilt));
                if is_constrained_loc
                    if maxfilt_in <= 0
                    elseif maxfilt_in < min(maxfun, 200)
                    elseif maxfilt < min(maxfilt_in, maxfun)
                    end
                end
            end

            % Validate ETA1 and ETA2
            if ~(eta1 >= 0 && eta1 < 1)
                % ETA1 = NaN falls into this case.

                % Take ETA2 into account if it has a valid value.
                if eta2 >= 0 && eta2 < 1
                    eta1 = eta2 / 7.0;
                else
                    eta1 = 0.1;
                end
            end

            if ~(eta2 >= eta1 && eta2 < 1)
                % ETA2 = NaN falls into this case.

                % Take ETA1 into account if it has a valid value.
                if eta1 >= 0 && eta1 < 1
                    eta2 = (eta1 + 2.0) / 3.0;
                else
                    eta2 = 0.7;
                end
            end

            % The following revision may update ETA1 slightly. It prevents ETA1 > ETA2 due to rounding
            % errors, which would not be accepted by the solvers.
            eta1 = min(eta1, eta2);

            % Validate GAMMA1 and GAMMA2
            if ~(gamma1 > 0 && gamma1 < 1)
                % GAMMA1 = NaN falls into this case.

                gamma1 = 0.5;
            end

            if ~(isfinite(gamma2) && gamma2 >= 1)
                % GAMMA2 = NaN falls into this case.

                gamma2 = 2.0;
            end

            % Validate RHOBEG and RHOEND

            rhobeg_in = rhobeg;


            % Revise the default values for RHOBEG/RHOEND according to the solver.
            if lower(solver) == "bobyqa"
                rhobeg_default = max(eps(1.0), min(1.0, min(xu - xl, [], 'all') / 4.0));
                rhoend_default = max(eps(1.0), min((1.0e-6 / 1.0) * rhobeg_default, 1.0e-6));
            else
                rhobeg_default = 1.0;
                rhoend_default = 1.0e-6;
            end

            if lower(solver) == "bobyqa"
                % Do NOT merge the IF below into the ELSEIF above! Otherwise, XU and XL may be accessed even if
                % the solver is not BOBYQA, because the logical evaluation is not short-circuit.
                if rhobeg > min(xu - xl, [], 'all') / 2.0
                    % Do NOT make this revision if RHOBEG not positive or not finite, because otherwise RHOBEG
                    % will get a huge value when XU or XL contains huge values that indicate unbounded variables.
                    rhobeg = min(xu - xl, [], 'all') / 4.0; % Here, we do not take RHOBEG_DEFAULT.

                end
            end

            if ~(isfinite(rhobeg) && rhobeg > 0)
                % RHOBEG = NaN falls into this case.
                % Take RHOEND into account if it has a valid value. We do not do this if the solver is BOBYQA,
                % which requires that RHOBEG <= (XU-XL)/2.
                if isfinite(rhoend) && rhoend > 0 && lower(solver) ~= "bobyqa"
                    rhobeg = max(10.0 * rhoend, rhobeg_default);
                else
                    rhobeg = rhobeg_default;
                end
            end

            if ~(isfinite(rhoend) && rhoend >= 0 && rhoend <= rhobeg)
                % RHOEND = NaN falls into this case.
                rhoend = max(eps(1.0), min((1.0e-6 / 1.0) * rhobeg, rhoend_default));
            end

            % For BOBYQA, revise X0 or RHOBEG so that the distance between X0 and the inactive bounds is at
            % least RHOBEG. If HONOUR_X0 == FALSE, revise X0 if needed; then revise RHOBEG if needed.
            % N.B.: We should do the same for LINCOA and COBYLA if we make them respect the bounds in the future.
            % %if (lower(solver) == 'bobyqa' .or. lower(solver) == 'lincoa' .or. lower(solver) == 'cobyla') then
            if lower(solver) == "bobyqa"
                % Revise X0 if allowed and needed.
                if ~honour_x0
                    x0_in(:) = x0; % Recorded to see whether X0 is really revised.
                    % N.B.: The following revision is valid only if XL <= X0 <= XU and RHOBEG <= MINVAL(XU-XL)/2,
                    % which should hold at this point due to the revision of RHOBEG and moderation of X0.
                    % The cases below are mutually exclusive in precise arithmetic as MINVAL(XU-XL) >= 2*RHOBEG.
                    mask00 = x0 <= xl + 0.5 * rhobeg; %Unsupported statement inside WHERE block: StatementLineBreak 1
                    x0(mask00) = xl(mask00); %Unsupported statement inside WHERE block: StatementLineBreak 1
                    mask01 = ~mask00 & x0 < xl + rhobeg; %Unsupported statement inside WHERE block: StatementLineBreak 1
                    x0(mask01) = xl(mask01) + rhobeg; %Unsupported statement inside WHERE block: StatementLineBreak 1

                    mask00 = x0 >= xu - 0.5 * rhobeg; %Unsupported statement inside WHERE block: StatementLineBreak 1
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
                    end
                end

                % Revise RHOBEG if needed.
                % N.B.: If X0 has been revised above (i.e., HONOUR_X0 is FALSE), then the following revision
                % is unnecessary in precise arithmetic. However, it may still be needed due to rounding errors.
                lbx(:) = (isfinite(xl) & x0 - xl <= eps(1.0) * max(1.0, abs(xl))); % X0 essentially equals XL
                ubx(:) = (isfinite(xu) & x0 - xu >= -eps(1.0) * max(1.0, abs(xu))); % X0 essentially equals XU
                x0(lbx) = xl(lbx);
                x0(ubx) = xu(ubx);
                rhobeg = max(eps(1.0), min([rhobeg; x0(find(~lbx)) - xl(find(~lbx)); xu(find(~ubx)) - x0(find(~ubx))], [], 'all'));
                if rhobeg_in - rhobeg > eps(1.0) * max(1.0, rhobeg_in)
                    rhoend = max(eps(1.0), min((rhoend / rhobeg_in) * rhobeg, rhoend)); % We do not revise RHOEND unless RHOBEG is truly revised.
                    if has_rhobeg
                    end
                end
            end

            % The following revision may update RHOBEG and RHOEND slightly. It particularly prevents
            % RHOEND > RHOBEG due to rounding errors, which would not be accepted by the solvers.
            rhobeg = max(rhobeg, eps(1.0));
            rhoend = min(max(rhoend, eps(1.0)), rhobeg);

            % Validate CTOL (it can be 0)
            if ~ismember('ctol', ipObj.UsingDefaults) || nargout >= 9
                if ~(ctol >= 0)
                    % CTOL = NaN falls into this case.

                    ctol = sqrt(eps(1.0));
                    if is_constrained_loc
                    end
                end
            end

            % Validate CWEIGHT (it can be +Inf)
            if ~ismember('cweight', ipObj.UsingDefaults) || nargout >= 10
                if ~(cweight >= 0)
                    % CWEIGHT = NaN falls into this case.

                    cweight = 1.0e8;
                    if is_constrained_loc
                    end
                end
            end

            %====================%
            %  Calculation ends  %
            %====================%

            % Postconditions


        end

    end
end