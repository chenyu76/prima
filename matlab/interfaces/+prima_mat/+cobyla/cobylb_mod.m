% TODO: Implement GETMODEL to get the model of the objective function and constraints, i.e., g and A.
classdef cobylb_mod
    %--------------------------------------------------------------------------------------------------%
    % This module performs the major calculations of COBYLA.
    %
    % Coded by Zaikun ZHANG (www.zhangzk.net) based on Powell's code and the COBYLA paper.
    %
    % N.B. (Zaikun 20220131): Powell's implementation of COBYLA uses RHO rather than DELTA as the
    % trust-region radius, and RHO is never increased. DELTA does not exist in Powell's COBYLA code.
    % Following the idea in Powell's other solvers (UOBYQA, ..., LINCOA), our code uses DELTA as the
    % trust-region radius, while RHO works a lower bound of DELTA and indicates the current resolution
    % of the algorithm. DELTA is updated in a classical way subject to DELTA >= RHO, whereas RHO is
    % updated as in Powell's COBYLA code and is never increased. The new implementation improves the
    % performance of COBYLA.
    %
    % Dedicated to the late Professor M. J. D. Powell FRS (1936--2015).
    %
    % Started: July 2021
    %
    % Last Modified: Wed 08 Apr 2026 06:38:14 PM CST
    %--------------------------------------------------------------------------------------------------%

    methods
        function [constr, f, x, nf, chist, conhist, cstrv, fhist, xhist, info] = cobylb(obj, calcfc, iprint, maxfilt, maxfun, amat, bvec, ctol, cweight, eta1, eta2, ftarget, gamma1, gamma2, rhobeg, rhoend, constr, f, x, chist, conhist, fhist, xhist, varargin)
            %--------------------------------------------------------------------------------------------------%
            % This subroutine performs the actual calculations of COBYLA.
            %
            % IPRINT, MAXFILT, MAXFUN, MAXHIST, CTOL, CWEIGHT, ETA1, ETA2, FTARGET, GAMMA1, GAMMA2, RHOBEG,
            % RHOEND, X, NF, F, XHIST, FHIST, CHIST, CONHIST, CSTRV, INFO and CALLBACK are identical to the corresponding
            % arguments in subroutine COBYLA.
            %--------------------------------------------------------------------------------------------------%


            checkexit_obj = prima_mat.common.checkexit_mod();

            evaluate_obj = prima_mat.common.evaluate_mod();
            history_obj = prima_mat.common.history_mod();

            message_obj = prima_mat.common.message_mod();

            ratio_obj = prima_mat.common.ratio_mod();
            redrho_obj = prima_mat.common.redrho_mod();
            selectx_obj = prima_mat.common.selectx_mod();

            % Solver-specific modules
            geometry_cobyla_obj = prima_mat.cobyla.geometry_cobyla_mod();
            initialize_cobyla_obj = prima_mat.cobyla.initialize_cobyla_mod();
            trustregion_cobyla_obj = prima_mat.cobyla.trustregion_cobyla_mod();
            update_cobyla_obj = prima_mat.cobyla.update_cobyla_mod();

            % On entry, [X, F, CONSTR] = [X0, F(X0), CONSTR(X0)]


            solver = "COBYLA";

            j = NaN;

            bad_trstep = false;
            adequate_geo = false;
            evaluated = false(numel(x) + 1, 1);
            improve_geo = false;
            reduce_rho = false;

            ximproved = false;
            A = NaN(numel(x), numel(constr)); % A contains the approximate gradient for the constraints

            cfilt = NaN(min(max(maxfilt, 1), maxfun), 1);
            confilt = NaN(numel(constr), numel(cfilt));
            conmat = NaN(numel(constr), numel(x) + 1);
            % Penalty parameter for constraint in merit function (PARMU in Powell's code)
            cval = NaN(numel(x) + 1, 1);
            d = NaN(size(x));
            delbar = NaN;

            distsq = NaN(numel(x) + 1, 1);
            dnorm = NaN;
            ffilt = NaN(size(cfilt));
            fval = NaN(numel(x) + 1, 1);
            g = NaN(size(x));

            % Predicted reduction in constraint violation
            % Predicted reduction in objective Function
            % Predicted reduction in merit function
            % Reduction ratio: ACTREM/PREREM

            sim = NaN(numel(x), numel(x) + 1);

            xfilt = NaN(numel(x), numel(cfilt));
            % CPENMIN is the minimum of the penalty parameter CPEN for the L-infinity constraint violation in
            % the merit function. Note that CPENMIN = 0 in Powell's implementation, which allows CPEN to be 0.
            % Here, we take CPENMIN > 0 so that CPEN is always positive. This avoids the situation where PREREM
            % becomes 0 when PREREF = 0 = CPEN. It brings two advantages as follows.
            % 1. If the trust-region subproblem solver works correctly and the trust-region center is not
            % optimal for the subproblem, then PREREM > 0 is guaranteed. This is because, in theory, PREREC >= 0
            % and MAX(PREREC, PREREF) > 0 , and the definition of CPEN in GETCPEN ensures that PREREM > 0.
            % 2. There is no need to revise ACTREM and PREREM when CPEN = 0 and F = FVAL(N+1) as in lines
            % 312--314 of Powell's cobylb.f code. Powell's code revises ACTREM to CVAL(N + 1) - CSTRV and PREREM
            % to PREREC in this case, which is crucial for feasibility problems.
            cpenmin = eps(1.0);

            m_lcon = numel(bvec);
            m = numel(constr);
            n = numel(x);

            %====================%
            % Calculation starts %
            %====================%

            % Initialize SIM, SIMI, FVAL, CONMAT, and CVAL, together with the history, NF, and EVALUATED.
            % After the initialization, SIM(:, N+1) holds the vertex of the initial simplex with the smallest
            % function value (regardless of the constraint violation), and SIM(:, 1:N) holds the displacements
            % from the other vertices to SIM(:, N+1). FVAL, CONMAT, and CVAL hold the function values,
            % constraint values, and constraint violations on the vertices in the order corresponding to SIM.
            [nf, chist, conhist, conmat, cval, fhist, fval, sim, simi, xhist, evaluated, subinfo] = initialize_cobyla_obj.initxfc(calcfc, iprint, maxfun, amat, bvec, constr, ctol, f, ftarget, rhobeg, x, chist, conhist, conmat, cval, fhist, fval, sim, xhist, evaluated);

            % Report the current best value, and check if user asks for early termination.
            terminate = false;
            ipObj = inputParser();
            addParameter(ipObj, 'callback_fcn', struct());
            parse(ipObj, varargin{:});
            callback_fcn = ipObj.Results.callback_fcn;
            if ~ismember('callback_fcn', ipObj.UsingDefaults)
                terminate = callback_fcn(sim(:, n + 1), fval(n + 1), nf, 0, 'cstrv', cval(n + 1), 'nlconstr', conmat(m_lcon + 1:m, n + 1));
                if terminate
                    subinfo = 30;
                end
            end

            % Initialize the filter, including XFILT, FFILT, CONFILT, CFILT, and NFILT.
            % N.B.: The filter is used only when selecting which iterate to return. It does not interfere with
            % the iterations. COBYLA is NOT a filter method but a trust-region method based on an L-infinity
            % merit function. Powell's implementation does not use a filter to select the iterate, possibly
            % returning a suboptimal iterate.
            [nfilt, cfilt, confilt, ffilt, xfilt] = initialize_cobyla_obj.initfilt(conmat, ctol, cweight, cval, fval, sim, evaluated, cfilt, confilt, ffilt, xfilt);

            % Check whether to return due to abnormal cases that may occur during the initialization.
            if subinfo ~= 0
                info = subinfo;
                % Return the best calculated values of the variables.
                % N.B. SELECTX and FINDPOLE choose X by different standards. One cannot replace the other.
                kopt = selectx_obj.selectx(ffilt(1:nfilt), cfilt(1:nfilt), cweight, ctol);
                x = xfilt(:, kopt);
                f = ffilt(kopt);
                constr = confilt(:, kopt);
                cstrv = cfilt(kopt);
                % Arrange CHIST, CONHIST, FHIST, and XHIST so that they are in the chronological order.
                [xhist, fhist, chist, conhist] = history_obj.rangehist(nf, xhist, fhist, 'chist', chist, 'conhist', conhist);
                % Print a return message according to IPRINT.
                message_obj.retmsg(solver, info, iprint, nf, f, x, 'cstrv', cstrv, 'constr', constr);

                return
            end

            % Set some more initial values.
            % We must initialize ACTREM and PREREM. Otherwise, when SHORTD = TRUE, compilers may raise a
            % run-time error that they are undefined. But their values will not be used: when SHORTD = FALSE,
            % they will be overwritten; when SHORTD = TRUE, the values are used only in BAD_TRSTEP, which is
            % TRUE regardless of ACTREM or PREREM. Similar for PREREC, PREREF, PREREM, RATIO, and JDROP_TR.
            % No need to initialize SHORTD unless MAXTR < 1, but some compilers may complain if we do not do it.
            % Our initialization of CPEN differs from Powell's in two ways. First, we use the ratio defined in
            % (13) of Powell's COBYLA paper to initialize CPEN. Second, we impose CPEN >= CPENMIN > 0. Powell's
            % code simply initializes CPEN to 0.
            rho = rhobeg;
            delta = rhobeg;
            cpen = max(cpenmin, min(1000.0, obj.fcratio(conmat, fval))); % Powell's code: CPEN = ZERO
            prerec = -realmax;
            preref = -realmax;
            prerem = -realmax;
            actrem = -realmax;
            shortd = false;
            trfail = false;
            ratio = -1.0;
            jdrop_tr = 0;
            jdrop_geo = 0;

            % If DELTA <= GAMMA3*RHO after an update, we set DELTA to RHO. GAMMA3 must be less than GAMMA2. The
            % reason is as follows. Imagine a very successful step with DENORM = the un-updated DELTA = RHO.
            % Then TRRAD will update DELTA to GAMMA2*RHO. If GAMMA3 >= GAMMA2, then DELTA will be reset to RHO,
            % which is not reasonable as D is very successful. See paragraph two of Sec. 5.2.5 in
            % T. M. Ragonneau's thesis: "Model-Based Derivative-Free Optimization Methods and Software".
            % According to test on 20230613, for COBYLA, this Powellful updating scheme of DELTA works slightly
            % better than setting directly DELTA = MAX(NEW_DELTA, RHO).
            gamma3 = max(1.0, min(0.75 * gamma2, 1.5));

            % MAXTR is the maximal number of trust-region iterations. Here, we set it to HUGE(MAXTR) - 1 so that
            % the algorithm will not terminate due to MAXTR. However, this may not be allowed in other languages
            % such as MATLAB. In that case, we can set MAXTR to 10*MAXFUN, which is unlikely to reach because
            % each trust-region iteration takes 1 or 2 function evaluations unless the trust-region step is short
            % or fails to reduce the trust-region model but the geometry step is not invoked.
            % N.B.: Do NOT set MAXTR to HUGE(MAXTR), as it may cause overflow and infinite cycling in the DO
            % loop. See
            % https://fortran-lang.discourse.group/t/loop-variable-reaching-integer-huge-causes-infinite-loop
            % https://fortran-lang.discourse.group/t/loops-dont-behave-like-they-should
            maxtr = intmax('int32') - 1; %%MATLAB: maxtr = 10 * maxfun;
            info = 20;

            % Begin the iterative procedure.
            % After solving a trust-region subproblem, we use three boolean variables to control the workflow.
            % SHORTD - Is the trust-region trial step too short to invoke a function evaluation?
            % IMPROVE_GEO - Will we improve the model after the trust-region iteration? If yes, a geometry step
            % will be taken, corresponding to the "Branch (Delta)" in the COBYLA paper.
            % REDUCE_RHO - Will we reduce rho after the trust-region iteration?
            % COBYLA never sets IMPROVE_GEO and REDUCE_RHO to TRUE simultaneously.
            for tr = 1:maxtr
                % Increase the penalty parameter CPEN, if needed, so that PREREM = PREREF + CPEN * PREREC > 0.
                % This is the first (out of two) update of CPEN, where CPEN increases or remains the same.
                % N.B.: CPEN and the merit function PHI = FVAL + CPEN*CVAL are used at three places only.
                % 1. In FINDPOLE/UPDATEPOLE, deciding the optimal vertex of the current simplex.
                % 2. After the trust-region trial step, calculating the reduction radio.
                % 3. In GEOSTEP, deciding the direction of the geometry step.
                % They do not appear explicitly in the trust-region subproblem, though the trust-region center
                % (i.e., the current optimal vertex) is defined by them.
                cpen = obj.getcpen(amat, bvec, conmat, cpen, cval, delta, fval, sim, simi);

                % Switch the best vertex of the current simplex to SIM(:, N + 1).
                [conmat, cval, fval, sim, simi, subinfo] = update_cobyla_obj.updatepole(cpen, conmat, cval, fval, sim, simi);
                % Check whether to exit due to damaging rounding in UPDATEPOLE.
                if subinfo == 7
                    info = subinfo;
                    break % Better action to take? Geometry step, or simply continue?

                end

                % Does the interpolation set have adequate geometry? It affects IMPROVE_GEO and REDUCE_RHO.
                adequate_geo = all(sum(sim(:, 1:n) .^ 2, 1) <= 4.0 * delta ^ 2, 'all');

                % Calculate the linear approximations to the objective and constraint functions.
                % N.B.: TRSTLP accesses A mostly by columns, so it is more reasonable to save A instead of A^T.
                % Zaikun 2023108: According to a test on 2023108, calculating G and A(:, M_LCON+1:M) by solving
                % the linear systems SIM^T*G = FVAL(1:N)-FVAL(N+1) and SIM^T*A = CONMAT(:, 1:N)-CONMAT(:, N+1)
                % does not seem to improve or worsen the performance of COBYLA in terms of the number of function
                % evaluations. The system was solved by SOLVE in LINALG_MOD based on a QR factorization of SIM
                % (not necessarily a good algorithm). No preconditioning or scaling was used.
                g(:) = simi.' * (fval(1:n) - fval(n + 1));
                A(:, 1:m_lcon) = amat;
                A(:, m_lcon + 1:m) = ((conmat(m_lcon + 1:m, 1:n) - conmat(m_lcon + 1:m, n + 1)) * simi).';
                %%MATLAB: A(:, m_lcon+1:m) = simi'*(conmat(m_lcon+1:m, 1:n) - conmat(m_lcon+1:m, n+1))' % Implicit expansion for subtraction

                % Calculate the trust-region trial step D. Note that D does NOT depend on CPEN.
                d(:) = trustregion_cobyla_obj.trstlp(A, -conmat(:, n + 1), delta, g);
                dnorm = min(delta, norm(d));

                % Is the trust-region trial step short? Note that we compare DNORM with RHO, not DELTA.
                % Powell's code essentially defines SHORTD by SHORTD = (DNORM < HALF * RHO). In our tests,
                % TENTH seems to work better than HALF or QUART, especially for linearly constrained problems.
                % Note that LINCOA has a slightly more sophisticated way of defining SHORTD, taking into account
                % whether D causes a change to the active set. Should we try the same here?
                shortd = (dnorm <= 0.1 * rho); % `<=` works better than `<` in case of underflow.

                % Predict the change to F (PREREF) and to the constraint violation (PREREC) due to D.
                % We have the following in precise arithmetic. They may fail to hold due to rounding errors.
                % 1. PREREC is the reduction of the L-infinity violation of the linearized constraints achieved
                % by D. It is nonnegative in theory; it is 0 iff CONMAT(1:M, N+1) <= 0, namely the trust-region
                % center satisfies the constraints.
                % 2. PREREF may be negative or 0, but it should be positive when PREREC = 0 and SHORTD is FALSE.
                % 3. Due to 2, in theory, MAXIMUM([PREREC, PREREF]) > 0 if SHORTD is FALSE.
                preref = -sum(d .* g, 'all'); % Can be negative.
                prerec = cval(n + 1) - max([0.0; conmat(:, n + 1) + A.' * d], [], 'all');

                % Evaluate PREREM, which is the predicted reduction in the merit function.
                % In theory, PREREM >= 0 and it is 0 iff CPEN = 0 = PREREF. This may not be true numerically.
                prerem = preref + cpen * prerec;
                trfail = (~(prerem > 1.0e-6 * min(cpen, 1.0) * rho)); % PREREM is tiny/negative or NaN.

                if shortd || trfail
                    % Reduce DELTA if D is short or D fails to render PREREM > 0. The latter can happen due to
                    % rounding errors. This seems important for performance.
                    delta = 0.1 * delta;
                    if delta <= gamma3 * rho
                        delta = rho; % Set DELTA to RHO when it is close to or below.

                    end
                else
                    % Calculate the next value of the objective and constraint functions.
                    % If X is close to one of the points in the interpolation set, then we do not evaluate the
                    % objective and constraints at X, assuming them to have the values at the closest point.
                    % N.B.: If this happens, do NOT include X into the filter, as F and CONSTR are inaccurate.
                    x = sim(:, n + 1) + d;
                    distsq(n + 1) = sum((x - sim(:, n + 1)) .^ 2, 'all');
                    distsq(1:n) = arrayfun(@(j) sum((x - (sim(:, n + 1) + sim(:, j))) .^ 2, 'all'), (1:n)'); % Implied do-loop
                    %%MATLAB: distsq(1:n) = sum((x - (sim(:,1:n) + sim(:, n+1)))**2, 1)  % Implicit expansion
                    [~, j] = min(distsq);
                    if distsq(j) <= (1.0e-4 * rhoend) ^ 2
                        f = fval(j);
                        constr = conmat(:, j);
                        cstrv = cval(j);
                    else
                        % Evaluate the objective and constraints at X, taking care of possible Inf/NaN values.
                        constr(1:m_lcon) = evaluate_obj.moderatec(amat.' * x - bvec); % Linear constraints
                        [f, constr_slice] = evaluate_obj.evaluatefc(calcfc, x, constr(m_lcon + 1:m)); constr(m_lcon + 1:m) = constr_slice; % Nonlinear constraints
                        % Note that EVALUATE moderates the nonlinear constraint values. Thus we also moderate the
                        % linear constraint values here to make CSTRV consistent.
                        cstrv = max([0.0; constr], [], 'all');
                        nf = nf + 1;
                        % Save X, F, CONSTR, CSTRV into the history.
                        [xhist, fhist, chist, conhist] = history_obj.savehist(nf, x, xhist, f, fhist, 'cstrv', cstrv, 'chist', chist, 'constr', constr, 'conhist', conhist);
                        % Save X, F, CONSTR, CSTRV into the filter.
                        [nfilt, cfilt, ffilt, xfilt, confilt] = selectx_obj.savefilt(cstrv, ctol, cweight, f, x, nfilt, cfilt, ffilt, xfilt, 'constr', constr, 'confilt', confilt);
                    end

                    % Print a message about the function/constraint evaluation according to IPRINT.
                    message_obj.fmsg(solver, "Trust region", iprint, nf, delta, f, x, 'cstrv', cstrv, 'constr', constr);

                    % Evaluate ACTREM, which is the actual reduction in the merit function.
                    actrem = (fval(n + 1) + cpen * cval(n + 1)) - (f + cpen * cstrv);

                    % Calculate the reduction ratio by REDRAT, which handles Inf/NaN carefully.
                    ratio = ratio_obj.redrat(actrem, prerem, eta1);

                    % Update DELTA. After this, DELTA < DNORM may hold.
                    % N.B.: 1. Powell's code uses RHO as the trust-region radius and updates it as follows.
                    % Reduce RHO to GAMMA1*RHO if ADEQUATE_GEO is TRUE and either SHORTD is TRUE or RATIO < ETA1,
                    % and then revise RHO to RHOEND if its new value is not more than GAMMA3*RHOEND; RHO remains
                    % unchanged in all other cases; in particular, RHO is never increased.
                    % 2. Our implementation uses DELTA as the trust-region radius, while using RHO as a lower
                    % bound for DELTA. DELTA is updated in a way that is typical for trust-region methods, and
                    % it is revised to RHO if its new value is not more than GAMMA3*RHO. RHO reflects the current
                    % resolution of the algorithm; its update is essentially the same as the update of RHO in
                    % Powell's code (see the definition of REDUCE_RHO below). Our implementation aligns with
                    % UOBYQA/NEWUOA/BOBYQA/LINCOA and improves the performance of COBYLA.
                    % 3. The same as Powell's code, we do not reduce RHO unless ADEQUATE_GEO is TRUE. This is
                    % also how Powell updated RHO in UOBYQA/NEWUOA/BOBYQA/LINCOA. What about we also use
                    % ADEQUATE_GEO == TRUE as a prerequisite for reducing DELTA? The argument would be that the
                    % bad (small) value of RATIO may be because of a bad geometry (and hence a bad model) rather
                    % than an improperly large DELTA, and it might be good to try improving the geometry first
                    % without reducing DELTA. However, according to a test on 20230206, it does not improve the
                    % performance if we skip the update of DELTA when ADEQUATE_GEO is FALSE and RATIO < 0.1.
                    % Therefore, we choose to update DELTA without checking ADEQUATE_GEO.
                    delta = trustregion_cobyla_obj.trrad(delta, dnorm, eta1, eta2, gamma1, gamma2, ratio);
                    if delta <= gamma3 * rho
                        delta = rho; % Set DELTA to RHO when it is close to or below.

                    end

                    % Is the newly generated X better than current best point?
                    ximproved = (actrem > 0); % If ACTREM is NaN, then XIMPROVED should & will be FALSE.

                    % Set JDROP_TR to the index of the vertex to be replaced with X. JDROP_TR = 0 means there
                    % is no good point to replace, and X will not be included into the simplex; in this case,
                    % the geometry of the simplex likely needs improvement, which will be handled below.
                    jdrop_tr = geometry_cobyla_obj.setdrop_tr(ximproved, d, delta, rho, sim, simi);

                    % Update SIM, SIMI, FVAL, CONMAT, and CVAL so that SIM(:, JDROP_TR) is replaced with D.
                    % UPDATEXFC does nothing if JDROP_TR == 0, as the algorithm decides to discard X.
                    [conmat, cval, fval, sim, simi, subinfo] = update_cobyla_obj.updatexfc(jdrop_tr, constr, cpen, cstrv, d, f, conmat, cval, fval, sim, simi);
                    % Check whether to exit due to damaging rounding in UPDATEXFC.
                    if subinfo == 7
                        info = subinfo;
                        break % Better action to take? Geometry step, or a RESCUE as in BOBYQA?

                    end

                    % Check whether to exit due to MAXFUN, FTARGET, etc.
                    subinfo = checkexit_obj.checkexit_con(maxfun, nf, cstrv, ctol, f, ftarget, x);
                    if subinfo ~= 0
                        info = subinfo;
                        break
                    end
                end % End of IF (SHORTD .OR. TRFAIL). The normal trust-region calculation ends.


                %----------------------------------------------------------------------------------------------%
                % Before the next trust-region iteration, we possibly improve the geometry of simplex or reduce
                % RHO according to IMPROVE_GEO and REDUCE_RHO. Now we decide these indicators.
                % N.B.: We must ensure that the algorithm does not set IMPROVE_GEO = TRUE at infinitely many
                % consecutive iterations without moving SIM(:, N+1) or reducing RHO. Otherwise, the algorithm
                % will get stuck in repetitive invocations of GEOSTEP. This is ensured by the following facts.
                % 1. If an iteration sets IMPROVE_GEO = TRUE, it must also reduce DELTA or set DELTA to RHO.
                % 2. If SIM(:, N+1) and RHO remains unchanged, then ADEQUATE_GEO will become TRUE after at
                % most N invocations of GEOSTEP.

                % BAD_TRSTEP: Is the last trust-region step bad?
                bad_trstep = (shortd || trfail || ratio <= 0 || jdrop_tr == 0);
                % IMPROVE_GEO: Should we take a geometry step to improve the geometry of the interpolation set?
                improve_geo = (bad_trstep && ~adequate_geo);
                % REDUCE_RHO: Should we enhance the resolution by reducing RHO?
                reduce_rho = (bad_trstep && adequate_geo && max(delta, dnorm) <= rho);

                % COBYLA never sets IMPROVE_GEO and REDUCE_RHO to TRUE simultaneously.
                % %call assert(.not. (improve_geo .and. reduce_rho), 'IMPROVE_GEO or REDUCE_RHO are not both TRUE', srname)

                % If SHORTD or TRFAIL is TRUE, then either IMPROVE_GEO or REDUCE_RHO is TRUE unless ADEQUATE_GEO
                % is TRUE and MAX(DELTA, DNORM) > RHO.
                % %call assert((.not. (shortd .or. trfail)) .or. (improve_geo .or. reduce_rho .or. &
                % %    & (adequate_geo .and. max(delta, dnorm) > rho)), 'If SHORTD or TRFAIL is TRUE, then &
                % %    & either IMPROVE_GEO or REDUCE_RHO is TRUE unless ADEQUATE_GEO is TRUE and MAX(DELTA, DNORM) > RHO', srname)
                %----------------------------------------------------------------------------------------------%

                % Comments on BAD_TRSTEP:
                % 1. Powell's definition of BAD_TRSTEP is as follows. The one used above seems to work better,
                % especially for linearly constrained problems due to the factor TENTH (= ETA1).
                % %bad_trstep = (shortd .or. actrem <= 0 .or. actrem < TENTH * prerem .or. jdrop_tr == 0)
                % Besides, Powell did not check PREREM > 0 in BAD_TRSTEP, which is reasonable to do but has
                % little impact upon the performance.
                % 2. NEWUOA/BOBYQA/LINCOA would define BAD_TRSTEP, IMPROVE_GEO, and REDUCE_RHO as follows. Two
                % different thresholds are used in BAD_TRSTEP. It outperforms Powell's version.
                % %bad_trstep = (shortd .or. trfail .or. ratio <= eta1 .or. jdrop_tr == 0)
                % %improve_geo = bad_trstep .and. .not. adequate_geo
                % %bad_trstep = (shortd .or. trfail .or. ratio <= 0 .or. jdrop_tr == 0)
                % %reduce_rho = bad_trstep .and. adequate_geo .and. max(delta, dnorm) <= rho
                % 3. Theoretically, JDROP_TR > 0 when ACTREM > 0 (guaranteed by RATIO > 0). However, in Powell's
                % implementation, JDROP_TR may be 0 even RATIO > 0 due to NaN. The modernized code has rectified
                % this in the function SETDROP_TR. After this rectification, we can indeed simplify the
                % definition of BAD_TRSTEP by removing the condition JDROP_TR == 0. We retain it for robustness.

                % Comments on REDUCE_RHO:
                % When SHORTD is TRUE, UOBYQA/NEWUOA/BOBYQA/LINCOA all set REDUCE_RHO to TRUE if the recent
                % models are sufficiently accurate according to certain criteria. See the paragraph around (37)
                % in the UOBYQA paper and the discussions about Box 14 in the NEWUOA paper. This strategy is
                % crucial for the performance of the solvers. However, as of 20221111, we have not managed to
                % make it work in COBYLA. As in NEWUOA, we recorded the errors of the recent models, and set
                % REDUCE_RHO to true if they are small (e.g., ALL(ABS(MODERR_REC) <= 0.1 * MAXVAL(ABS(A))*RHO) or
                % ALL(ABS(MODERR_REC) <= RHO**2)) when SHORTD is TRUE. It made little impact on the performance.


                % Since COBYLA never sets IMPROVE_GEO and REDUCE_RHO to TRUE simultaneously, the following
                % two blocks are exchangeable: IF (IMPROVE_GEO) ... END IF and IF (REDUCE_RHO) ... END IF.

                % Improve the geometry of the simplex by removing a point and adding a new one.
                % If the current interpolation set has adequate geometry, then we skip the geometry step.
                % The code has a small difference from Powell's original code here: If the current geometry
                % is adequate, then we will continue with a new trust-region iteration; however, at the
                % beginning of the iteration, CPEN may be updated, which may alter the pole point SIM(:, N+1)
                % by UPDATEPOLE; the quality of the interpolation point depends on SIM(:, N + 1), meaning
                % that the same interpolation set may have good or bad geometry with respect to different
                % "poles"; if the geometry turns out bad with the new pole, the original COBYLA code will
                % take a geometry step, yet our code here will NOT do it but continue to take a trust-region
                % step. The argument is this: even if the geometry step is not skipped in the first place, the
                % geometry may turn out bad again after the pole is altered due to an update to CPEN; should
                % we take another geometry step in that case? If no, why should we do it here? Indeed, this
                % distinction makes no practical difference for CUTEst problems with at most 100 variables
                % and 5000 constraints, while the algorithm framework is simplified.
                if improve_geo && ~all(sum(sim(:, 1:n) .^ 2, 1) <= 4.0 * delta ^ 2, 'all')
                    % Before the geometry step, UPDATEPOLE has been called either implicitly by UPDATEXFC or
                    % explicitly after CPEN is updated, so that SIM(:, N + 1) is the optimal vertex.

                    % Decide a vertex to drop from the simplex. It will be replaced with SIM(:, N + 1) + D to
                    % improve the geometry of the simplex.
                    % N.B.: 1. COBYLA never sets JDROP_GEO = N + 1.
                    % 2. The following JDROP_GEO comes from UOBYQA/NEWUOA/BOBYQA/LINCOA.
                    % 3. In Powell's original algorithm, the geometry of the simplex is considered acceptable
                    % iff the distance between any vertex and the pole is at most 2.1*DELTA, and the distance
                    % between any vertex and the opposite face of the simplex is at least 0.25*DELTA, as
                    % specified in (14) of the COBYLA paper. Correspondingly, JDROP_GEO is set to the index of
                    % the vertex with the largest distance to the pole provided that the distance is larger than
                    % 2.1*DELTA, or the vertex with the smallest distance to the opposite face of the simplex,
                    % in which case the distance must be less than 0.25*DELTA, as the current simplex does not
                    % have acceptable geometry (see (15)--(16) of the COBYLA paper). Once JDROP_GEO is set, the
                    % algorithm replaces SIM(:, JDROP_GEO) with D specified in (17) of the COBYLA paper, which
                    % is orthogonal to the face opposite to SIM(:, JDROP_GEO) and has a length of 0.5*DELTA,
                    % intending to improve the geometry of the simplex as per (14).
                    % 4. Powell's geometry-improving procedure outlined above has an intrinsic flaw: it may lead
                    % to infinite cycling, as was observed in a test on 20240320. In this test, the geometry-
                    % improving point introduced in the previous iteration was replaced with the trust-region
                    % trial point in the current iteration, which was then replaced with the same geometry-
                    % improving point in the next iteration, and so on. In this process, the simplex alternated
                    % between two configurations, neither of which had acceptable geometry. Thus RHO was never
                    % reduced, leading to infinite cycling. (N.B.: Our implementation uses DELTA as the trust
                    % region radius, with RHO being its lower bound. When the infinite cycling occurred in this
                    % test, DELTA = RHO and it could not be reduced due to the requirement that DELTA >= RHO.)
                    jdrop_geo = fortran.maxloc(sum(sim(:, 1:n) .^ 2, 1), 'dim', 1);

                    % Calculate the geometry step D.
                    delbar = 0.5 * delta;
                    d(:) = geometry_cobyla_obj.geostep(jdrop_geo, amat, bvec, conmat, cpen, cval, delbar, fval, simi);

                    % Calculate the next value of the objective and constraint functions.
                    % If X is close to one of the points in the interpolation set, then we do not evaluate the
                    % objective and constraints at X, assuming them to have the values at the closest point.
                    % N.B.:
                    % 1. If this happens, do NOT include X into the filter, as F and CONSTR are inaccurate.
                    % 2. In precise arithmetic, the geometry improving step ensures that the distance between X
                    % and any interpolation point is at least DELBAR, yet X may be close to them due to
                    % rounding. In an experiment with single precision on 20240317, X = SIM(:, N+1) occurred.
                    x = sim(:, n + 1) + d;
                    distsq(n + 1) = sum((x - sim(:, n + 1)) .^ 2, 'all');
                    distsq(1:n) = arrayfun(@(j) sum((x - (sim(:, n + 1) + sim(:, j))) .^ 2, 'all'), (1:n)'); % Implied do-loop
                    %%MATLAB: distsq(1:n) = sum((x - (sim(:,1:n) + sim(:, n+1)))**2, 1)  % Implicit expansion
                    [~, j] = min(distsq);
                    if distsq(j) <= (1.0e-4 * rhoend) ^ 2
                        f = fval(j);
                        constr = conmat(:, j);
                        cstrv = cval(j);
                    else
                        % Evaluate the objective and constraints at X, taking care of possible Inf/NaN values.
                        constr(1:m_lcon) = evaluate_obj.moderatec(amat.' * x - bvec); % Linear constraints
                        [f, constr_slice] = evaluate_obj.evaluatefc(calcfc, x, constr(m_lcon + 1:m)); constr(m_lcon + 1:m) = constr_slice; % Nonlinear constraints
                        % Note that EVALUATE moderates the nonlinear constraint values. Thus we also moderate the
                        % linear constraint values here to make CSTRV consistent.
                        cstrv = max([0.0; constr], [], 'all');
                        nf = nf + 1;
                        % Save X, F, CONSTR, CSTRV into the history.
                        [xhist, fhist, chist, conhist] = history_obj.savehist(nf, x, xhist, f, fhist, 'cstrv', cstrv, 'chist', chist, 'constr', constr, 'conhist', conhist);
                        % Save X, F, CONSTR, CSTRV into the filter.
                        [nfilt, cfilt, ffilt, xfilt, confilt] = selectx_obj.savefilt(cstrv, ctol, cweight, f, x, nfilt, cfilt, ffilt, xfilt, 'constr', constr, 'confilt', confilt);
                    end

                    % Print a message about the function/constraint evaluation according to IPRINT.
                    message_obj.fmsg(solver, "Geometry", iprint, nf, delta, f, x, 'cstrv', cstrv, 'constr', constr);
                    % Update SIM, SIMI, FVAL, CONMAT, and CVAL so that SIM(:, JDROP_GEO) is replaced with D.
                    [conmat, cval, fval, sim, simi, subinfo] = update_cobyla_obj.updatexfc(jdrop_geo, constr, cpen, cstrv, d, f, conmat, cval, fval, sim, simi);
                    % Check whether to exit due to damaging rounding in UPDATEXFC.
                    if subinfo == 7
                        info = subinfo;
                        break % Better action to take? Geometry step, or simply continue?

                    end

                    % Check whether to exit due to MAXFUN, FTARGET, etc.
                    subinfo = checkexit_obj.checkexit_con(maxfun, nf, cstrv, ctol, f, ftarget, x);
                    if subinfo ~= 0
                        info = subinfo;
                        break
                    end
                end % End of IF (IMPROVE_GEO). The procedure of improving geometry ends.

                % The calculations with the current RHO are complete. Enhance the resolution of the algorithm
                % by reducing RHO; update DELTA and CPEN at the same time.
                if reduce_rho
                    if rho <= rhoend
                        info = 0;
                        break
                    end
                    delta = max(0.5 * rho, redrho_obj.redrho(rho, rhoend));
                    rho = redrho_obj.redrho(rho, rhoend);
                    % The second (out of two) update of CPEN, where CPEN decreases or remains the same.
                    % Powell's code: CPEN = MIN(CPEN, FCRATIO(FVAL, CONMAT)), which may set CPEN to 0.
                    cpen = max(cpenmin, min(cpen, obj.fcratio(conmat, fval)));
                    % Print a message about the reduction of RHO according to IPRINT.
                    message_obj.rhomsg(solver, iprint, nf, delta, fval(n + 1), rho, sim(:, n + 1), 'cstrv', cval(n + 1), 'constr', conmat(:, n + 1), 'cpen', cpen);
                    % Switch the best vertex of the current simplex to SIM(:, N + 1).
                    [conmat, cval, fval, sim, simi, subinfo] = update_cobyla_obj.updatepole(cpen, conmat, cval, fval, sim, simi);
                    % Check whether to exit due to damaging rounding in UPDATEPOLE.
                    if subinfo == 7
                        info = subinfo;
                        break % Better action to take? Geometry step, or simply continue?

                    end
                end % End of IF (REDUCE_RHO). The procedure of reducing RHO ends.

                % Report the current best value, and check if user asks for early termination.
                if ~ismember('callback_fcn', ipObj.UsingDefaults)
                    terminate = callback_fcn(sim(:, n + 1), fval(n + 1), nf, tr, 'cstrv', cval(n + 1), 'nlconstr', conmat(m_lcon + 1:m, n + 1));
                    if terminate
                        info = 30;
                        break
                    end
                end

            end % End of DO TR = 1, MAXTR. The iterative procedure ends.

            % Return from the calculation, after trying the last trust-region step if it has not been tried yet.
            % Ensure that D has not been updated after SHORTD == TRUE occurred, or the code below is incorrect.
            x = sim(:, n + 1) + d;
            if info == 0 && shortd && norm(x - sim(:, n + 1)) > 1.0e-3 * rhoend && nf < maxfun
                constr(1:m_lcon) = evaluate_obj.moderatec(amat.' * x - bvec); % Linear constraints
                [f, constr_slice] = evaluate_obj.evaluatefc(calcfc, x, constr(m_lcon + 1:m)); constr(m_lcon + 1:m) = constr_slice; % Nonlinear constraints
                % Note that EVALUATE moderates the nonlinear constraint values. Thus we also moderate the linear
                % constraint values here to make CSTRV consistent.
                cstrv = max([0.0; constr], [], 'all');
                nf = nf + 1;
                % Save X, F, CONSTR, CSTRV into the history.
                [xhist, fhist, chist, conhist] = history_obj.savehist(nf, x, xhist, f, fhist, 'cstrv', cstrv, 'chist', chist, 'constr', constr, 'conhist', conhist);
                % Save X, F, CONSTR, CSTRV into the filter.
                [nfilt, cfilt, ffilt, xfilt, confilt] = selectx_obj.savefilt(cstrv, ctol, cweight, f, x, nfilt, cfilt, ffilt, xfilt, 'constr', constr, 'confilt', confilt);
                % Print a message about the function/constraint evaluation according to IPRINT.
                % Zaikun 20230512: DELTA has been updated. RHO is only indicative here. TO BE IMPROVED.
                message_obj.fmsg(solver, "Trust region", iprint, nf, rho, f, x, 'cstrv', cstrv, 'constr', constr);
            end

            % Return the best calculated values of the variables.
            % N.B. SELECTX and FINDPOLE choose X by different standards. One cannot replace the other.
            kopt = selectx_obj.selectx(ffilt(1:nfilt), cfilt(1:nfilt), max(cpen, cweight), ctol);
            x = xfilt(:, kopt);
            f = ffilt(kopt);
            constr = confilt(:, kopt);
            cstrv = cfilt(kopt);

            % Arrange CHIST, CONHIST, FHIST, and XHIST so that they are in the chronological order.
            [xhist, fhist, chist, conhist] = history_obj.rangehist(nf, xhist, fhist, 'chist', chist, 'conhist', conhist);

            % Print a return message according to IPRINT.
            message_obj.retmsg(solver, info, iprint, nf, f, x, 'cstrv', cstrv, 'constr', constr);
            %====================%
            %  Calculation ends  %
            %====================%


        end
        function cpen = getcpen(~, amat, bvec, conmat_in, cpen_in, cval_in, delta, fval_in, sim_in, simi_in)
            %--------------------------------------------------------------------------------------------------%
            % This function gets the penalty parameter CPEN so that PREREM = PREREF + CPEN * PREREC > 0.
            % See the discussions around equation (9) of the COBYLA paper.
            %--------------------------------------------------------------------------------------------------%


            % Solver-specific modules
            trustregion_cobyla_obj = prima_mat.cobyla.trustregion_cobyla_mod();
            update_cobyla_obj = prima_mat.cobyla.update_cobyla_mod();

            A = NaN(size(sim_in, 1), size(conmat_in, 1));
            conmat = NaN(size(conmat_in, 1), size(conmat_in, 2));

            d = NaN(size(sim_in, 1), 1);

            g = NaN(size(sim_in, 1), 1);

            sim = NaN(size(sim_in, 1), size(sim_in, 2));

            m_lcon = numel(bvec);
            m = size(conmat, 1);
            n = size(sim, 1);

            %====================%
            % Calculation starts %
            %====================%

            % Copy the inputs.
            conmat = conmat_in;
            cpen = cpen_in;
            cval = cval_in;
            fval = fval_in;
            sim = sim_in;
            simi = simi_in;

            % Initialize INFO, PREREF, and PREREC, which are needed in the postconditions.

            preref = 0.0;
            prerec = 0.0;

            % Increase CPEN if necessary to ensure PREREM > 0. Branch back for the next loop if this change
            % alters the optimal vertex of the current simplex. Note the following.
            % 1. In each loop, CPEN is changed only if PREREC > 0 > PREREF, in which case PREREM is guaranteed
            % positive after the update. Note that PREREC >= 0 and MAX(PREREC, PREREF) > 0 in theory. If this
            % holds numerically as well, then CPEN is not changed only if PREREC = 0 or PREREF >= 0, in which
            % case PREREM is currently positive, explaining why CPEN needs no update.
            % 2. Even without an upper bound for the loop counter, the loop can occur at most N+1 times. This is
            % because the update of CPEN does not decrease CPEN, and hence it can make vertex J (J <= N) become
            % the new optimal vertex only if CVAL(J) is less than CVAL(N+1), which can happen at most N times.
            % See the paragraph below (9) in the COBYLA paper. After the "correct" optimal vertex is found,
            % one more loop is needed to calculate CPEN, and hence the loop can occur at most N+1 times.
            for iter = 1:n + 1
                % Switch the best vertex of the current simplex to SIM(:, N + 1).
                [conmat, cval, fval, sim, simi, info] = update_cobyla_obj.updatepole(cpen, conmat, cval, fval, sim, simi);
                % Check whether to exit due to damaging rounding in UPDATEPOLE.
                if info == 7
                    break
                end

                % Calculate the linear approximations to the objective and constraint functions.
                g(:) = simi.' * (fval(1:n) - fval(n + 1));
                A(:, 1:m_lcon) = amat;
                A(:, m_lcon + 1:m) = ((conmat(m_lcon + 1:m, 1:n) - conmat(m_lcon + 1:m, n + 1)) * simi).';
                %%MATLAB: A(:, m_lcon+1:m) = simi'*(conmat(m_lcon+1:m, 1:n) - conmat(m_lcon+1:m, n+1))' % Implicit expansion for subtraction

                % Calculate the trust-region trial step D. Note that D does NOT depend on CPEN.
                d(:) = trustregion_cobyla_obj.trstlp(A, -conmat(:, n + 1), delta, g);

                % Predict the change to F (PREREF) and to the constraint violation (PREREC) due to D.
                preref = -sum(d .* g, 'all'); % Can be negative.
                prerec = cval(n + 1) - max([0.0; conmat(:, n + 1) + A.' * d], [], 'all');

                if ~(prerec > 0 && preref < 0)
                    % PREREC <= 0 or PREREF >= 0 or either is NaN.
                    break
                end

                % Powell's code defines BARMU = -PREREF / PREREC, and CPEN is increased to 2*BARMU if and
                % only if it is currently less than 1.5*BARMU, a very "Powellful" scheme. In our implementation,
                % however, we set CPEN directly to the maximum between its current value and 2*BARMU while
                % handling possible overflow. This simplifies the scheme without worsening the performance.
                cpen = max(cpen, min(-2.0 * (preref / prerec), realmax));

                if update_cobyla_obj.findpole(cpen, cval, fval) == n + 1
                    break
                end
            end

            %====================%
            %  Calculation ends  %
            %====================%


        end
        function r = fcratio(~, conmat, fval)
            %--------------------------------------------------------------------------------------------------%
            % This function calculates the ratio between the "typical change" of F and that of CONSTR.
            % See equations (12)--(13) in Section 3 of the COBYLA paper for the definition of the ratio.
            %--------------------------------------------------------------------------------------------------%


            %====================%
            % Calculation starts %
            %====================%

            % N.B.: In the original version of COBYLA, Powell proposed the ratio for constraints in the form of
            % CONSTR(X) >= 0, but the constraints we consider here are CONSTR(X) <= 0. Hence we need to change
            % the sign of the constraints before defining CMIN and CMAX.
            cmin = min(-conmat, [], 2);
            cmax = max(-conmat, [], 2);
            fmin = min(fval, [], 'all');
            fmax = max(fval, [], 'all');
            r = 0.0;
            if any(cmin < 0.5 * cmax, 'all') && fmin < fmax
                denom = min(fortran.merge('tsource', max(cmax, 0.0) - cmin, 'fsource', realmax, 'mask', (cmin < 0.5 * cmax)), [], 'all');
                % Powell mentioned the following alternative in Section 4 of his COBYLA paper. According to a
                % test on 20230610, it does not make much difference to the performance.
                % %denom = maxval(max(cmax, ZERO) - cmin, mask=(cmin < HALF * cmax))
                r = (fmax - fmin) / denom;
            end

            %====================%
            %  Calculation ends  %
            %====================%


        end

    end
end