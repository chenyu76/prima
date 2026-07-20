% TODO:
% 1. Improve RESCUE so that it accepts an [XNEW, FNEW] that is not interpolated yet, or even accepts
% [XRESERVE, FRESERVE], which contains points that have been evaluated.
%
classdef bobyqb_mod
    %--------------------------------------------------------------------------------------------------%
    % This module performs the major calculations of BOBYQA.
    %
    % Coded by Zaikun ZHANG (www.zhangzk.net) based on Powell's code and the BOBYQA paper.
    %
    % N.B. (Zaikun 20230312): In Powell's code, the strategy concerning RESCUE is a bit complex.
    %
    % 1. Suppose that a trust-region step D is calculated. Powell's code sets KNEW_TR before evaluating
    % F at the trial point XOPT+D, assuming that the value of F at this point is not better than the
    % current FOPT. With this KNEW_TR, the denominator of the update is calculated. If this denominator
    % is sufficiently large, then evaluate F at XOPT+D, recalculate KNEW_TR if the function value turns
    % out better than FOPT, and perform the update to include XOPT+D in the interpolation. If the
    % denominator is not sufficiently large, then RESCUE is called, and another trust-region step is
    % taken immediately after, discarding the previously calculated trust-region step D.
    %
    % 2. Suppose that a geometry step D is calculated. Then KNEW_GEO must have been set before. Powell's
    % code then calculates the denominator of the update. If the denominator is sufficiently large, then
    % evaluate F at XOPT+D, and perform the update. If the denominator is not sufficiently large, then
    % RESCUE is called; if RESCUE does not evaluate F at any new point (allowed by Powell's code but not
    % ours), then take a new geometry step, or else take a trust-region step, discarding the previously
    % calculated geometry step D in both cases.
    %
    % 3. If it turns out necessary to call RESCUE again, but no new function value has been evaluated
    % after the last RESCUE, then Powell's code will terminate.
    %
    % Dedicated to the late Professor M. J. D. Powell FRS (1936--2015).
    %
    % Started: February 2022
    %
    % Last Modified: Wed 08 Apr 2026 06:38:40 PM CST
    %--------------------------------------------------------------------------------------------------%

    methods
        function [x, nf, f, fhist, xhist, info] = bobyqb(obj, calfun, iprint, maxfun, npt, eta1, eta2, ftarget, gamma1, gamma2, rhobeg, rhoend, xl, xu, x, fhist, xhist, varargin)
            %--------------------------------------------------------------------------------------------------%
            % This subroutine performs the major calculations of BOBYQA.
            %
            % IPRINT, MAXFUN, MAXHIST, NPT, ETA1, ETA2, FTARGET, GAMMA1, GAMMA2, RHOBEG, RHOEND, XL, XU, X, NF,
            % F, FHIST, XHIST, and INFO are identical to the corresponding arguments in subroutine BOBYQA.
            %
            % XBASE holds a shift of origin that should reduce the contributions from rounding errors to values
            %   of the model and Lagrange functions.
            % SL and SU hold XL - XBASE and XU - XBASE, respectively.
            % XOPT is the displacement from XBASE of the best vector of variables so far (i.e., the one provides
            %   the least calculated F so far). XOPT satisfies SL(I) <= XOPT(I) <= SU(I), with appropriate
            %   equalities when XOPT is on a constraint boundary. FOPT = F(XOPT + XBASE). However, we do not
            %   save XOPT and FOPT explicitly, because XOPT = XPT(:, KOPT) and FOPT = FVAL(KOPT), which is
            %   explained below.
            % [XPT, FVAL, KOPT] describes the interpolation set:
            % XPT contains the interpolation points relative to XBASE, each COLUMN for a point; FVAL holds the
            %   values of F at the interpolation points; KOPT is the index of XOPT in XPT.
            % [GOPT, HQ, PQ] describes the quadratic model: GOPT will hold the gradient of the quadratic model
            %   at XBASE + XOPT; HQ will hold the explicit second order derivatives of the quadratic model; PQ
            %   will contain the parameters of the implicit second order derivatives of the quadratic model.
            % [BMAT, ZMAT] describes the matrix H in the BOBYQA paper (eq. 2.7), which is the inverse of
            %   the coefficient matrix of the KKT system for the least-Frobenius norm interpolation problem:
            % ZMAT will hold a factorization of the leading NPT*NPT submatrix of H, the factorization being
            %   OMEGA = ZMAT*ZMAT^T, which provides both the correct rank and positive semi-definiteness. BMAT
            %   will hold the last N ROWs of H except for the (NPT+1)th column. Note that the (NPT + 1)th row
            %   and column of H are not saved as they are unnecessary for the calculation.
            % D is reserved for trial steps from XOPT. It is chosen by subroutine TRSBOX or GEOSTEP. Usually
            %   XBASE + XOPT + D is the vector of variables for the next call of CALFUN.
            %--------------------------------------------------------------------------------------------------%

            % Common modules
            checkexit_obj = prima_mat.common.checkexit_mod();
            consts_obj = prima_mat.common.consts_mod();
            debug_obj = prima_mat.common.debug_mod(); %, wassert, validate
            evaluate_obj = prima_mat.common.evaluate_mod();
            history_obj = prima_mat.common.history_mod();
            infnan_obj = prima_mat.common.infnan_mod();
            infos_obj = prima_mat.common.infos_mod();
            linalg_obj = prima_mat.common.linalg_mod();
            message_obj = prima_mat.common.message_mod();
            prima_mat.common.pintrf_mod();
            powalg_obj = prima_mat.common.powalg_mod(); %, errquad
            ratio_obj = prima_mat.common.ratio_mod();
            redrho_obj = prima_mat.common.redrho_mod();
            shiftbase_obj = prima_mat.common.shiftbase_mod();
            xinbd_obj = prima_mat.common.xinbd_mod();

            % Solver-specific modules
            geometry_bobyqa_obj = prima_mat.bobyqa.geometry_bobyqa_mod();
            initialize_bobyqa_obj = prima_mat.bobyqa.initialize_bobyqa_mod();
            rescue_obj = prima_mat.bobyqa.rescue_mod();
            trustregion_bobyqa_obj = prima_mat.bobyqa.trustregion_bobyqa_mod();
            update_bobyqa_obj = prima_mat.bobyqa.update_bobyqa_mod();


            % Inputs
            % N.B.: INTENT cannot be specified if a dummy procedure is not a POINTER



            % XL(N)
            % XU(N)

            % In-outputs
            % X(N)

            % Outputs



            % FHIST(MAXFHIST)
            % XHIST(N, MAXXHIST)

            % Local variables
            solver = "BOBYQA";
            srname = "BOBYQB";
            ij = NaN(2, max(0, npt - 2 * numel(x) - 1));


            accurate_mod = false;
            adequate_geo = false;
            bad_trstep = false;
            close_itpset = false;
            improve_geo = false;
            reduce_rho = false;


            small_trrad = false;

            to_rescue = false;

            ximproved = false;
            bmat = NaN(numel(x), npt + numel(x));

            d = NaN(numel(x), 1);
            delbar = NaN;

            den = NaN(npt, 1);
            distsq = NaN(npt, 1);
            dnorm = NaN;
            dnorm_rec = NaN(2, 1); % Powell's implementation: DNORM_REC(3)

            fval = NaN(npt, 1);

            gopt = NaN(numel(x), 1);
            hq = NaN(numel(x));
            moderr = NaN;
            moderr_rec = NaN(numel(dnorm_rec), 1);
            pq = NaN(npt, 1);


            sl = NaN(numel(x), 1);
            su = NaN(numel(x), 1);
            vlag = NaN(npt + numel(x), 1);
            xbase = NaN(numel(x), 1);
            xdrop = NaN(numel(x), 1);
            xosav = NaN(numel(x), 1);
            xpt = NaN(numel(x), npt);
            zmat = NaN(npt, npt - numel(x) + -1);
            trtol = 1.0e-2; % Convergence tolerance of trust-region subproblem solver

            % Sizes.
            n = numel(x);
            maxxhist = size(xhist, 2);
            maxfhist = numel(fhist);
            maxhist = max(maxxhist, maxfhist);

            % Preconditions
            if consts_obj.DEBUGGING
                debug_obj.assert(abs(iprint) <= 3, "IPRINT is 0, 1, -1, 2, -2, 3, or -3", srname);
                debug_obj.assert(n >= 1, "N >= 1", srname);
                debug_obj.assert(npt >= n + 2, "NPT >= N+2", srname);
                debug_obj.assert(maxfun >= npt + 1, "MAXFUN >= NPT+1", srname);
                debug_obj.assert(eta1 >= 0 && eta1 <= eta2 && eta2 < 1, "0 <= ETA1 <= ETA2 < 1", srname);
                debug_obj.assert(gamma1 > 0 && gamma1 < 1 && gamma2 > 1, "0 < GAMMA1 < 1 < GAMMA2", srname);
                debug_obj.assert(rhobeg >= rhoend && rhoend > 0, "RHOBEG >= RHOEND > 0", srname);
                debug_obj.assert(numel(xl) == n && numel(xu) == n, "SIZE(XL) == N == SIZE(XU)", srname);
                debug_obj.assert(all(rhobeg <= (xu - xl) ./ consts_obj.TWO, 'all'), "RHOBEG <= MINVAL(XU-XL)/2", srname);
                debug_obj.assert(all(infnan_obj.is_finite(x), 'all'), "X is finite", srname);
                debug_obj.assert(all(x >= xl & (x <= xl | x - xl >= rhobeg), 'all'), "X == XL or X - XL >= RHOBEG", srname);
                debug_obj.assert(all(x <= xu & (x >= xu | xu - x >= rhobeg), 'all'), "X == XU or XU - X >= RHOBEG", srname);
                debug_obj.assert(maxhist >= 0 && maxhist <= maxfun, "0 <= MAXHIST <= MAXFUN", srname);
                debug_obj.assert(size(xhist, 1) == n && maxxhist * (maxxhist - maxhist) == 0, "SIZE(XHIST, 1) == N, SIZE(XHIST, 2) == 0 or MAXHIST", srname);
                debug_obj.assert(maxfhist * (maxfhist - maxhist) == 0, "SIZE(FHIST) == 0 or MAXHIST", srname);
            end

            %====================%
            % Calculation starts %
            %====================%

            % Initialize XBASE, XPT, SL, SU, FVAL, and KOPT, together with the history, NF, and IJ.
            [x, ij, kopt, nf, fhist, fval, sl, su, xbase, xhist, xpt, subinfo] = initialize_bobyqa_obj.initxf(calfun, iprint, maxfun, ftarget, rhobeg, xl, xu, x, ij, fhist, fval, sl, su, xbase, xhist, xpt);

            % Report the current best value, and check if user asks for early termination.
            terminate = false;
            ipObj = inputParser();
            addParameter(ipObj, 'callback_fcn', struct());
            parse(ipObj, varargin{:});
            callback_fcn = ipObj.Results.callback_fcn;
            if ~ismember('callback_fcn', ipObj.UsingDefaults)
                terminate = callback_fcn(xbase + xpt(:, kopt), fval(kopt), nf, 0);
                if terminate
                    subinfo = infos_obj.CALLBACK_TERMINATE;
                end
            end

            % Initialize X and F according to KOPT.
            x(:) = xinbd_obj.xinbd(xbase, xpt(:, kopt), xl, xu, sl, su); % In precise arithmetic, X = XBASE + XOPT.
            f = fval(kopt);

            % Finish the initialization if INITXF completed normally and CALLBACK did not request termination;
            % otherwise, do not proceed, as XPT etc may be uninitialized, leading to errors or exceptions.
            if subinfo == infos_obj.INFO_DFT
                % Initialize [BMAT, ZMAT], representing inverse of KKT matrix of the interpolation system.
                [bmat, zmat] = initialize_bobyqa_obj.inith(ij, xpt, bmat, zmat);

                % Initialize the quadratic represented by [GOPT, HQ, PQ], so that its gradient at XBASE+XOPT is
                % GOPT; its Hessian is HQ + sum_{K=1}^NPT PQ(K)*XPT(:, K)*XPT(:, K)'.
                [gopt, hq, pq] = initialize_bobyqa_obj.initq(ij, fval, xpt, gopt, hq, pq);
                if ~(all(infnan_obj.is_finite(gopt), 'all') && all(infnan_obj.is_finite(hq), 'all') && all(infnan_obj.is_finite(pq), 'all'))
                    subinfo = infos_obj.NAN_INF_MODEL;
                end
            end

            % Check whether to return due to abnormal cases that may occur during the initialization.
            if subinfo ~= infos_obj.INFO_DFT
                info = subinfo;
                % Arrange FHIST and XHIST so that they are in the chronological order.
                [xhist, fhist] = history_obj.rangehist(nf, xhist, fhist);
                % Print a return message according to IPRINT.
                message_obj.retmsg(solver, info, iprint, nf, f, x);
                % Postconditions
                if consts_obj.DEBUGGING
                    debug_obj.assert(nf <= maxfun, "NF <= MAXFUN", srname);
                    debug_obj.assert(numel(x) == n && ~any(infnan_obj.is_nan_sp(x), 'all'), "SIZE(X) == N, X does not contain NaN", srname);
                    debug_obj.assert(all(x >= xl, 'all') && all(x <= xu, 'all'), "XL <= X <= XU", srname);
                    debug_obj.assert(~(infnan_obj.is_nan_sp(f) || infnan_obj.is_posinf(f)), "F is not NaN/+Inf", srname);
                    debug_obj.assert(size(xhist, 1) == n && size(xhist, 2) == maxxhist, "SIZE(XHIST) == [N, MAXXHIST]", srname);
                    debug_obj.assert(~any(infnan_obj.is_nan_sp(xhist(:, 1:min(nf, maxxhist))), 'all'), "XHIST does not contain NaN", srname);
                    % The last calculated X can be Inf (finite + finite can be Inf numerically).
                    for k = 1:min(nf, maxxhist)
                        debug_obj.assert(all(xhist(:, k) >= xl, 'all') && all(xhist(:, k) <= xu, 'all'), "XL <= XHIST <= XU", srname);
                    end
                    debug_obj.assert(numel(fhist) == maxfhist, "SIZE(FHIST) == MAXFHIST", srname);
                    debug_obj.assert(~any(infnan_obj.is_nan_sp(fhist(1:min(nf, maxfhist))) | infnan_obj.is_posinf(fhist(1:min(nf, maxfhist))), 'all'), "FHIST does not contain NaN/+Inf", srname);
                    debug_obj.assert(~any(fhist(1:min(nf, maxfhist)) < f, 'all'), "F is the smallest in FHIST", srname);
                end
                return
            end

            % Set some more initial values.
            % We must initialize RATIO. Otherwise, when SHORTD = TRUE, compilers may raise a run-time error that
            % RATIO is undefined. But its value will not be used: when SHORTD = FALSE, its value will be
            % overwritten; when SHORTD = TRUE, its value is used only in BAD_TRSTEP, which is TRUE regardless of
            % RATIO. Similar for KNEW_TR.
            % No need to initialize SHORTD unless MAXTR < 1, but some compilers may complain if we do not do it.
            rho = rhobeg;
            delta = rho;
            ebound = consts_obj.ZERO;
            rescued = false;
            shortd = false;

            ratio = -consts_obj.ONE;
            dnorm_rec(:) = consts_obj.REALMAX;
            moderr_rec(:) = consts_obj.REALMAX;
            knew_tr = 0;
            knew_geo = 0;
            itest = 0;

            % If DELTA <= GAMMA3*RHO after an update, we set DELTA to RHO. GAMMA3 must be less than GAMMA2. The
            % reason is as follows. Imagine a very successful step with DENORM = the un-updated DELTA = RHO.
            % Then TRRAD will update DELTA to GAMMA2*RHO. If GAMMA3 >= GAMMA2, then DELTA will be reset to RHO,
            % which is not reasonable as D is very successful. See paragraph two of Sec. 5.2.5 in
            % T. M. Ragonneau's thesis: "Model-Based Derivative-Free Optimization Methods and Software".
            % According to test on 20230613, for BOBYQA, this Powellful updating scheme of DELTA works better
            % than setting directly DELTA = MAX(NEW_DELTA, RHO).
            gamma3 = max(consts_obj.ONE, min(0.75 * gamma2, 1.5));

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
            info = infos_obj.MAXTR_REACHED;

            % Begin the iterative procedure.
            % After solving a trust-region subproblem, we use three boolean variables to control the workflow.
            % SHORTD: Is the trust-region trial step too short to invoke a function evaluation?
            % IMPROVE_GEO: Should we improve the geometry?
            % REDUCE_RHO: Should we reduce rho?
            % BOBYQA never sets IMPROVE_GEO and REDUCE_RHO to TRUE simultaneously.
            for tr = 1:maxtr
                % Generate the next trust region step D.
                [crvmin, d] = trustregion_bobyqa_obj.trsbox(delta, gopt, hq, pq, sl, su, trtol, xpt(:, kopt), xpt, d);
                dnorm = min(delta, linalg_obj.p_norm(d));
                shortd = (dnorm <= consts_obj.HALF * rho); % `<=` works better than `<` in case of underflow.

                % Set QRED to the reduction of the quadratic model when the move D is made from XOPT. QRED
                % should be positive. If it is nonpositive due to rounding errors, we will not take this step.
                qred = -powalg_obj.quadinc_d0(d, xpt, gopt, pq, 'hq', hq); % QRED = Q(XOPT) - Q(XOPT + D)
                trfail = (~(qred > 1.0e-6 * rho ^ 2)); % QRED is tiny/negative or NaN.

                % When D is short, make a choice between reducing RHO and improving the geometry depending
                % on whether or not our work with the current RHO seems complete. RHO is reduced if the
                % errors in the quadratic model at the recent interpolation points compare favourably
                % with predictions of likely improvements to the model within distance HALF*RHO of XOPT.
                % Why do we reduce RHO when SHORTD is true and the entries of MODERR_REC and DNORM_REC are all
                % small? The reason is well explained by the BOBYQA paper in the paragraphs surrounding
                % (6.8)--(6.11). Roughly speaking, in this case, a trust-region step is unlikely to decrease the
                % objective function according to some estimations. This suggests that the current trust-region
                % center may be an approximate local minimizer up to the current "resolution" of the algorithm.
                % When this occurs, the algorithm takes the view that the work for the current RHO is complete,
                % and hence it will reduce RHO, which will enhance the resolution of the algorithm in general.
                if shortd || trfail
                    delta = consts_obj.TENTH * delta;
                    if delta <= gamma3 * rho
                        delta = rho; % Set DELTA to RHO when it is close to or below.

                    end
                    % Evaluate EBOUND. It will be used as a bound to test if the entries of MODERR_REC are small.
                    ebound = obj.errbd(crvmin, d, gopt, hq, moderr_rec, pq, rho, sl, su, xpt(:, kopt), xpt);
                else
                    % Calculate the next value of the objective function.
                    x(:) = xinbd_obj.xinbd(xbase, xpt(:, kopt) + d, xl, xu, sl, su); % X = XBASE + XOPT + D without rounding.
                    f = evaluate_obj.evaluatef(calfun, x);
                    nf = nf + 1;
                    rescued = false; % Set RESCUED to FALSE after evaluating F at a new point.

                    % Print a message about the function evaluation according to IPRINT.
                    message_obj.fmsg(solver, "Trust region", iprint, nf, delta, f, x);
                    % Save X, F into the history.
                    [xhist, fhist] = history_obj.savehist(nf, x, xhist, f, fhist);

                    % Check whether to exit
                    subinfo = checkexit_obj.checkexit_unc(maxfun, nf, f, ftarget, x);
                    if subinfo ~= infos_obj.INFO_DFT
                        info = subinfo;
                        break
                    end

                    % Update DNORM_REC and MODERR_REC.
                    % DNORM_REC records the DNORM of the recent function evaluations with the current RHO.
                    dnorm_rec(:) = [reshape(dnorm_rec(2:numel(dnorm_rec)), [], 1); dnorm];
                    % MODERR is the error of the current model in predicting the change in F due to D.
                    % MODERR_REC records the prediction errors of the recent models with the current RHO.
                    moderr = f - fval(kopt) + qred;
                    moderr_rec(:) = [reshape(moderr_rec(2:numel(moderr_rec)), [], 1); moderr];

                    % Calculate the reduction ratio by REDRAT, which handles Inf/NaN carefully.
                    ratio = ratio_obj.redrat(fval(kopt) - f, qred, eta1);

                    % Update DELTA. After this, DELTA < DNORM may hold.
                    delta = trustregion_bobyqa_obj.trrad(delta, dnorm, eta1, eta2, gamma1, gamma2, ratio);
                    if delta <= gamma3 * rho
                        delta = rho; % Set DELTA to RHO when it is close to or below.

                    end

                    % Is the newly generated X better than current best point?
                    ximproved = (f < fval(kopt));

                    % Call RESCUE if rounding errors have damaged the denominator corresponding to D.
                    % RESCUE is invoked sometimes though not often after a trust-region step, and it does
                    % improve the performance, especially when pursing high-precision solutions.
                    vlag(:) = powalg_obj.calvlag_lfqint(kopt, bmat, d, xpt, zmat);
                    den(:) = powalg_obj.calden(kopt, bmat, d, xpt, zmat);
                    to_rescue = (ximproved && ~(infnan_obj.is_finite(sum(abs(vlag), 'all')) && any(den > max(vlag(1:npt) .^ 2, [], 'all'), 'all')));
                    % Below are some alternatives conditions for calling RESCUE. They perform fairly well.
                    % %to_rescue = .false.  ! Do not call RESCUE at all.
                    % %to_rescue = (ximproved .and. .not. any(den > 0.25_RP * maxval(vlag(1:npt)**2)))
                    % %to_rescue = (ximproved .and. .not. any(den > HALF * maxval(vlag(1:npt)**2)))
                    % %to_rescue = (.not. any(den > HALF * maxval(vlag(1:npt)**2)))  ! Powell's code.
                    % %to_rescue = (.not. any(den > maxval(vlag(1:npt)**2)))
                    if to_rescue
                        if rescued
                            info = infos_obj.DAMAGING_ROUNDING; % The last RESCUE did not improve the situation.
                            break
                        end
                        [kopt, nf, fhist, fval, gopt, hq, pq, sl, su, xbase, xhist, xpt, bmat, zmat, subinfo] = rescue_obj.rescue(calfun, solver, iprint, maxfun, delta, ftarget, xl, xu, kopt, nf, fhist, fval, gopt, hq, pq, sl, su, xbase, xhist, xpt, bmat, zmat);
                        if subinfo ~= infos_obj.INFO_DFT
                            info = subinfo;
                            break
                        end
                        rescued = true;
                        dnorm_rec(:) = consts_obj.REALMAX;
                        moderr_rec(:) = consts_obj.REALMAX;

                        % RESCUE shifts XBASE to the best point before RESCUE. Update D, MODERR, and XIMPROVED.
                        % Do NOT calculate QRED according to this D, as it is not really a trust region step.
                        % Note that QRED will be used afterward for defining IMPROVE_GEO and REDUCE_RHO.
                        d(:) = max(sl, min(su, d)) - xpt(:, kopt);
                        moderr = f - fval(kopt) - powalg_obj.quadinc_d0(d, xpt, gopt, pq, 'hq', hq);
                        ximproved = (f < fval(kopt));
                    end

                    % Set KNEW_TR to the index of the interpolation point to be replaced with XOPT + D.
                    % KNEW_TR will ensure that the geometry of XPT is "good enough" after the replacement.
                    knew_tr = geometry_bobyqa_obj.setdrop_tr(kopt, ximproved, bmat, d, delta, rho, xpt, zmat);

                    % Update [BMAT, ZMAT] (representing H in the BOBYQA paper), [GQ, HQ, PQ] (the quadratic
                    % model), and [FVAL, XPT, KOPT, FOPT, XOPT] so that XPT(:, KNEW_TR) becomes XOPT + D. If
                    % KNEW_TR = 0, the updating subroutines will do essentially nothing, as the algorithm
                    % decides not to include XOPT + D into XPT.
                    if knew_tr > 0
                        xdrop(:) = xpt(:, knew_tr);
                        xosav(:) = xpt(:, kopt);
                        [bmat, zmat] = update_bobyqa_obj.updateh(knew_tr, kopt, d, xpt, bmat, zmat);
                        [kopt, fval, xpt] = update_bobyqa_obj.updatexf(knew_tr, ximproved, f, max(sl, min(su, xosav + d)), kopt, fval, xpt);
                        [gopt, hq, pq] = update_bobyqa_obj.updateq(knew_tr, ximproved, bmat, d, moderr, xdrop, xosav, xpt, zmat, gopt, hq, pq);
                        % Try whether to replace the new quadratic model with the alternative model, namely the
                        % least Frobenius norm interpolant.
                        [itest, gopt, hq, pq] = update_bobyqa_obj.tryqalt(bmat, fval - fval(kopt), ratio, sl, su, xpt(:, kopt), xpt, zmat, itest, gopt, hq, pq);
                        if ~(all(infnan_obj.is_finite(gopt), 'all') && all(infnan_obj.is_finite(hq), 'all') && all(infnan_obj.is_finite(pq), 'all'))
                            info = infos_obj.NAN_INF_MODEL;
                            break
                        end
                    end
                end % End of IF (SHORTD .OR. TRFAIL). The normal trust-region calculation ends.


                %----------------------------------------------------------------------------------------------%
                % Before the next trust-region iteration, we may improve the geometry of XPT or reduce RHO
                % according to IMPROVE_GEO and REDUCE_RHO, which in turn depend on the following indicators.
                % N.B.: We must ensure that the algorithm does not set IMPROVE_GEO = TRUE at infinitely many
                % consecutive iterations without moving XOPT or reducing RHO. Otherwise, the algorithm will get
                % stuck in repetitive invocations of GEOSTEP. To this end, make sure the following.
                % 1. The threshold for CLOSE_ITPSET is at least DELBAR, the trust region radius for GEOSTEP.
                % Normally, DELBAR <= DELTA <= the threshold (In Powell's UOBYQA, DELBAR = RHO < the threshold).
                % 2. If an iteration sets IMPROVE_GEO = TRUE, it must also reduce DELTA or set DELTA to RHO.

                % ACCURATE_MOD: Are the recent models sufficiently accurate? Used only if SHORTD is TRUE.
                accurate_mod = all(abs(moderr_rec) <= ebound, 'all') && all(dnorm_rec <= rho, 'all');
                % CLOSE_ITPSET: Are the interpolation points close to XOPT?
                distsq(:) = sum((xpt - xpt(:, kopt)) .^ 2, 1);
                %%MATLAB: distsq = sum((xpt - xpt(:, kopt)).^2)  % Implicit expansion
                close_itpset = all(distsq <= max(delta ^ 2, (consts_obj.TEN * rho) ^ 2), 'all');
                % Below are some alternative definitions of CLOSE_ITPSET.
                % N.B.: The threshold for CLOSE_ITPSET is at least DELBAR, the trust region radius for GEOSTEP.
                % %close_itpset = all(distsq <= max((TWO * delta)**2, (TEN * rho)**2))  ! Powell's code.
                % %close_itpset = all(distsq <= 4.0_RP * delta**2)  ! Powell's NEWUOA code.
                % %close_itpset = all(distsq <= max(delta**2, 4.0_RP * rho**2))  ! Powell's LINCOA code.
                % ADEQUATE_GEO: Is the geometry of the interpolation set "adequate"?
                % N.B. (Zaikun 20240314): Even if RESCUE has just been called (RESCUED = TRUE), the geometry may
                % still be inadequate/improvable if XPT contains points far away from XOPT.
                adequate_geo = (shortd && accurate_mod) || close_itpset;
                % SMALL_TRRAD: Is the trust-region radius small? This indicator seems not impactive in practice.
                small_trrad = (max(delta, dnorm) <= rho); % Powell's code. See also (6.7) of the BOBYQA paper.
                %small_trrad = (delsav <= rho)  ! Behaves the same as Powell's version. DELSAV = unupdated DELTA.

                % IMPROVE_GEO and REDUCE_RHO are defined as follows.
                % N.B.: If SHORTD is TRUE at the very first iteration, then REDUCE_RHO will be set to TRUE.
                % Powell's code does not have TRFAIL in BAD_TRSTEP; it terminates if TRFAIL is TRUE.

                % BAD_TRSTEP (for IMPROVE_GEO): Is the last trust-region step bad?
                bad_trstep = (shortd || trfail || ratio <= eta1 || knew_tr == 0);
                improve_geo = bad_trstep && ~adequate_geo; % See the text above (6.7) of the BOBYQA paper.
                % BAD_TRSTEP (for REDUCE_RHO): Is the last trust-region step bad?
                bad_trstep = (shortd || trfail || ratio <= 0 || knew_tr == 0);
                reduce_rho = bad_trstep && adequate_geo && small_trrad; % See (6.7) of the BOBYQA paper.
                % Zaikun 20221111: What if RESCUE has been called? Is it still reasonable to use RATIO?
                % Zaikun 20221127: If RESCUE has been called, then KNEW_TR may be 0 even if RATIO > 0.

                % Equivalently, REDUCE_RHO can be set as follows. It shows that REDUCE_RHO is TRUE in two cases.
                % %bad_trstep = (shortd .or. trfail .or. ratio <= 0 .or. knew_tr == 0)
                % %reduce_rho = (shortd .and. accurate_mod) .or. (bad_trstep .and. close_itpset .and. small_trrad)

                % With REDUCE_RHO properly defined, we can also set IMPROVE_GEO as follows.
                % %bad_trstep = (shortd .or. trfail .or. ratio <= eta1 .or. knew_tr == 0)
                % %improve_geo = bad_trstep .and. (.not. reduce_rho) .and. (.not. close_itpset)

                % With IMPROVE_GEO properly defined, we can also set REDUCE_RHO as follows.
                % %bad_trstep = (shortd .or. trfail .or. ratio <= 0 .or. knew_tr == 0)
                % %reduce_rho = bad_trstep .and. (.not. improve_geo) .and. small_trrad

                % BOBYQA never sets IMPROVE_GEO and REDUCE_RHO to TRUE simultaneously.
                %call assert(.not. (improve_geo .and. reduce_rho), 'IMPROVE_GEO and REDUCE_RHO are not both TRUE', srname)
                %
                % If SHORTD or TRFAIL is TRUE, then either IMPROVE_GEO or REDUCE_RHO is TRUE unless CLOSE_ITPSET
                % is TRUE but SMALL_TRRAD is FALSE.
                %call assert((.not. (shortd .or. trfail)) .or. (improve_geo .or. reduce_rho .or. &
                %    & (close_itpset .and. .not. small_trrad)), 'If SHORTD or TRFAIL is TRUE, then either &
                %    & IMPROVE_GEO or REDUCE_RHO is TRUE unless CLOSE_ITPSET is TRUE but SMALL_TRRAD is FALSE', srname)
                %----------------------------------------------------------------------------------------------%


                % Since IMPROVE_GEO and REDUCE_RHO are never TRUE simultaneously, the following two blocks are
                % exchangeable: IF (IMPROVE_GEO) ... END IF and IF (REDUCE_RHO) ... END IF.

                % Improve the geometry of the interpolation set by removing a point and adding a new one.
                if improve_geo
                    % XPT(:, KNEW_GEO) will become XOPT + D below. KNEW_GEO /= KOPT unless there is a bug.
                    knew_geo = fortran.maxloc(distsq, 'dim', 1);

                    % Set DELBAR, which will be used as the trust-region radius for the geometry-improving
                    % scheme GEOSTEP. Note that DELTA has been updated before arriving here.
                    delbar = max(min(consts_obj.TENTH * sqrt(max(distsq, [], 'all')), delta), rho); % Powell's code
                    %delbar = rho  ! Powell's UOBYQA code
                    %delbar = max(min(TENTH * sqrt(maxval(distsq)), HALF * delta), rho)  ! Powell's NEWUOA code
                    %delbar = max(TENTH * delta, rho)  ! Powell's LINCOA code

                    % Find D so that the geometry of XPT will be improved when XPT(:, KNEW_GEO) becomes XOPT + D.
                    d(:) = geometry_bobyqa_obj.geostep(knew_geo, kopt, bmat, delbar, sl, su, xpt, zmat);

                    % Call RESCUE if rounding errors have damaged the denominator corresponding to D.
                    % 1. This does make a difference, yet RESCUE seems not invoked often after a geometry step.
                    % 2. In Powell's implementation, it may happen that RESCUE only recalculates [BMAT, ZMAT]
                    % without introducing any new point into XPT. In that case, GEOSTEP will have to be called
                    % after RESCUE, without which the code may encounter an infinite cycling. We have modified
                    % RESCUE so that it introduces at least one new point into XPT and there is no need to call
                    % GEOSTEP afterward. This improves the performance a bit and simplifies the flow of the code.
                    % 3. It is tempting to incorporate XOPT+D into the interpolation even if RESCUE is called.
                    % However, this cannot be done without recalculating KNEW_GEO, as XPT has been changed by
                    % RESCUE, so that it is invalid to replace XPT(:, KNEW_GEO) with XOPT+D anymore. With a new
                    % KNEW_GEO, the step D will become improper as it was chosen according to the old KNEW_GEO.
                    vlag(:) = powalg_obj.calvlag_lfqint(kopt, bmat, d, xpt, zmat);
                    den(:) = powalg_obj.calden(kopt, bmat, d, xpt, zmat);
                    to_rescue = (~(infnan_obj.is_finite(sum(abs(vlag), 'all')) && den(knew_geo) > consts_obj.HALF * vlag(knew_geo) ^ 2));
                    if to_rescue
                        if rescued
                            info = infos_obj.DAMAGING_ROUNDING; % The last RESCUE did not improve the situation.
                            break
                        end
                        [kopt, nf, fhist, fval, gopt, hq, pq, sl, su, xbase, xhist, xpt, bmat, zmat, subinfo] = rescue_obj.rescue(calfun, solver, iprint, maxfun, delta, ftarget, xl, xu, kopt, nf, fhist, fval, gopt, hq, pq, sl, su, xbase, xhist, xpt, bmat, zmat);
                        if subinfo ~= infos_obj.INFO_DFT
                            info = subinfo;
                            break
                        end
                        rescued = true;
                        dnorm_rec(:) = consts_obj.REALMAX;
                        moderr_rec(:) = consts_obj.REALMAX;
                    else
                        % Calculate the next value of the objective function.
                        x(:) = xinbd_obj.xinbd(xbase, xpt(:, kopt) + d, xl, xu, sl, su); % X = XBASE + XOPT + D without rounding.
                        f = evaluate_obj.evaluatef(calfun, x);
                        nf = nf + 1;
                        rescued = false; % Set RESCUED to FALSE after evaluating F at a new point.

                        % Print a message about the function evaluation according to IPRINT.
                        message_obj.fmsg(solver, "Geometry", iprint, nf, delbar, f, x);
                        % Save X, F into the history.
                        [xhist, fhist] = history_obj.savehist(nf, x, xhist, f, fhist);

                        % Check whether to exit
                        subinfo = checkexit_obj.checkexit_unc(maxfun, nf, f, ftarget, x);
                        if subinfo ~= infos_obj.INFO_DFT
                            info = subinfo;
                            break
                        end

                        % Update DNORM_REC and MODERR_REC.
                        % DNORM_REC records the DNORM of the recent function evaluations with the current RHO.
                        % Powell's code does not update DNORM. Therefore, DNORM is the length of the last
                        % trust-region trial step, inconsistent with MODERR_REC. The same problem exists in NEWUOA.
                        dnorm = min(delbar, linalg_obj.p_norm(d));
                        dnorm_rec(:) = [reshape(dnorm_rec(2:numel(dnorm_rec)), [], 1); dnorm];
                        % MODERR is the error of the current model in predicting the change in F due to D.
                        % MODERR_REC records the prediction errors of the recent models with the current RHO.
                        moderr = f - fval(kopt) - powalg_obj.quadinc_d0(d, xpt, gopt, pq, 'hq', hq); % QRED = Q(XOPT) - Q(XOPT + D)
                        moderr_rec(:) = [reshape(moderr_rec(2:numel(moderr_rec)), [], 1); moderr];

                        % Is the newly generated X better than current best point?
                        ximproved = (f < fval(kopt));

                        % Update [BMAT, ZMAT] (represents H in the BOBYQA paper), [FVAL, XPT, KOPT, FOPT, XOPT],
                        % and [GQ, HQ, PQ] (the quadratic model), so that XPT(:, KNEW_GEO) becomes XOPT + D.
                        xdrop(:) = xpt(:, knew_geo);
                        xosav(:) = xpt(:, kopt);
                        [bmat, zmat] = update_bobyqa_obj.updateh(knew_geo, kopt, d, xpt, bmat, zmat);
                        [kopt, fval, xpt] = update_bobyqa_obj.updatexf(knew_geo, ximproved, f, max(sl, min(su, xosav + d)), kopt, fval, xpt);
                        [gopt, hq, pq] = update_bobyqa_obj.updateq(knew_geo, ximproved, bmat, d, moderr, xdrop, xosav, xpt, zmat, gopt, hq, pq);
                        if ~(all(infnan_obj.is_finite(gopt), 'all') && all(infnan_obj.is_finite(hq), 'all') && all(infnan_obj.is_finite(pq), 'all'))
                            info = infos_obj.NAN_INF_MODEL;
                            break
                        end
                    end
                end % End of IF (IMPROVE_GEO). The procedure of improving geometry ends.

                % The calculations with the current RHO are complete. Enhance the resolution of the algorithm
                % by reducing RHO; update DELTA at the same time.
                if reduce_rho
                    if rho <= rhoend
                        info = infos_obj.SMALL_TR_RADIUS;
                        break
                    end
                    delta = max(consts_obj.HALF * rho, redrho_obj.redrho(rho, rhoend));
                    rho = redrho_obj.redrho(rho, rhoend);
                    % Print a message about the reduction of RHO according to IPRINT.
                    message_obj.rhomsg(solver, iprint, nf, delta, fval(kopt), rho, xbase + xpt(:, kopt));
                    % DNORM_REC and MODERR_REC are corresponding to the recent function evaluations with
                    % the current RHO. Update them after reducing RHO.
                    dnorm_rec(:) = consts_obj.REALMAX;
                    moderr_rec(:) = consts_obj.REALMAX;
                end % End of IF (REDUCE_RHO). The procedure of reducing RHO ends.

                % Shift XBASE if XOPT may be too far from XBASE.
                % Powell's original criteria for shifting XBASE is as follows.
                % 1. After a trust region step that is not short, shift XBASE if SUM(XOPT**2) >= 1.0E3*DNORM**2.
                % In this case, it seems quite important for the performance to recalculate QRED.
                % 2. Before a geometry step, shift XBASE if SUM(XOPT**2) >= 1.0E3*DELBAR**2.
                if sum(xpt(:, kopt) .^ 2, 'all') >= 1000.0 * delta ^ 2
                    % Other possible criteria: SUM(XOPT**2) >= 1.0E4*DELTA**2, SUM(XOPT**2) >= 1.0E4*RHO**2.
                    sl(:) = min(sl - xpt(:, kopt), consts_obj.ZERO);
                    su(:) = max(su - xpt(:, kopt), consts_obj.ZERO);
                    [xbase, xpt, bmat, hq] = shiftbase_obj.shiftbase_lfqint(kopt, xbase, xpt, zmat, bmat, pq, hq);
                    xbase(:) = max(xl, min(xu, xbase));
                end

                % Report the current best value, and check if user asks for early termination.
                if ~ismember('callback_fcn', ipObj.UsingDefaults)
                    terminate = callback_fcn(xbase + xpt(:, kopt), fval(kopt), nf, tr);
                    if terminate
                        info = infos_obj.CALLBACK_TERMINATE;
                        break
                    end
                end

            end % End of DO TR = 1, MAXTR. The iterative procedure ends.

            % Return from the calculation, after trying the Newton-Raphson step if it has not been tried yet.
            if info == infos_obj.SMALL_TR_RADIUS && shortd && dnorm > consts_obj.TENTH * rhoend && nf < maxfun
                x(:) = xinbd_obj.xinbd(xbase, xpt(:, kopt) + d, xl, xu, sl, su); % In precise arithmetic, X = XBASE + XOPT + D.
                f = evaluate_obj.evaluatef(calfun, x);
                nf = nf + 1;
                % Print a message about the function evaluation according to IPRINT.
                % Zaikun 20230512: DELTA has been updated. RHO is only indicative here. TO BE IMPROVED.
                message_obj.fmsg(solver, "Trust region", iprint, nf, rho, f, x);
                % Save X, F into the history.
                [xhist, fhist] = history_obj.savehist(nf, x, xhist, f, fhist);
            end

            % Choose the [X, F] to return: either the current [X, F] or [XBASE + XOPT, FOPT].
            if fval(kopt) < f || infnan_obj.is_nan_sp(f)
                x(:) = xinbd_obj.xinbd(xbase, xpt(:, kopt), xl, xu, sl, su); % In precise arithmetic, X = XBASE + XOPT.
                f = fval(kopt);
            end

            % Arrange FHIST and XHIST so that they are in the chronological order.
            [xhist, fhist] = history_obj.rangehist(nf, xhist, fhist);

            % Print a return message according to IPRINT.
            message_obj.retmsg(solver, info, iprint, nf, f, x);

            %====================%
            %  Calculation ends  %
            %====================%

            % Postconditions
            if consts_obj.DEBUGGING
                debug_obj.assert(nf <= maxfun, "NF <= MAXFUN", srname);
                debug_obj.assert(numel(x) == n && ~any(infnan_obj.is_nan_sp(x), 'all'), "SIZE(X) == N, X does not contain NaN", srname);
                debug_obj.assert(all(x >= xl, 'all') && all(x <= xu, 'all'), "XL <= X <= XU", srname);
                debug_obj.assert(~(infnan_obj.is_nan_sp(f) || infnan_obj.is_posinf(f)), "F is not NaN/+Inf", srname);
                debug_obj.assert(size(xhist, 1) == n && size(xhist, 2) == maxxhist, "SIZE(XHIST) == [N, MAXXHIST]", srname);
                debug_obj.assert(~any(infnan_obj.is_nan_sp(xhist(:, 1:min(nf, maxxhist))), 'all'), "XHIST does not contain NaN", srname);
                % The last calculated X can be Inf (finite + finite can be Inf numerically).
                for k = 1:min(nf, maxxhist)
                    debug_obj.assert(all(xhist(:, k) >= xl, 'all') && all(xhist(:, k) <= xu, 'all'), "XL <= XHIST <= XU", srname);
                end
                debug_obj.assert(numel(fhist) == maxfhist, "SIZE(FHIST) == MAXFHIST", srname);
                debug_obj.assert(~any(infnan_obj.is_nan_sp(fhist(1:min(nf, maxfhist))) | infnan_obj.is_posinf(fhist(1:min(nf, maxfhist))), 'all'), "FHIST does not contain NaN/+Inf", srname);
                debug_obj.assert(~any(fhist(1:min(nf, maxfhist)) < f, 'all'), "F is the smallest in FHIST", srname);
            end

        end
        function ebound = errbd(~, crvmin, d, gopt, hq, moderr_rec, pq, rho, sl, su, xopt, xpt)
            %--------------------------------------------------------------------------------------------------%
            % This function defines EBOUND, which will be used as a bound to test whether the errors in recent
            % models are sufficiently small. See the elaboration on pages 30--31 of the BOBYQA paper, in the
            % paragraphs surrounding (6.8)--(6.11).
            %--------------------------------------------------------------------------------------------------%

            % Common modules
            consts_obj = prima_mat.common.consts_mod();
            debug_obj = prima_mat.common.debug_mod();
            infnan_obj = prima_mat.common.infnan_mod();
            linalg_obj = prima_mat.common.linalg_mod();
            powalg_obj = prima_mat.common.powalg_mod();


            % Inputs



            % Outputs
            ebound = NaN;

            % Local variables
            srname = "ERRBD";


            bfirst = NaN(numel(d), 1);
            bsecond = NaN(numel(d), 1);
            gnew = NaN(numel(d), 1);
            xnew = NaN(numel(d), 1);

            % Sizes
            n = size(xpt, 1);
            npt = size(xpt, 2);

            % Preconditions
            if consts_obj.DEBUGGING
                debug_obj.assert(n >= 1 && npt >= n + 2, "N >= 1, NPT >= N + 2", srname);
                debug_obj.assert(crvmin >= 0, "CRVMIN >= 0", srname);
                debug_obj.assert(numel(d) == n && all(infnan_obj.is_finite(d), 'all'), "SIZE(D) == N, D is finite", srname);
                debug_obj.assert(numel(gopt) == n, "SIZE(GOPT) == N", srname);
                debug_obj.assert(size(hq, 1) == n && linalg_obj.issymmetric(hq), "HQ is n-by-n and symmetric", srname);
                debug_obj.assert(numel(pq) == npt, "SIZE(PQ) == NPT", srname);
                debug_obj.assert(rho > 0, "RHO > 0", srname);
                debug_obj.assert(numel(sl) == n && numel(su) == n, "SIZE(SL) == N == SIZE(SU)", srname);
                debug_obj.assert(numel(xopt) == n && all(infnan_obj.is_finite(xopt), 'all'), "SIZE(XOPT) == N, XOPT is finite", srname);
                debug_obj.assert(all(infnan_obj.is_finite(xpt), 'all'), "XPT is finite", srname);
                debug_obj.assert(all(xopt >= sl & xopt <= su, 'all'), "SL <= XOPT <= SU", srname);
                debug_obj.assert(all(xpt >= sl & xpt <= su, 'all'), "SL <= XPT <= SU", srname);
            end

            %====================%
            % Calculation starts %
            %====================%

            xnew(:) = xopt + d;
            gnew(:) = gopt + powalg_obj.hess_mul(d, xpt, pq, 'hq', hq);
            bfirst(:) = max(abs(moderr_rec), [], 'all');
            bfirst(linalg_obj.trueloc(xnew <= sl)) = gnew(linalg_obj.trueloc(xnew <= sl)) * rho;
            bfirst(linalg_obj.trueloc(xnew >= su)) = -gnew(linalg_obj.trueloc(xnew >= su)) * rho;
            bsecond(:) = consts_obj.HALF * (linalg_obj.diag(hq) + linalg_obj.matprod21(xpt .^ 2, pq)) * rho ^ 2;
            ebound = min(max(bfirst, bfirst + bsecond), [], 'all');
            if crvmin > 0
                ebound = min(ebound, 0.125 * crvmin * rho ^ 2);
            end

            %====================%
            %  Calculation ends  %
            %====================%

        end

    end
end