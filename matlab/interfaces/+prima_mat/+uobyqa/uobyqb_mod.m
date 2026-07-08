classdef uobyqb_mod
    %--------------------------------------------------------------------------------------------------%
    % This module performs the major calculations of UOBYQA.
    %
    % Coded by Zaikun ZHANG (www.zhangzk.net) based on Powell's code and the UOBYQA paper.
    %
    % Dedicated to the late Professor M. J. D. Powell FRS (1936--2015).
    %
    % Started: February 2022
    %
    % Last Modified: Wed 08 Apr 2026 06:39:00 PM CST
    %--------------------------------------------------------------------------------------------------%

    methods
        function [x, nf, f, fhist, xhist, info] = uobyqb(~, calfun, iprint, maxfun, eta1, eta2, ftarget, gamma1, gamma2, rhobeg, rhoend, x, fhist, xhist, varargin)
            %--------------------------------------------------------------------------------------------------%
            % This subroutine performs the major calculations of UOBYQA.
            %
            % The arguments N, X, RHOBEG, RHOEND, IPRINT and MAXFUN are identical to the corresponding arguments
            % in subroutine UOBYQA.
            %
            % XBASE will contain a shift of origin that reduces the contributions from rounding errors to values
            %   of the model and Lagrange functions.
            % XBASE holds a shift of origin that should reduce the contributions from rounding errors to values
            %   of the model and Lagrange functions.
            % XOPT is the displacement from XBASE of the best vector of variables so far (i.e., the one provides
            %   the least calculated F so far). FOPT = F(XOPT + XBASE). However, we do not save XOPT and FOPT
            %   explicitly, because XOPT = XPT(:, KOPT) and FOPT = FVAL(KOPT), which is explained below.
            % [XPT, FVAL, KOPT] describes the interpolation set:
            % XPT contains the interpolation points relative to XBASE, each COLUMN for a point; FVAL holds the
            %   values of F at the interpolation points; KOPT is the index of XOPT in XPT.
            % PQ will contain the parameters of the quadratic model.
            % PL will contain the parameters of the Lagrange functions.
            % D is reserved for trial steps from XOPT. It is chosen by subroutine TRSTEP or GEOSTEP. Usually
            %   XBASE + XOPT + D is the vector of variables for the next call of CALFUN.
            %--------------------------------------------------------------------------------------------------%

            % Common modules
            checkexit_obj = prima_mat.common.checkexit_mod();
            consts_obj = prima_mat.common.consts_mod();
            debug_obj = prima_mat.common.debug_mod();
            evaluate_obj = prima_mat.common.evaluate_mod();
            history_obj = prima_mat.common.history_mod();
            infnan_obj = prima_mat.common.infnan_mod();
            infos_obj = prima_mat.common.infos_mod();
            linalg_obj = prima_mat.common.linalg_mod();
            memory_obj = prima_mat.common.memory_mod();
            message_obj = prima_mat.common.message_mod();
            prima_mat.common.pintrf_mod();
            powalg_obj = prima_mat.common.powalg_mod();
            ratio_obj = prima_mat.common.ratio_mod();
            redrho_obj = prima_mat.common.redrho_mod();
            shiftbase_obj = prima_mat.common.shiftbase_mod();

            % Solver-specific modules
            geometry_uobyqa_obj = prima_mat.uobyqa.geometry_uobyqa_mod();
            initialize_uobyqa_obj = prima_mat.uobyqa.initialize_uobyqa_mod();
            trustregion_uobyqa_obj = prima_mat.uobyqa.trustregion_uobyqa_mod();
            update_uobyqa_obj = prima_mat.uobyqa.update_uobyqa_mod();


            % Inputs
            % N.B.: INTENT cannot be specified if a dummy procedure is not a POINTER



            % In-outputs
            % X(N)

            % Outputs



            % FHIST(MAXFHIST)
            % XHIST(N, MAXXHIST)

            % Local variables
            solver = "UOBYQA";
            srname = "UOBYQB";
            k = NaN;


            accurate_mod = false;
            adequate_geo = false;
            bad_trstep = false;

            improve_geo = false;
            reduce_rho = false;

            small_trrad = false;


            ximproved = false;

            d = NaN(numel(x), 1);

            delbar = NaN;

            distsq = NaN((numel(x) + 1) * (numel(x) + 2) / 2, 1);

            dnorm_rec = NaN(2, 1); % Powell's implementation: DNORM_REC(3)
            fval = NaN(numel(distsq), 1);
            g = NaN(numel(x), 1);

            h = NaN(numel(x));
            moderr = NaN;
            moderr_rec = NaN(numel(dnorm_rec), 1);
            pq = NaN(numel(distsq) + -1, 1);


            xbase = NaN(numel(x), 1);
            xdrop = NaN(numel(x), 1);
            xpt = NaN(numel(x), numel(distsq));
            pl = NaN;
            trtol = 1.0e-2; % Convergence tolerance of trust-region subproblem solver

            % Sizes.
            n = fix(numel(x));
            npt = (n + 1) * (n + 2) / 2;
            debug_obj.validate(npt > 0, "NPT > 0", srname); % Validate that NPT does not overflow.
            maxxhist = size(xhist, 2);
            maxfhist = fix(numel(fhist));
            maxhist = max(maxxhist, maxfhist);

            % Preconditions.
            if consts_obj.DEBUGGING
                debug_obj.assert(abs(iprint) <= 3, "IPRINT is 0, 1, -1, 2, -2, 3, or -3", srname);
                debug_obj.assert(n >= 1, "N >= 1", srname);
                debug_obj.assert(maxfun >= npt + 1, "MAXFUN >= NPT + 1", srname);
                debug_obj.assert(rhobeg >= rhoend && rhoend > 0, "RHOBEG >= RHOEND > 0", srname);
                debug_obj.assert(all(infnan_obj.is_finite(x), 'all'), "X is finite", srname);
                debug_obj.assert(eta1 >= 0 && eta1 <= eta2 && eta2 < 1, "0 <= ETA1 <= ETA2 < 1", srname);
                debug_obj.assert(gamma1 > 0 && gamma1 < 1 && gamma2 > 1, "0 < GAMMA1 < 1 < GAMMA2", srname);
                debug_obj.assert(maxhist >= 0 && maxhist <= maxfun, "0 <= MAXHIST <= MAXFUN", srname);
                debug_obj.assert(maxfhist * (maxfhist - maxhist) == 0, "SIZE(FHIST) == 0 or MAXHIST", srname);
                debug_obj.assert(size(xhist, 1) == n && maxxhist * (maxxhist - maxhist) == 0, "SIZE(XHIST, 1) == N, SIZE(XHIST, 2) == 0 or MAXHIST", srname);
            end

            %====================%
            % Calculation starts %
            %====================%

            % Initialize XBASE, XPT, FVAL, and KOPT, together with the history and NF.
            [kopt, nf, fhist, fval, xbase, xhist, xpt, subinfo] = initialize_uobyqa_obj.initxf(calfun, iprint, maxfun, ftarget, rhobeg, x, fhist, fval, xbase, xhist, xpt);

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
            x(:) = xbase + xpt(:, kopt);
            f = fval(kopt);

            % Finish the initialization if INITXF completed normally and CALLBACK did not request termination;
            % otherwise, do not proceed, as XPT etc may be uninitialized, leading to errors or exceptions.
            if subinfo == infos_obj.INFO_DFT
                % Initialize the Lagrange polynomials represented by PL. Allocate memory for it first. In
                % general, to make the implementation simple and straightforward, we use automatic arrays rather
                % than allocable ones whenever possible. However, PL is an exception, as its size is O(N^4). If
                % SAFEALLOC fails, an informative error will be raised, which is preferred to a silent or
                % ambiguous failure.
                pl = memory_obj.alloc_rmatrix_sp(pl, npt - 1, npt);
                pl = initialize_uobyqa_obj.initl(xpt, pl);

                % Initialize the quadratic model represented by PQ.
                pq = initialize_uobyqa_obj.initq(fval, xpt, pq);
                if ~(all(infnan_obj.is_finite(pq), 'all'))
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
                    debug_obj.assert(numel(x) == n && ~any(infnan_obj.is_nan(x), 'all'), "SIZE(X) == N, X does not contain NaN", srname);
                    debug_obj.assert(~(infnan_obj.is_nan_sp(f) || infnan_obj.is_posinf(f)), "F is not NaN/+Inf", srname);
                    debug_obj.assert(size(xhist, 1) == n && size(xhist, 2) == maxxhist, "SIZE(XHIST) == [N, MAXXHIST]", srname);
                    debug_obj.assert(~any(infnan_obj.is_nan(xhist(:, 1:min(nf, maxxhist))), 'all'), "XHIST does not contain NaN", srname);
                    % The last calculated X can be Inf (finite + finite can be Inf numerically).
                    debug_obj.assert(numel(fhist) == maxfhist, "SIZE(FHIST) == MAXFHIST", srname);
                    debug_obj.assert(~any(infnan_obj.is_nan(fhist(1:min(nf, maxfhist))) | infnan_obj.is_posinf(fhist(1:min(nf, maxfhist))), 'all'), "FHIST does not contain NaN/+Inf", srname);
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
            shortd = false;

            ratio = -consts_obj.ONE;
            ddmove = -consts_obj.ONE;
            dnorm_rec(:) = consts_obj.REALMAX;
            moderr_rec(:) = consts_obj.REALMAX;
            knew_tr = 0;
            knew_geo = 0;

            % If DELTA <= GAMMA3*RHO after an update, we set DELTA to RHO. GAMMA3 must be less than GAMMA2. The
            % reason is as follows. Imagine a very successful step with DENORM = the un-updated DELTA = RHO.
            % Then TRRAD will update DELTA to GAMMA2*RHO. If GAMMA3 >= GAMMA2, then DELTA will be reset to RHO,
            % which is not reasonable as D is very successful. See paragraph two of Sec. 5.2.5 in
            % T. M. Ragonneau's thesis: "Model-Based Derivative-Free Optimization Methods and Software".
            % According to test on 20230613, for UOBYQA, this Powellful updating scheme of DELTA works better
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
            % UOBYQA never sets IMPROVE_GEO and REDUCE_RHO to TRUE simultaneously.
            for tr = 1:maxtr
                % CLOSE_ITPSET: Are the interpolation points close to XOPT? It affects IMPROVE_GEO, REDUCE_RHO.
                % N.B. (Zaikun 20240331): In Powell's algorithms, CLOSE_ITPSET is defined after XPT is updated
                % according to the trust-region trial step.
                distsq(:) = sum(fortran.power((xpt - fortran.spread(xpt(:, kopt), 'dim', 2, 'ncopies', npt)), 2), 1);
                %%MATLAB: distsq = sum((xpt - xpt(:, kopt)).^2)  % Implicit expansion
                close_itpset = all(distsq <= 4.0 * fortran.power(delta, 2), 'all'); % Powell's NEWUOA code.
                % Below are some alternative definitions of CLOSE_ITPSET.
                % N.B.: The threshold for CLOSE_ITPSET is at least DELBAR, the trust region radius for GEOSTEP.
                % %close_itpset = all(distsq <= 4.0_RP * rho**2)  ! Powell's code.
                % %close_itpset = all(distsq <= max((TWO * delta)**2, (TEN * rho)**2))  ! Powell's BOBYQA code.
                % %close_itpset = all(distsq <= max(delta**2, 4.0_RP * rho**2))  ! Powell's LINCOA code.

                % Generate trust region step D, and also calculate a lower bound on the Hessian of Q.
                g(:) = pq(1:n) + linalg_obj.smat_mul_vec(pq(n + 1:npt - 1), xpt(:, kopt));
                h(:, :) = linalg_obj.vec2smat(pq(n + 1:npt - 1));
                [d, crvmin] = trustregion_uobyqa_obj.trstep(delta, g, h, trtol, d);
                dnorm = min(delta, linalg_obj.p_norm(d));

                % Check whether D is too short to invoke a function evaluation.
                shortd = (dnorm <= consts_obj.HALF * rho); % `<=` works better than `<` in case of underflow.

                % Set QRED to the reduction of the quadratic model when the move D is made from XOPT. QRED
                % should be positive. If it is nonpositive due to rounding errors, we will not take this step.
                qred = -powalg_obj.quadinc_ghv(pq, d, xpt(:, kopt)); % QRED = Q(XOPT) - Q(XOPT + D)
                trfail = (~(qred > 1.0e-6 * fortran.power(rho, 2))); % QRED is tiny/negative or NaN.

                if shortd || trfail
                    % Powell's code does not reduce DELTA as follows. This comes from NEWUOA and works well.
                    delta = consts_obj.TENTH * delta;
                    if delta <= gamma3 * rho
                        delta = rho; % Set DELTA to RHO when it is close to or below.

                    end
                else
                    % Calculate the next value of the objective function.
                    % If X is close to one of the points in the interpolation set, then we do not evaluate the
                    % objective function X, assuming it to have the value at the closest point.
                    x(:) = xbase + (xpt(:, kopt) + d);
                    distsq(:) = reshape(sum(fortran.power((x - (xbase + xpt(:, 1:npt))), 2), 1), [], 1); % Implied do-loop
                    %%MATLAB: distsq = sum((x - (xbase + xpt))**2, 1)  % Implicit expansion
                    k = fix(fortran.minloc(distsq, 'dim', 1));
                    if distsq(k) <= fortran.power((1.0e-4 * rhoend), 2)
                        f = fval(k);
                    else
                        % Evaluate the objective function at X, taking care of possible Inf/NaN values.
                        f = evaluate_obj.evaluatef(calfun, x);
                        nf = nf + 1;
                        % Save X and F into the history.
                        [xhist, fhist] = history_obj.savehist(nf, x, xhist, f, fhist);
                    end

                    % Print a message about the function evaluation according to IPRINT.
                    message_obj.fmsg(solver, "Trust region", iprint, nf, delta, f, x);

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
                    delta = trustregion_uobyqa_obj.trrad(delta, dnorm, eta1, eta2, gamma1, gamma2, ratio);
                    if delta <= gamma3 * rho
                        delta = rho; % Set DELTA to RHO when it is close to or below.

                    end

                    % Is the newly generated X better than current best point?
                    ximproved = (f < fval(kopt));

                    % Set KNEW to the index of the next interpolation point to be deleted.
                    knew_tr = geometry_uobyqa_obj.setdrop_tr(kopt, ximproved, d, pl, rho, xpt);

                    % DDMOVE is norm square of DMOVE in the UOBYQA paper. See Steps 6--7 in Sec. 5 of the paper.
                    ddmove = consts_obj.ZERO;
                    if knew_tr > 0
                        xdrop(:) = xpt(:, knew_tr);
                        % Update PL, PQ, XPT, FVAL, and KOPT so that XPT(:, KNEW_TR) becomes XOPT + D.
                        [kopt, fval, pl, pq, xpt] = update_uobyqa_obj.update(knew_tr, d, f, moderr, kopt, fval, pl, pq, xpt);
                        if ~(all(infnan_obj.is_finite(pq), 'all'))
                            info = infos_obj.NAN_INF_MODEL;
                            break
                        end
                        ddmove = sum(fortran.power((xdrop - xpt(:, kopt)), 2), 'all'); % KOPT is updated.

                    end

                    % Check whether to exit
                    subinfo = checkexit_obj.checkexit_unc(maxfun, nf, f, ftarget, x);
                    if subinfo ~= infos_obj.INFO_DFT
                        info = subinfo;
                        break
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
                accurate_mod = all(abs(moderr_rec) <= 0.125 * crvmin * fortran.power(rho, 2), 'all') && all(dnorm_rec <= rho, 'all');
                % ADEQUATE_GEO: Is the geometry of the interpolation set "adequate"?
                adequate_geo = (shortd && accurate_mod) || close_itpset;
                % SMALL_TRRAD: Is the trust-region radius small? This indicator seems not impactive in practice.
                small_trrad = (max(delta, dnorm) <= rho); % Behaves the same as Powell's version.
                %small_trrad = (dnorm <= rho)  ! Powell's code.

                % Comments on ACCURATE_MOD:
                % 1. ACCURATE_MOD is needed only when SHORTD is TRUE.
                % 2. In Powell's UOBYQA code, ACCURATE_MOD is defined according to (28), (37), and (38) in the
                % UOBYQA paper (see also (32) of Powell 2001: "On the Lagrange functions of quadratic models
                % that are defined by interpolation"). As elaborated in Sec. 3 of the paper (also Sec. 4 of
                % Powell 2001), the idea is to test whether the current model is sufficiently accurate by
                % checking whether the interpolation error bound in (28) is (sufficiently) small. If the bound
                % is small, then set ACCURATE_MOD to TRUE. Otherwise, it identifies a "bad" interpolation point
                % that makes a significant contribution to the bound, with a preference to the interpolation
                % points that are far away from the current trust-region center. Such a point will be replaced
                % with a new point obtained by the geometry step. If all the interpolation points are close
                % enough to the trust-region center, then they are all considered to be good.
                % 3. Our implementation defines ACCURATE_MOD by a method from NEWUOA and BOBYQA, which is also
                % reflected in LINCOA. It sets ACCURATE_MOD to TRUE if recent model errors and step lengths are
                % all small. In addition, it identifies a "bad" interpolation point by simply taking the
                % farthest point from the current trust region center, unless they are all close enough to the
                % center. This implementation is much simpler and less costly in terms of flops yet it performs
                % almost the same as Powell's original implementation.

                % Powell's original definition of IMPROVE_GEO and REDUCE_RHO:
                % %bad_trstep = (shortd .or. knew_tr == 0 .or. (ratio <= 0 .and. dnorm <= 2.0_RP*rho .and. ddmove <= 4.0_RP * rho**2))
                % %improve_geo = bad_trstep .and. .not. (shortd .and. accurate_mod) .and. .not. close_itpset
                % %reduce_rho = bad_trstep .and. dnorm <= rho .and. .not. improve_geo

                % IMPROVE_GEO and REDUCE_RHO are defined as follows.
                % N.B.: If SHORTD is TRUE at the very first iteration, then REDUCE_RHO will be set to TRUE.
                % Powell's code does not have TRFAIL in BAD_TRSTEP; it terminates if TRFAIL is TRUE.

                % BAD_TRSTEP (for IMPROVE_GEO): Is the last trust-region step bad? For UOBYQA, it is CRUCIAL to
                % include DMOVE <= 4.0_RP*RHO**2 in the definition of BAD_TRSTEP for IMPROVE_GEO.
                bad_trstep = (shortd || trfail || (ratio <= eta1 && ddmove <= 4.0 * fortran.power(delta, 2)) || knew_tr == 0);
                %bad_trstep = (shortd .or. trfail .or. ratio <= eta1 .or. knew_tr == 0)  ! Works poorly!
                improve_geo = bad_trstep && ~adequate_geo;
                % BAD_TRSTEP (for REDUCE_RHO): Is the last trust-region step bad?
                bad_trstep = (shortd || trfail || ratio <= 0 || knew_tr == 0); % Performs better than the one below from Powell.
                %bad_trstep = (shortd .or. trfail .or. (ratio <= 0 .and. ddmove <= 4.0_RP * delta**2) .or. knew_tr == 0)
                reduce_rho = bad_trstep && adequate_geo && small_trrad;

                % Equivalently, REDUCE_RHO can be set as follows. It shows that REDUCE_RHO is TRUE in two cases.
                % %bad_trstep = (shortd .or. trfail .or. (ratio <= 0 .and. ddmove <= 4.0_RP * delta**2) .or. knew_tr == 0)
                % %reduce_rho = (shortd .and. accurate_mod) .or. (bad_trstep .and. close_itpset .and. small_trrad)

                % With REDUCE_RHO properly defined, we can also set IMPROVE_GEO as follows.
                % %bad_trstep = (shortd .or. trfail .or. (ratio <= eta1 .and. ddmove <= 4.0_RP * delta**2) .or. knew_tr == 0)
                % %improve_geo = bad_trstep .and. (.not. reduce_rho) .and. (.not. close_itpset)

                % With IMPROVE_GEO properly defined, we can also set REDUCE_RHO as follows.
                % %bad_trstep = (shortd .or. trfail .or. (ratio <= 0 .and. ddmove <= 4.0_RP * delta**2) .or. knew_tr == 0)
                % %reduce_rho = bad_trstep .and. (.not. improve_geo) .and. small_trrad

                % UOBYQA never sets IMPROVE_GEO and REDUCE_RHO to TRUE simultaneously.
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
                    distsq(:) = sum(fortran.power((xpt - fortran.spread(xpt(:, kopt), 'dim', 2, 'ncopies', npt)), 2), 1);
                    %%MATLAB: distsq = sum((xpt - xpt(:, kopt)).^2)  % Implicit expansion
                    knew_geo = fix(fortran.maxloc(distsq, 'dim', 1));

                    % DELBAR is the trust-region radius for the geometry improvement subproblem.
                    % Powell's UOBYQA code sets DELBAR = RHO, but NEWUOA/BOBYQA/LINCOA all take DELTA and/or
                    % DISTSQ into consideration.
                    delbar = rho; % Powell's code
                    %delbar = max(min(TENTH * sqrt(maxval(distsq)), HALF * delta), rho)  ! Powell's NEWUOA code
                    %delbar = max(TENTH * delta, rho)  ! Powell's LINCOA code
                    %delbar = max(min(TENTH * sqrt(maxval(distsq)), delta), rho)  ! Powell's BOBYQA code

                    d(:) = geometry_uobyqa_obj.geostep(knew_geo, kopt, delbar, pl, xpt);

                    % Calculate the next value of the objective function.
                    % If X is close to one of the points in the interpolation set, then we do not evaluate the
                    % objective function X, assuming it to have the value at the closest point.
                    x(:) = xbase + (xpt(:, kopt) + d);
                    distsq(:) = reshape(sum(fortran.power((x - (xbase + xpt(:, 1:npt))), 2), 1), [], 1); % Implied do-loop
                    %%MATLAB: distsq = sum((x - (xbase + xpt))**2, 1)  % Implicit expansion
                    k = fix(fortran.minloc(distsq, 'dim', 1));
                    if distsq(k) <= fortran.power((1.0e-4 * rhoend), 2)
                        f = fval(k);
                    else
                        % Evaluate the objective function at X, taking care of possible Inf/NaN values.
                        f = evaluate_obj.evaluatef(calfun, x);
                        nf = nf + 1;
                        % Save X and F into the history.
                        [xhist, fhist] = history_obj.savehist(nf, x, xhist, f, fhist);
                    end

                    % Print a message about the function evaluation according to IPRINT.
                    message_obj.fmsg(solver, "Geometry", iprint, nf, delbar, f, x);

                    % Update DNORM_REC and MODERR_REC.
                    % DNORM_REC records the DNORM of the recent function evaluations with the current RHO.
                    dnorm = min(delbar, linalg_obj.p_norm(d)); % In theory, DNORM = DELBAR in this case.
                    dnorm_rec(:) = [reshape(dnorm_rec(2:numel(dnorm_rec)), [], 1); dnorm];
                    % MODERR is the error of the current model in predicting the change in F due to D.
                    % MODERR_REC records the prediction errors of the recent models with the current RHO.
                    moderr = f - fval(kopt) - powalg_obj.quadinc_ghv(pq, d, xpt(:, kopt)); % QUADINC = Q(XOPT + D) - Q(XOPT)
                    moderr_rec(:) = [reshape(moderr_rec(2:numel(moderr_rec)), [], 1); moderr];

                    % Update PL, PQ, XPT, FVAL, and KOPT so that XPT(:, KNEW_GEO) becomes XOPT + D.
                    [kopt, fval, pl, pq, xpt] = update_uobyqa_obj.update(knew_geo, d, f, moderr, kopt, fval, pl, pq, xpt);
                    if ~(all(infnan_obj.is_finite(pq), 'all'))
                        info = infos_obj.NAN_INF_MODEL;
                        break
                    end

                    % Check whether to exit
                    subinfo = checkexit_obj.checkexit_unc(maxfun, nf, f, ftarget, x);
                    if subinfo ~= infos_obj.INFO_DFT
                        info = subinfo;
                        break
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

                % Shifting XBASE to the best point so far, and make the corresponding changes to the gradients
                % of the Lagrange functions and the quadratic model. Powell's implementation does this each time
                % after RHO is reduced. Our implementation aligns with NEWUOA/BOBYQA/LINCOA.
                if sum(fortran.power(xpt(:, kopt), 2), 'all') >= 1000.0 * fortran.power(delta, 2)
                    [pl, pq, xbase, xpt] = shiftbase_obj.shiftbase_qint(kopt, pl, pq, xbase, xpt);
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

            % Deallocate PL. F2003 automatically deallocate local ALLOCATABLE variables at exit, yet we prefer
            % to deallocate them immediately when they finish their jobs.


            % Return from the calculation, after trying the Newton-Raphson step if it has not been tried yet.
            % Ensure that D has not been updated after SHORTD == TRUE occurred, or the code below is incorrect.
            x(:) = xbase + (xpt(:, kopt) + d);
            if info == infos_obj.SMALL_TR_RADIUS && shortd && linalg_obj.p_norm(x - (xbase + xpt(:, kopt))) > consts_obj.TENTH * rhoend && nf < maxfun
                f = evaluate_obj.evaluatef(calfun, x);
                nf = nf + 1;
                % Save X, F into the history.
                [xhist, fhist] = history_obj.savehist(nf, x, xhist, f, fhist);
                % Print a message about the function evaluation according to IPRINT.
                % Zaikun 20230512: DELTA has been updated. RHO is only indicative here. TO BE IMPROVED.
                message_obj.fmsg(solver, "Trust region", iprint, nf, rho, f, x);
                if f < fval(kopt)
                    xpt(:, kopt) = xpt(:, kopt) + d;
                    fval(kopt) = f;
                end
            end

            % Choose the [X, F] to return.
            x(:) = xbase + xpt(:, kopt);
            f = fval(kopt);

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
                debug_obj.assert(numel(x) == n && ~any(infnan_obj.is_nan(x), 'all'), "SIZE(X) == N, X does not contain NaN", srname);
                debug_obj.assert(~(infnan_obj.is_nan_sp(f) || infnan_obj.is_posinf(f)), "F is not NaN/+Inf", srname);
                debug_obj.assert(size(xhist, 1) == n && size(xhist, 2) == maxxhist, "SIZE(XHIST) == [N, MAXXHIST]", srname);
                debug_obj.assert(~any(infnan_obj.is_nan(xhist(:, 1:min(nf, maxxhist))), 'all'), "XHIST does not contain NaN", srname);
                % The last calculated X can be Inf (finite + finite can be Inf numerically).
                debug_obj.assert(numel(fhist) == maxfhist, "SIZE(FHIST) == MAXFHIST", srname);
                debug_obj.assert(~any(infnan_obj.is_nan(fhist(1:min(nf, maxfhist))) | infnan_obj.is_posinf(fhist(1:min(nf, maxfhist))), 'all'), "FHIST does not contain NaN/+Inf", srname);
                debug_obj.assert(~any(fhist(1:min(nf, maxfhist)) < f, 'all'), "F is the smallest in FHIST", srname);
            end

        end

    end
end