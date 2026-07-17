%TODO:
% 1. Check whether it is possible to change the definition of RESCON, RESNEW, RESTMP, RESACT so that
% we do not need to encode information into their signs.
%
classdef lincob_mod
    %--------------------------------------------------------------------------------------------------%
    % This module performs the major calculations of LINCOA.
    %
    % Coded by Zaikun ZHANG (www.zhangzk.net) based on Powell's code and the paper
    %
    % M. J. D. Powell, On fast trust region methods for quadratic models with linear constraints,
    % Math. Program. Comput., 7:237--267, 2015
    %
    % Dedicated to the late Professor M. J. D. Powell FRS (1936--2015).
    %
    % Started: February 2022
    %
    % Last Modified: Wed 08 Apr 2026 06:38:26 PM CST
    %--------------------------------------------------------------------------------------------------%

    methods
        function [x, nf, chist, cstrv, f, fhist, xhist, info] = lincob(~, calfun, iprint, maxfilt, maxfun, npt, Aeq, Aineq, amat, beq, bineq, bvec, ctol, cweight, eta1, eta2, ftarget, gamma1, gamma2, rhobeg, rhoend, xl, xu, x, chist, fhist, xhist, varargin)
            %--------------------------------------------------------------------------------------------------%
            % This subroutine performs the actual calculations of LINCOA.
            %
            % The arguments IPRINT, MAXFILT, MAXFUN, MAXHIST, NPT, AEQ, AINEQ, BEQ, BINEQ, CTOL, CWEIGHT, ETA1,
            % ETA2, FTARGET, GAMMA1, GAMMA2, RHOBEG, RHOEND, X, NF, F, XHIST, FHIST, CHIST, CSTRV and INFO are
            % identical to the corresponding arguments in subroutine LINCOA.
            % AMAT is a matrix whose columns are the constraint gradients, scaled so that they have unit length.
            % BVEC contains on entry the right hand sides of the constraints, scaled as above.
            % XBASE holds a shift of origin that should reduce the contributions from rounding errors to values
            %   of the model and Lagrange functions.
            % XOPT is the displacement from XBASE of the feasible vector of variables that provides the least
            %   calculated F so far, this vector being the current trust region centre. FOPT = F(XOPT + XBASE).
            %   However, we do not save XOPT and FOPT explicitly, because XOPT = XPT(:, KOPT) and
            %   FOPT = FVAL(KOPT), which is explained below.
            % [XPT, FVAL, KOPT] describes the interpolation set:
            % XPT contains the interpolation points relative to XBASE, each COLUMN for a point; FVAL holds the
            %   values of F at the interpolation points; KOPT is the index of XOPT in XPT.
            % [GOPT, HQ, PQ] describes the quadratic model: GOPT will hold the gradient of the quadratic model
            %   at XBASE + XOPT; HQ will hold the explicit second order derivatives of the quadratic model; PQ
            %   will contain the parameters of the implicit second order derivatives of the quadratic model.
            % [BMAT, ZMAT, IDZ] describes the matrix H in the NEWUOA paper (eq. 3.12), which is the inverse of
            %   the coefficient matrix of the KKT system for the least-Frobenius norm interpolation problem:
            %   ZMAT will hold a factorization of the leading NPT*NPT submatrix of H, the factorization being
            %   ZMAT*Diag(DZ)*ZMAT^T with DZ(1:IDZ-1)=-1, DZ(IDZ:NPT-N-1)=1. BMAT will hold the last N ROWs of H
            %   except for the (NPT+1)th column. Note that the (NPT + 1)th row and column of H are not saved as
            %   they are unnecessary for the calculation.
            % D is reserved for trial steps from XOPT. It is chosen by subroutine TRSTEP or GEOSTEP. Usually
            %   XBASE + XOPT + D is the vector of variables for the next call of CALFUN.
            % IACT is an integer array for the indices of the active constraints.
            % RESCON holds information about the constraint residuals at the current trust region center XOPT.
            %   1. If if B(J) - AMAT(:, J)^T*XOPT <= DELTA, then RESCON(J) = B(J) - AMAT(:, J)^T*XOPT. Note that
            %   RESCON >= 0 in this case, because the algorithm keeps XOPT to be feasible.
            %   2. Otherwise, RESCON(J) is a negative value that B(J) - AMAT(:,J)^T*XOPT >= |RESCON(J)| >= DELTA.
            %   RESCON can be updated without calculating the constraints that are far from being active, so
            %   that we only need to evaluate the constraints that are nearly active.
            % QFAC is the orthogonal part of the QR factorization of the matrix of active constraint gradients,
            %   these gradients being ordered in accordance with IACT. When NACT is less than N, columns are
            %   appended to QFAC to complete an N by N orthogonal matrix, which is important for keeping
            %   calculated steps sufficiently close to the boundaries of the active constraints.
            % RFAC is the upper triangular part of this QR factorization.
            %--------------------------------------------------------------------------------------------------%

            % Generic models
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
            selectx_obj = prima_mat.common.selectx_mod();
            shiftbase_obj = prima_mat.common.shiftbase_mod();

            % Solver-specific modules
            geometry_lincoa_obj = prima_mat.lincoa.geometry_lincoa_mod();
            initialize_lincoa_obj = prima_mat.lincoa.initialize_lincoa_mod();
            trustregion_lincoa_obj = prima_mat.lincoa.trustregion_lincoa_mod();
            update_lincoa_obj = prima_mat.lincoa.update_lincoa_mod();


            % Inputs
            % N.B.: INTENT cannot be specified if a dummy procedure is not a POINTER



            % Aeq(Meq, N)
            % Aineq(Mineq, N)
            % AMAT(N, M)
            % Beq(Meq)
            % Bineq(Mineq)
            % BVEC(M)



            % In-outputs
            % X(N)

            % Outputs


            % CHIST(MAXCHIST)


            % FHIST(MAXFHIST)
            % XHIST(N, MAXXHIST)

            % Local variables
            solver = "LINCOA";
            srname = "LINCOB";
            iact = NaN(numel(bvec), 1);
            idz = NaN;
            ij = NaN(2, max(0, fix(npt - 2 * numel(x) - 1)));


            accurate_mod = false;
            adequate_geo = false;
            bad_trstep = false;
            close_itpset = false;
            evaluated = false(npt, 1);
            feasible = false;
            improve_geo = false;
            qalt_better = false(3, 1);
            reduce_rho = false;

            small_trrad = false;


            ximproved = false;
            b = NaN(numel(bvec), 1);
            bmat = NaN(numel(x), npt + numel(x));
            cfilt = NaN(maxfilt, 1);
            constr = NaN(nnz(xl > -consts_obj.BOUNDMAX) + nnz(xu < consts_obj.BOUNDMAX) + 2 * numel(beq) + numel(bineq), 1);
            constr_leq = NaN(numel(beq), 1);
            cval = NaN(npt, 1);
            d = NaN(numel(x), 1);
            delbar = NaN;

            distsq = NaN(npt, 1);
            dnorm = NaN;
            dnorm_rec = NaN(3, 1); % Powell's implementation: DNORM_REC(5)
            ffilt = NaN(maxfilt, 1);
            fval = NaN(npt, 1);
            galt = NaN(numel(x), 1);

            gopt = NaN(numel(x), 1);
            hq = NaN(numel(x));
            moderr = NaN;
            moderr_alt = NaN;
            pq = NaN(npt, 1);
            pqalt = NaN(npt, 1);
            qfac = NaN(numel(x));


            rescon = NaN(numel(bvec), 1);
            rfac = NaN(numel(x));

            xbase = NaN(numel(x), 1);
            xdrop = NaN(numel(x), 1);
            xfilt = NaN(numel(x), maxfilt);
            xosav = NaN(numel(x), 1);
            xpt = NaN(numel(x), npt);
            zmat = NaN(npt, npt - numel(x) + -1);
            trtol = 1.0e-2; % Convergence tolerance of trust-region subproblem solver

            % Sizes.
            m = fix(numel(bvec));
            n = fix(numel(x));
            maxxhist = size(xhist, 2);
            maxfhist = fix(numel(fhist));
            maxchist = fix(numel(chist));
            maxhist = fix(max(maxxhist, max(maxfhist, maxchist)));

            % Preconditions
            if consts_obj.DEBUGGING
                debug_obj.assert(abs(iprint) <= 3, "IPRINT is 0, 1, -1, 2, -2, 3, or -3", srname);
                debug_obj.assert(m >= 0, "M >= 0", srname);
                debug_obj.assert(n >= 1, "N >= 1", srname);
                debug_obj.assert(npt >= n + 2, "NPT >= N+2", srname);
                debug_obj.assert(maxfun >= npt + 1, "MAXFUN >= NPT+1", srname);
                debug_obj.assert(size(Aeq, 1) == numel(beq) && size(Aeq, 2) == n, "SIZE(Aeq) == [SIZE(Beq), N]", srname);
                debug_obj.assert(size(Aineq, 1) == numel(bineq) && size(Aineq, 2) == n, "SIZE(Aineq) == [SIZE(Bineq), N]", srname);
                debug_obj.assert(size(amat, 1) == n && size(amat, 2) == m, "SIZE(AMAT) == [N, M]", srname);
                debug_obj.assert(eta1 >= 0 && eta1 <= eta2 && eta2 < 1, "0 <= ETA1 <= ETA2 < 1", srname);
                debug_obj.assert(gamma1 > 0 && gamma1 < 1 && gamma2 > 1, "0 < GAMMA1 < 1 < GAMMA2", srname);
                debug_obj.assert(rhobeg >= rhoend && rhoend > 0, "RHOBEG >= RHOEND > 0", srname);
                debug_obj.assert(all(infnan_obj.is_finite(x), 'all'), "X is finite", srname);
                debug_obj.assert(numel(xl) == n && numel(xu) == n, "SIZE(XL) == N == SIZE(XU)", srname);
                debug_obj.assert(maxfilt >= min(consts_obj.MIN_MAXFILT, maxfun) && maxfilt <= maxfun, "MIN(MIN_MAXFILT, MAXFUN) <= MAXFILT <= MAXFUN", srname);
                debug_obj.assert(maxhist >= 0 && maxhist <= maxfun, "0 <= MAXHIST <= MAXFUN", srname);
                debug_obj.assert(size(xhist, 1) == n && maxxhist * (maxxhist - maxhist) == 0, "SIZE(XHIST, 1) == N, SIZE(XHIST, 2) == 0 or MAXHIST", srname);
                debug_obj.assert(maxfhist * (maxfhist - maxhist) == 0, "SIZE(FHIST) == 0 or MAXHIST", srname);
                debug_obj.assert(maxchist * (maxchist - maxhist) == 0, "SIZE(CHIST) == 0 or MAXHIST", srname);
            end

            %====================%
            % Calculation starts %
            %====================%

            % IXL and IXU are the indices of the nontrivial lower and upper bounds, respectively.
            memory_obj.alloc_ivector(fix(nnz(xl > -consts_obj.BOUNDMAX))); % Removable in F2003.
            memory_obj.alloc_ivector(fix(nnz(xu < consts_obj.BOUNDMAX))); % Removable in F2003.
            ixl = linalg_obj.trueloc(xl > -consts_obj.BOUNDMAX);
            ixu = linalg_obj.trueloc(xu < consts_obj.BOUNDMAX);

            % Initialize B, XBASE, XPT, FVAL, CVAL, and KOPT, together with the history, NF, IJ, and EVALUATED.
            b(:) = bvec;
            [b, ij, kopt, nf, chist, cval, fhist, fval, xbase, xhist, xpt, evaluated, subinfo] = initialize_lincoa_obj.initxf(calfun, iprint, maxfun, Aeq, Aineq, amat, beq, bineq, ctol, ftarget, rhobeg, xl, xu, x, b, ij, chist, cval, fhist, fval, xbase, xhist, xpt, evaluated);

            % Report the current best value, and check if user asks for early termination.
            terminate = false;
            ipObj = inputParser();
            addParameter(ipObj, 'callback_fcn', struct());
            parse(ipObj, varargin{:});
            callback_fcn = ipObj.Results.callback_fcn;
            if ~ismember('callback_fcn', ipObj.UsingDefaults)
                terminate = callback_fcn(xbase + xpt(:, kopt), fval(kopt), nf, 0, 'cstrv', cval(kopt));
                if terminate
                    subinfo = infos_obj.CALLBACK_TERMINATE;
                end
            end

            % Initialize X, F, CONSTR, and CSTRV according to KOPT.
            % N.B.: We must set CONSTR and CSTRV. Otherwise, if REDUCE_RHO is TRUE after the very first
            % iteration due to SHORTD, then RHOMSG will be called with CONSTR and CSTRV uninitialized.
            x(:) = xbase + xpt(:, kopt);
            f = fval(kopt);
            constr_leq(:) = linalg_obj.matprod21(Aeq, x) - beq;
            constr(:) = [reshape(xl(ixl) - x(ixl), [], 1); reshape(x(ixu) - xu(ixu), [], 1); reshape(-constr_leq, [], 1); reshape(constr_leq, [], 1); reshape(linalg_obj.matprod21(Aineq, x) - bineq, [], 1)];
            cstrv = linalg_obj.maximum1([consts_obj.ZERO; reshape(constr, [], 1)]);

            % Initialize the filter, including XFILT, FFILT, CONFILT, CFILT, and NFILT.
            % N.B.: The filter is used only when selecting which iterate to return. It does not interfere with
            % the iterations. LINCOA is NOT a filter method but a trust-region method. All the trust-region
            % iterates are supposed to be feasible, but can be infeasible due to rounding errors; the
            % geometry-improving iterates are not necessarily feasible. Powell's implementation does not use a
            % filter to select the iterate, possibly returning a suboptimal iterate.
            nfilt = 0;
            for k = 1:npt
                if evaluated(k)
                    [nfilt, cfilt, ffilt, xfilt] = selectx_obj.savefilt(cval(k), ctol, cweight, fval(k), xbase + xpt(:, k), nfilt, cfilt, ffilt, xfilt);
                end
            end

            % Finish the initialization if INITXF completed normally and CALLBACK did not request termination;
            % otherwise, do not proceed, as XPT etc may be uninitialized, leading to errors or exceptions.
            if subinfo == infos_obj.INFO_DFT
                % Initialize [BMAT, ZMAT, IDZ], representing inverse of KKT matrix of the interpolation system.
                [idz, bmat, zmat] = initialize_lincoa_obj.inith(ij, xpt, bmat, zmat);

                % Initialize the quadratic represented by [GOPT, HQ, PQ], so that its gradient at XBASE+XOPT is
                % GOPT; its Hessian is HQ + sum_{K=1}^NPT PQ(K)*XPT(:, K)*XPT(:, K)'.
                hq = repmat(consts_obj.ZERO, size(hq));
                pq(:) = powalg_obj.omega_mul(idz, zmat, fval);
                gopt(:) = linalg_obj.matprod21(bmat(:, 1:npt), fval) + powalg_obj.hess_mul(xpt(:, kopt), xpt, pq);
                pqalt(:) = pq;
                galt(:) = gopt;
                if ~(all(infnan_obj.is_finite(gopt), 'all') && all(infnan_obj.is_finite(hq), 'all') && all(infnan_obj.is_finite(pq), 'all'))
                    subinfo = infos_obj.NAN_INF_MODEL;
                end
            end

            % Check whether to return due to abnormal cases that may occur during the initialization.
            if subinfo ~= infos_obj.INFO_DFT
                info = subinfo;
                % Return the best calculated values of the variables. If CTOL > 0, the KOPT decided by SELECTX
                % may not be the same as the one by INITXF.
                kopt = selectx_obj.selectx(ffilt(1:nfilt), cfilt(1:nfilt), cweight, ctol);
                x(:) = xfilt(:, kopt);
                f = ffilt(kopt);
                constr_leq(:) = linalg_obj.matprod21(Aeq, x) - beq;
                constr(:) = [reshape(xl(ixl) - x(ixl), [], 1); reshape(x(ixu) - xu(ixu), [], 1); reshape(-constr_leq, [], 1); reshape(constr_leq, [], 1); reshape(linalg_obj.matprod21(Aineq, x) - bineq, [], 1)];
                cstrv = linalg_obj.maximum1([consts_obj.ZERO; reshape(constr, [], 1)]);
                message_obj.retmsg(solver, info, iprint, nf, f, x, 'cstrv', cstrv, 'constr', constr);
                % Arrange CHIST, FHIST, and XHIST so that they are in the chronological order.
                [xhist, fhist, chist] = history_obj.rangehist(nf, xhist, fhist, 'chist', chist);
                % Postconditions
                if consts_obj.DEBUGGING
                    debug_obj.assert(nf <= maxfun, "NF <= MAXFUN", srname);
                    debug_obj.assert(numel(x) == n && ~any(infnan_obj.is_nan_sp(x), 'all'), "SIZE(X) == N, X does not contain NaN", srname);
                    debug_obj.assert(~(infnan_obj.is_nan_sp(f) || infnan_obj.is_posinf(f)), "F is not NaN/+Inf", srname);
                    debug_obj.assert(size(xhist, 1) == n && size(xhist, 2) == maxxhist, "SIZE(XHIST) == [N, MAXXHIST]", srname);
                    debug_obj.assert(~any(infnan_obj.is_nan_sp(xhist(:, 1:min(nf, maxxhist))), 'all'), "XHIST does not contain NaN", srname);
                    % The last calculated X can be Inf (finite + finite can be Inf numerically).
                    debug_obj.assert(numel(fhist) == maxfhist, "SIZE(FHIST) == MAXFHIST", srname);
                    debug_obj.assert(~any(infnan_obj.is_nan_sp(fhist(1:min(nf, maxfhist))) | infnan_obj.is_posinf(fhist(1:min(nf, maxfhist))), 'all'), "FHIST does not contain NaN/+Inf", srname);
                    debug_obj.assert(numel(chist) == maxchist, "SIZE(CHIST) == MAXCHIST", srname);
                    debug_obj.assert(~any(chist(1:min(nf, maxchist)) < 0 | infnan_obj.is_nan_sp(chist(1:min(nf, maxchist))) | infnan_obj.is_posinf(chist(1:min(nf, maxchist))), 'all'), "CHIST does not contain negative values or NaN/+Inf", srname);
                    nhist = min([nf, maxfhist, maxchist], [], 'all');
                    debug_obj.assert(~any(selectx_obj.isbetter10(fhist(1:nhist), chist(1:nhist), f, cstrv, ctol), 'all'), "No point in the history is better than X", srname);
                end
                return
            end

            % Initialize RESCON.
            rescon(:) = max(b - linalg_obj.matprod12(xpt(:, kopt), amat), consts_obj.ZERO);
            rescon(linalg_obj.trueloc(rescon >= rhobeg)) = -rescon(linalg_obj.trueloc(rescon >= rhobeg));
            %%MATLAB: rescon(rescon >= rhobeg) = -rescon(rescon >= rhobeg)

            % Set some more initial values.
            % We must initialize RATIO. Otherwise, when SHORTD = TRUE, compilers may raise a run-time error that
            % RATIO is undefined. But its value will not be used: when SHORTD = FALSE, its value will be
            % overwritten; when SHORTD = TRUE, its value is used only in BAD_TRSTEP, which is TRUE regardless of
            % RATIO. Similar for KNEW_TR.
            % No need to initialize SHORTD unless MAXTR < 1, but some compilers may complain if we do not do it.
            rho = rhobeg;
            delta = rho;
            ratio = -consts_obj.ONE;
            dnorm_rec(:) = consts_obj.REALMAX;
            shortd = false;

            qalt_better(:) = false;
            knew_tr = 0;
            knew_geo = 0;
            qfac(:, :) = linalg_obj.eye1(n);
            rfac = repmat(consts_obj.ZERO, size(rfac));
            nact = 0;
            iact(:) = linalg_obj.linspace_i(1, m, m);

            % If DELTA <= GAMMA3*RHO after an update, we set DELTA to RHO. GAMMA3 must be less than GAMMA2. The
            % reason is as follows. Imagine a very successful step with DENORM = the un-updated DELTA = RHO.
            % Then TRRAD will update DELTA to GAMMA2*RHO. If GAMMA3 >= GAMMA2, then DELTA will be reset to RHO,
            % which is not reasonable as D is very successful. See paragraph two of Sec. 5.2.5 in
            % T. M. Ragonneau's thesis: "Model-Based Derivative-Free Optimization Methods and Software".
            % According to test on 20230613, for LINCOA, this Powellful updating scheme of DELTA works evidently
            % better than setting directly DELTA = MAX(NEW_DELTA, RHO).
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
            % LINCOA never sets IMPROVE_GEO and REDUCE_RHO to TRUE simultaneously.
            for tr = 1:maxtr
                % Generate the next trust region step D by calling TRSTEP. Note that D is feasible.
                [iact, nact, qfac, rfac, d, ngetact] = trustregion_lincoa_obj.trstep(amat, delta, gopt, hq, pq, rescon, trtol, xpt, iact, nact, qfac, rfac, d);
                dnorm = min(delta, linalg_obj.p_norm(d));

                % A trust region step is applied whenever its length is at least 0.5*DELTA. It is also
                % applied if its length is at least 0.1999*DELTA and if a line search of TRSTEP has caused a
                % change to the active set, indicated by NGETACT >= 2 (note that NGETACT is at least 1).
                % Otherwise, the trust region step is considered too short to try.
                % N.B. The magic number 0.1999 seems to be related to the fact that a linear constraint is
                % considered nearly active if the point under consideration is within 0.2*DELTA to the boundary
                % of the constraint. See the subroutine GETACT and Section 3 of Powell (2015) for more details.
                % `<=` works better than `<` in case of underflow.
                shortd = ((dnorm <= consts_obj.HALF * delta && ngetact < 2) || dnorm <= 0.1999 * delta);
                %------------------------------------------------------------------------------------------%
                % The SHORTD defined above needs NGETACT, which relies on Powell's trust region subproblem
                % solver. If a different subproblem solver is used, we can take the following SHORTD adopted
                % from UOBYQA, NEWUOA and BOBYQA.
                % %SHORTD = (DNORM < HALF * RHO)
                %------------------------------------------------------------------------------------------%

                % DNORM_REC records the DNORM of recent trust-region iterations. It will be used to decide
                % whether we should improve the geometry of the interpolation set or reduce RHO when SHORTD
                % is TRUE. Note that it does not record the geometry steps.
                dnorm_rec(:) = [reshape(dnorm_rec(2:numel(dnorm_rec)), [], 1); dnorm];

                % In some cases, we reset DNORM_REC to REALMAX. This indicates a preference of improving the
                % geometry of the interpolation set to reducing RHO in the subsequent three or more iterations.
                % This is important for the performance of LINCOA.
                % Zaikun 20230609: This does not exist in NEWUOA/BOBYQA/UOBYQA. Try it!
                if delta > rho || ~shortd
                    % Another possibility: IF (DELTA > RHO) THEN
                    dnorm_rec(:) = consts_obj.REALMAX;
                end

                % Set QRED to the reduction of the quadratic model when the move D is made from XOPT. QRED
                % should be positive. If it is nonpositive due to rounding errors, we will not take this step.
                qred = -powalg_obj.quadinc_d0(d, xpt, gopt, pq, 'hq', hq); % QRED = Q(XOPT) - Q(XOPT + D)
                trfail = (~(qred > 1.0e-6 * fortran.power(rho, 2))); % QRED is tiny/negative or NaN.

                if shortd || trfail
                    % In this case, do nothing but reducing DELTA. Afterward, DELTA < DNORM may occur.
                    % N.B.: 1. This value of DELTA will be discarded if REDUCE_RHO turns out TRUE later.
                    % 2. Powell's code does not shrink DELTA when TRFAIL is TRUE (i.e., when VQUAD >= 0 in
                    % Powell's code, where VQUAD = -QRED). Consequently, the algorithm may be stuck in an
                    % infinite cycling, because both REDUCE_RHO and IMPROVE_GEO may end up with FALSE in this
                    % case, which did happen in tests.
                    % 3. The factor HALF works better than TENTH (used in NEWUOA/BOBYQA), 0.2, and 0.7.
                    delta = consts_obj.HALF * delta;
                    if delta <= gamma3 * rho
                        delta = rho; % Set DELTA to RHO when it is close to or below.

                    end
                else
                    % Calculate the next value of the objective function.
                    x(:) = xbase + (xpt(:, kopt) + d);
                    f = evaluate_obj.evaluatef(calfun, x);
                    nf = nf + 1;

                    % Evaluate the constraints. They are used only for printing messages.
                    constr_leq(:) = linalg_obj.matprod21(Aeq, x) - beq;
                    constr(:) = [reshape(xl(ixl) - x(ixl), [], 1); reshape(x(ixu) - xu(ixu), [], 1); reshape(-constr_leq, [], 1); reshape(constr_leq, [], 1); reshape(linalg_obj.matprod21(Aineq, x) - bineq, [], 1)];
                    cstrv = linalg_obj.maximum1([consts_obj.ZERO; reshape(constr, [], 1)]);

                    % Print a message about the function evaluation according to IPRINT.
                    message_obj.fmsg(solver, "Trust region", iprint, nf, delta, f, x, 'cstrv', cstrv, 'constr', constr);
                    % Save X, F, CSTRV into the history.
                    [xhist, fhist, chist] = history_obj.savehist(nf, x, xhist, f, fhist, 'cstrv', cstrv, 'chist', chist);
                    % Save X, F, CSTRV into the filter.
                    [nfilt, cfilt, ffilt, xfilt] = selectx_obj.savefilt(cstrv, ctol, cweight, f, x, nfilt, cfilt, ffilt, xfilt);

                    % Check whether to exit.
                    subinfo = checkexit_obj.checkexit_con(maxfun, nf, cstrv, ctol, f, ftarget, x);
                    if subinfo ~= infos_obj.INFO_DFT
                        info = subinfo;
                        break
                    end

                    % QALT_BETTER is a boolean array indicating whether the recent few (three) alternative
                    % models are more accurate in predicting the function value at XOPT + D.
                    % N.B.: Do NOT change the "<" in the comparison to "<="; otherwise, the result will not be
                    % reasonable if the two values being compared are both ZERO or INF.
                    moderr = f - fval(kopt) + qred;
                    moderr_alt = f - fval(kopt) - powalg_obj.quadinc_d0(d, xpt, galt, pqalt);
                    qalt_better(:) = [reshape(qalt_better(2:numel(qalt_better)), [], 1); abs(moderr_alt) < consts_obj.TENTH * abs(moderr)];

                    % Calculate the reduction ratio by REDRAT, which handles Inf/NaN carefully.
                    ratio = ratio_obj.redrat(fval(kopt) - f, qred, eta1);

                    % Update DELTA. After this, DELTA < DNORM may hold.
                    % The new DELTA lies in [GAMMA1*DNORM, GAMMA2*DNORM].
                    delta = trustregion_lincoa_obj.trrad(delta, dnorm, eta1, eta2, gamma1, gamma2, ratio);
                    if delta <= gamma3 * rho
                        delta = rho; % Set DELTA to RHO when it is close to or below.

                    end

                    % Is the newly generated X better than current best point?
                    ximproved = (f < fval(kopt));

                    % Set KNEW_TR to the index of the interpolation point to be replaced with XNEW = XOPT + D.
                    % KNEW_TR will ensure that the geometry of XPT is "good enough" after the replacement.
                    knew_tr = geometry_lincoa_obj.setdrop_tr(idz, kopt, ximproved, bmat, d, delta, rho, xpt, zmat);
                    if knew_tr > 0
                        % Update [BMAT, ZMAT, IDZ] (represents H in the NEWUOA paper), [XPT, FVAL, KOPT] and
                        % [GOPT, HQ, PQ] (the quadratic model), so that XPT(:, KNEW_TR) becomes XNEW = XOPT + D.
                        xdrop(:) = xpt(:, knew_tr);
                        xosav(:) = xpt(:, kopt);
                        [idz, bmat, zmat] = powalg_obj.updateh(knew_tr, kopt, d, xpt, idz, bmat, zmat);
                        [kopt, fval, xpt] = update_lincoa_obj.updatexf(knew_tr, ximproved, f, xosav + d, kopt, fval, xpt);
                        [gopt, hq, pq] = update_lincoa_obj.updateq(idz, knew_tr, ximproved, bmat, d, moderr, xdrop, xosav, xpt, zmat, gopt, hq, pq);

                        % Establish the alternative model, namely the least Frobenius norm interpolant. Replace
                        % the current model with the alternative model if the recent few (three) alternative
                        % models are more accurate in predicting the function value of XOPT + D.
                        [qalt_better, gopt, pq, hq, galt, pqalt] = update_lincoa_obj.tryqalt(idz, bmat, fval - fval(kopt), xpt(:, kopt), xpt, zmat, qalt_better, gopt, pq, hq, galt, pqalt);
                        if ~(all(infnan_obj.is_finite(gopt), 'all') && all(infnan_obj.is_finite(hq), 'all') && all(infnan_obj.is_finite(pq), 'all'))
                            info = infos_obj.NAN_INF_MODEL;
                            break
                        end

                        % Update RESCON if XOPT is changed.
                        % Zaikun 20221115: Shouldn't we do it after DELTA is updated?
                        rescon = update_lincoa_obj.updateres(ximproved, amat, b, delta, linalg_obj.p_norm(d), xpt(:, kopt), rescon);
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
                % N.B.: The ACCURATE_MOD here plays a similar role as the variable with the same name in UOBYQA,
                % NEWUOA, and BOBYQA. However, the definition of ACCURATE_MOD here is different from that in
                % those solvers, which do not only check whether DNORM is small in recent iterations, but also
                % verify a curvature condition that really indicates that recent models are sufficiently
                % accurate. Here, however, we are not really sure whether they are accurate or not. Therefore,
                % ACCURATE_MOD is not the best name, but we keep it to align with the other solvers.
                accurate_mod = all(dnorm_rec <= rho, 'all') || all(dnorm_rec(2:numel(dnorm_rec)) <= 0.2 * rho, 'all');
                % Powell's version (note that size(dnorm_rec) = 5 in his implementation):
                %accurate_mod = all(dnorm_rec <= HALF * rho) .or. all(dnorm_rec(3:size(dnorm_rec)) <= TENTH * rho)
                % CLOSE_ITPSET: Are the interpolation points close to XOPT?
                distsq(:) = fortran.sum(fortran.power((xpt - xpt(:, kopt)), 2), 1);
                %%MATLAB: distsq = sum((xpt - xpt(:, kopt)).^2)  % Implicit expansion
                close_itpset = all(distsq <= 4.0 * fortran.power(delta, 2), 'all'); % Powell's NEWUOA code.
                % Below are some alternative definitions of CLOSE_ITPSET.
                % N.B.: The threshold for CLOSE_ITPSET is at least DELBAR, the trust region radius for GEOSTEP.
                % %close_itpset = all(distsq <= 4.0_RP * rho**2)  ! Powell's UOBYQA code.
                % %close_itpset = all(distsq <= max(delta**2, 4.0_RP * rho**2))  ! Powell's code.
                % %close_itpset = all(distsq <= max((TWO * delta)**2, (TEN * rho)**2))  ! Powell's BOBYQA code.
                % ADEQUATE_GEO: Is the geometry of the interpolation set "adequate"?
                adequate_geo = (shortd && accurate_mod) || close_itpset;
                % SMALL_TRRAD: Is the trust-region radius small? This indicator seems not impactive in practice.
                small_trrad = (max(delta, dnorm) <= rho); % Behaves the same as Powell's version.
                %small_trrad = (delsav <= rho)  ! Powell's code. DELSAV = unupdated DELTA.

                % IMPROVE_GEO and REDUCE_RHO are defined as follows.
                % N.B.: If SHORTD is TRUE at the very first iteration, then REDUCE_RHO will be set to TRUE.

                % BAD_TRSTEP (for IMPROVE_GEO): Is the last trust-region step bad?
                bad_trstep = (shortd || trfail || ratio <= eta1 || knew_tr == 0);
                improve_geo = bad_trstep && ~adequate_geo;
                % BAD_TRSTEP (for REDUCE_RHO): Is the last trust-region step bad?
                bad_trstep = (shortd || trfail || ratio <= 0 || knew_tr == 0);
                reduce_rho = bad_trstep && adequate_geo && small_trrad;

                % Equivalently, REDUCE_RHO can be set as follows. It shows that REDUCE_RHO is TRUE in two cases.
                % %bad_trstep = (shortd .or. trfail .or. ratio <= 0 .or. knew_tr == 0)
                % %reduce_rho = (shortd .and. accurate_mod) .or. (bad_trstep .and. close_itpset .and. small_trrad)

                % With REDUCE_RHO properly defined, we can also set IMPROVE_GEO as follows.
                % %bad_trstep = (shortd .or. trfail .or. ratio <= eta1 .or. knew_tr == 0)
                % %improve_geo = bad_trstep .and. (.not. reduce_rho) .and. (.not. close_itpset)

                % With IMPROVE_GEO properly defined, we can also set REDUCE_RHO as follows.
                % %bad_trstep = (shortd .or. trfail .or. ratio <= 0 .or. knew_tr == 0)
                % %reduce_rho = bad_trstep .and. (.not. improve_geo) .and. small_trrad

                % LINCOA never sets IMPROVE_GEO and REDUCE_RHO to TRUE simultaneously.
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

                if improve_geo
                    % XPT(:, KNEW_GEO) will become  XOPT + D below. KNEW_GEO /= KOPT unless there is a bug.
                    knew_geo = fix(fortran.maxloc(distsq, 'dim', 1));

                    % Set DELBAR, which will be used as the trust-region radius for the geometry-improving
                    % scheme GEOSTEP. Note that DELTA has been updated before arriving here.
                    delbar = max(consts_obj.TENTH * delta, rho); % Powell's code
                    %delbar = rho  ! Powell's UOBYQA code
                    %delbar = max(min(TENTH * sqrt(maxval(distsq)), HALF * delta), rho)  ! Powell's NEWUOA code
                    %delbar = max(min(TENTH * sqrt(maxval(distsq)), delta), rho)  ! Powell's BOBYQA code
                    % Find D so that the geometry of XPT will be improved when XPT(:, KNEW_GEO) becomes XOPT + D.
                    [feasible, d] = geometry_lincoa_obj.geostep(iact, idz, knew_geo, kopt, nact, amat, bmat, delbar, qfac, rescon, xpt, zmat, d);

                    % Calculate the next value of the objective function.
                    x(:) = xbase + (xpt(:, kopt) + d);
                    f = evaluate_obj.evaluatef(calfun, x);
                    nf = nf + 1;

                    % Evaluate the constraints. They are used only for printing messages.
                    constr_leq(:) = linalg_obj.matprod21(Aeq, x) - beq;
                    constr(:) = [reshape(xl(ixl) - x(ixl), [], 1); reshape(x(ixu) - xu(ixu), [], 1); reshape(-constr_leq, [], 1); reshape(constr_leq, [], 1); reshape(linalg_obj.matprod21(Aineq, x) - bineq, [], 1)];
                    cstrv = linalg_obj.maximum1([consts_obj.ZERO; reshape(constr, [], 1)]);

                    % Print a message about the function evaluation according to IPRINT.
                    message_obj.fmsg(solver, "Geometry", iprint, nf, delbar, f, x, 'cstrv', cstrv, 'constr', constr);
                    % Save X, F, CSTRV into the history.
                    [xhist, fhist, chist] = history_obj.savehist(nf, x, xhist, f, fhist, 'cstrv', cstrv, 'chist', chist);
                    % Save X, F, CSTRV into the filter.
                    [nfilt, cfilt, ffilt, xfilt] = selectx_obj.savefilt(cstrv, ctol, cweight, f, x, nfilt, cfilt, ffilt, xfilt);

                    % Check whether to exit.
                    subinfo = checkexit_obj.checkexit_con(maxfun, nf, cstrv, ctol, f, ftarget, x);
                    if subinfo ~= infos_obj.INFO_DFT
                        info = subinfo;
                        break
                    end

                    % QALT_BETTER is a boolean array indicating whether the recent few (three) alternative
                    % models are more accurate in predicting the function value at XOPT + D.
                    % Powell's code takes XOPT + D into account only if it is feasible.
                    % N.B.: Do NOT change the "<" in the comparison to "<="; otherwise, the result will not be
                    % reasonable if the two values being compared are both ZERO or INF.
                    moderr = f - fval(kopt) - powalg_obj.quadinc_d0(d, xpt, gopt, pq, 'hq', hq);
                    moderr_alt = f - fval(kopt) - powalg_obj.quadinc_d0(d, xpt, galt, pqalt);
                    qalt_better(:) = [reshape(qalt_better(2:numel(qalt_better)), [], 1); abs(moderr_alt) < consts_obj.TENTH * abs(moderr)];

                    % Is the newly generated X better than current best point?
                    ximproved = (f < fval(kopt) && feasible);

                    % Update [BMAT, ZMAT, IDZ] (represents H in the NEWUOA paper), [XPT, FVAL, KOPT] and
                    % [GOPT, HQ, PQ] (the quadratic model), so that XPT(:, KNEW_GEO) becomes XNEW = XOPT + D.
                    xdrop(:) = xpt(:, knew_geo);
                    xosav(:) = xpt(:, kopt);
                    [idz, bmat, zmat] = powalg_obj.updateh(knew_geo, kopt, d, xpt, idz, bmat, zmat);
                    [kopt, fval, xpt] = update_lincoa_obj.updatexf(knew_geo, ximproved, f, xosav + d, kopt, fval, xpt);
                    [gopt, hq, pq] = update_lincoa_obj.updateq(idz, knew_geo, ximproved, bmat, d, moderr, xdrop, xosav, xpt, zmat, gopt, hq, pq);

                    % Establish the alternative model, namely the least Frobenius norm interpolant. Replace the
                    % current model with the alternative model if the recent few (three) alternative models are
                    % more accurate in predicting the function value of XOPT + D.
                    % N.B.: Powell's code does this only if XOPT + D is feasible.
                    [qalt_better, gopt, pq, hq, galt, pqalt] = update_lincoa_obj.tryqalt(idz, bmat, fval - fval(kopt), xpt(:, kopt), xpt, zmat, qalt_better, gopt, pq, hq, galt, pqalt);
                    if ~(all(infnan_obj.is_finite(gopt), 'all') && all(infnan_obj.is_finite(hq), 'all') && all(infnan_obj.is_finite(pq), 'all'))
                        info = infos_obj.NAN_INF_MODEL;
                        break
                    end

                    % Update RESCON. Zaikun 20221115: Currently, UPDATERES does not update RESCON if XIMPROVED
                    % is FALSE. Shouldn't we do it whenever DELTA is updated? Have we MISUNDERSTOOD RESCON?
                    rescon = update_lincoa_obj.updateres(ximproved, amat, b, delta, linalg_obj.p_norm(d), xpt(:, kopt), rescon);
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
                    message_obj.rhomsg(solver, iprint, nf, delta, fval(kopt), rho, xbase + xpt(:, kopt), 'cstrv', cstrv, 'constr', constr);
                    % DNORM_REC is corresponding to the latest function evaluations with the current RHO.
                    % Update it after reducing RHO.
                    dnorm_rec(:) = consts_obj.REALMAX;
                end % End of IF (REDUCE_RHO). The procedure of reducing RHO ends.

                % Shift XBASE if XOPT may be too far from XBASE.
                % Powell's original criterion for shifting XBASE: before a trust region step or a geometry step,
                % shift XBASE if SUM(XOPT**2) >= 1.0E3*DELTA**2.
                if fortran.sum(fortran.power(xpt(:, kopt), 2), 'all') >= 1000.0 * fortran.power(delta, 2)
                    % Other possible criteria: SUM(XOPT**2) >= 1.0E4*DELTA**2, SUM(XOPT**2) >= 1.0E3*RHO**2.
                    b(:) = b - linalg_obj.matprod12(xpt(:, kopt), amat);
                    [xbase, xpt, bmat, hq] = shiftbase_obj.shiftbase_lfqint(kopt, xbase, xpt, zmat, bmat, pq, hq, 'idz', idz);
                    % SHIFTBASE shifts XBASE to XBASE + XOPT and XOPT to 0.
                    pqalt(:) = powalg_obj.omega_mul(idz, zmat, fval - fval(kopt));
                    galt(:) = linalg_obj.matprod21(bmat(:, 1:npt), fval - fval(kopt)) + powalg_obj.hess_mul(xpt(:, kopt), xpt, pqalt);
                end

                % Report the current best value, and check if user asks for early termination.
                if ~ismember('callback_fcn', ipObj.UsingDefaults)
                    % FIXME: CVAL(KOP) is WRONG! CVAL is not updated.
                    terminate = callback_fcn(xbase + xpt(:, kopt), fval(kopt), nf, tr, 'cstrv', cval(kopt));
                    if terminate
                        info = infos_obj.CALLBACK_TERMINATE;
                        break
                    end
                end

            end % End of DO TR = 1, MAXTR. The iterative procedure ends.

            % Return from the calculation, after trying the Newton-Raphson step if it has not been tried yet.
            if info == infos_obj.SMALL_TR_RADIUS && shortd && dnorm > consts_obj.TENTH * rhoend && nf < maxfun
                x(:) = xbase + (xpt(:, kopt) + d);
                f = evaluate_obj.evaluatef(calfun, x);
                nf = nf + 1;
                constr_leq(:) = linalg_obj.matprod21(Aeq, x) - beq;
                constr(:) = [reshape(xl(ixl) - x(ixl), [], 1); reshape(x(ixu) - xu(ixu), [], 1); reshape(-constr_leq, [], 1); reshape(constr_leq, [], 1); reshape(linalg_obj.matprod21(Aineq, x) - bineq, [], 1)];
                cstrv = linalg_obj.maximum1([consts_obj.ZERO; reshape(constr, [], 1)]);
                % Print a message about the function evaluation according to IPRINT.
                % Zaikun 20230512: DELTA has been updated. RHO is only indicative here. TO BE IMPROVED.
                message_obj.fmsg(solver, "Trust region", iprint, nf, rho, f, x, 'cstrv', cstrv, 'constr', constr);
                % Save X, F, CSTRV into the history.
                [xhist, fhist, chist] = history_obj.savehist(nf, x, xhist, f, fhist, 'cstrv', cstrv, 'chist', chist);
                % Save X, F, CSTRV into the filter.
                [nfilt, cfilt, ffilt, xfilt] = selectx_obj.savefilt(cstrv, ctol, cweight, f, x, nfilt, cfilt, ffilt, xfilt);
            end

            % Return the best calculated values of the variables.
            kopt = selectx_obj.selectx(ffilt(1:nfilt), cfilt(1:nfilt), cweight, ctol);
            x(:) = xfilt(:, kopt);
            f = ffilt(kopt);
            constr_leq(:) = linalg_obj.matprod21(Aeq, x) - beq;
            constr(:) = [reshape(xl(ixl) - x(ixl), [], 1); reshape(x(ixu) - xu(ixu), [], 1); reshape(-constr_leq, [], 1); reshape(constr_leq, [], 1); reshape(linalg_obj.matprod21(Aineq, x) - bineq, [], 1)];
            cstrv = linalg_obj.maximum1([consts_obj.ZERO; reshape(constr, [], 1)]);

            % Deallocate IXL and IXU as they have finished their job.


            % Arrange CHIST, FHIST, and XHIST so that they are in the chronological order.
            [xhist, fhist, chist] = history_obj.rangehist(nf, xhist, fhist, 'chist', chist);

            % Print a return message according to IPRINT.
            message_obj.retmsg(solver, info, iprint, nf, f, x, 'cstrv', cstrv, 'constr', constr);

            %====================%
            %  Calculation ends  %
            %====================%

            % Postconditions
            if consts_obj.DEBUGGING
                debug_obj.assert(nf <= maxfun, "NF <= MAXFUN", srname);
                debug_obj.assert(numel(x) == n && ~any(infnan_obj.is_nan_sp(x), 'all'), "SIZE(X) == N, X does not contain NaN", srname);
                debug_obj.assert(~(infnan_obj.is_nan_sp(f) || infnan_obj.is_posinf(f)), "F is not NaN/+Inf", srname);
                debug_obj.assert(size(xhist, 1) == n && size(xhist, 2) == maxxhist, "SIZE(XHIST) == [N, MAXXHIST]", srname);
                debug_obj.assert(~any(infnan_obj.is_nan_sp(xhist(:, 1:min(nf, maxxhist))), 'all'), "XHIST does not contain NaN", srname);
                % The last calculated X can be Inf (finite + finite can be Inf numerically).
                debug_obj.assert(numel(fhist) == maxfhist, "SIZE(FHIST) == MAXFHIST", srname);
                debug_obj.assert(~any(infnan_obj.is_nan_sp(fhist(1:min(nf, maxfhist))) | infnan_obj.is_posinf(fhist(1:min(nf, maxfhist))), 'all'), "FHIST does not contain NaN/+Inf", srname);
                debug_obj.assert(numel(chist) == maxchist, "SIZE(CHIST) == MAXCHIST", srname);
                debug_obj.assert(~any(chist(1:min(nf, maxchist)) < 0 | infnan_obj.is_nan_sp(chist(1:min(nf, maxchist))) | infnan_obj.is_posinf(chist(1:min(nf, maxchist))), 'all'), "CHIST does not contain negative values or NaN/+Inf", srname);
                nhist = min([nf, maxfhist, maxchist], [], 'all');
                debug_obj.assert(~any(selectx_obj.isbetter10(fhist(1:nhist), chist(1:nhist), f, cstrv, ctol), 'all'), "No point in the history is better than X", srname);
            end

        end

    end
end