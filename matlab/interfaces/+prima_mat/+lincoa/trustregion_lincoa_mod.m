classdef trustregion_lincoa_mod
    %--------------------------------------------------------------------------------------------------%
    % This module provides subroutines concerning the trust-region calculations of LINCOA.
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
    % Last Modified: Saturday, March 09, 2024 PM12:09:39
    %--------------------------------------------------------------------------------------------------%

    methods
        function [iact, nact, qfac, rfac, s, ngetact] = trstep(~, amat, delta, gopt_in, hq_in, pq_in, rescon, tol, xpt, iact, nact, qfac, rfac, s, varargin)
            %--------------------------------------------------------------------------------------------------%
            % This subroutine solves
            %       minimize Q(XOPT + D)  s.t. ||D|| <= DELTA, AMAT^T*D <= B.
            % It is assumed that D = 0 is feasible, namely B >= 0 except for rounding errors. See Powell 2015
            % for details.
            %
            % AMAT, B, XPT, GOPT, HQ, PQ, NACT, IACT, RESCON, QFAC and RFAC are the same as the terms with these
            % names in LINCOB.
            %
            % S is the total calculated step so far from the trust region centre, its final value being given by
            %   the sequence of CG iterations, which terminate if the trust region boundary is reached.
            % G is always the gradient of the model at the current S.
            % D is the search direction of each line search.
            % RESCON: If RESCON(J) is negative, then |RESCON(J)| must be no less than the trust region radius,
            %   so that the J-th constraint can be ignored.
            % RESNEW: A negative value of RESNEW(J) indicates that the J-th constraint does not restrict the CG
            %   steps of the current trust region calculation, a zero value of RESNEW(J) indicates that the J-th
            %   constraint is active, and otherwise RESNEW(J) is set to the greater of TINYCV and the actual
            %   residual of the J-th constraint for the current S.
            % RESACT holds the residuals of the active constraints, which may be positive.
            %--------------------------------------------------------------------------------------------------%

            % Common modules
            consts_obj = prima_mat.common.consts_mod();
            debug_obj = prima_mat.common.debug_mod();
            infnan_obj = prima_mat.common.infnan_mod();
            linalg_obj = prima_mat.common.linalg_mod();
            powalg_obj = prima_mat.common.powalg_mod();

            % Solver-specific modules
            getact_obj = prima_mat.lincoa.getact_mod();


            % Inputs
            % AMAT(N, M)

            % GOPT_IN(N)
            % HQ_IN(N, N)
            % PQ_IN(NPT)
            % RESCON(M)

            % XPT(N, NPT)

            % In-outputs
            % IACT(M); Will be updated in GETACT
            % Will be updated in GETACT
            % QFAC(N, N); Will be updated in GETACT
            % RFAC(N, N); Will be updated in GETACT

            % Outputs
            % S(N)


            % Local variables
            srname = "TRSTEP";


            jsav = NaN;


            ad = NaN(size(amat, 2), 1);
            alpha = NaN;
            alphm = NaN;
            alpht = NaN;
            beta = NaN;
            d = NaN(numel(gopt_in), 1);


            dg = NaN;
            dhd = NaN;
            dproj = NaN(numel(gopt_in), 1);
            ds = NaN;
            frac = NaN(size(amat, 2), 1);
            g = NaN(numel(gopt_in), 1);
            gamma = NaN;
            gopt = NaN(numel(gopt_in), 1);
            hd = NaN(numel(gopt_in), 1);
            hq = NaN(size(hq_in, 1), size(hq_in, 2));

            orthtol = NaN;
            pg = NaN(numel(gopt_in), 1);
            pq = NaN(numel(pq_in), 1);
            psd = NaN(numel(gopt_in), 1);

            resact = NaN(size(amat, 2), 1);
            resid = NaN;
            resnew = NaN(size(amat, 2), 1);
            restmp = NaN(size(amat, 2), 1);
            sold = NaN(numel(s), 1);
            sqrtd = NaN;


            % Sizes.
            m = size(amat, 2);
            n = fix(numel(gopt_in));
            npt = fix(numel(pq_in));

            % Preconditions
            if consts_obj.DEBUGGING
                debug_obj.assert(m >= 0, "M >= 0", srname);
                debug_obj.assert(n >= 1 && npt >= n + 2, "N >= 1, NPT >= N + 2", srname);
                debug_obj.assert(delta > 0, "DELTA > 0", srname);
                debug_obj.assert(size(amat, 1) == n && size(amat, 2) == m, "SIZE(AMAT) == [N, M]", srname);
                debug_obj.assert(size(hq_in, 1) == n && linalg_obj.issymmetric(hq_in), "HQ is n-by-n and symmetric", srname);
                debug_obj.assert(numel(rescon) == m, "SIZE(RESCON) == M", srname);
                debug_obj.assert(size(xpt, 1) == n && size(xpt, 2) == npt, "SIZE(XPT) == [N, NPT]", srname);
                debug_obj.assert(all(infnan_obj.is_finite(xpt), 'all'), "XPT is finite", srname);
                debug_obj.assert(nact >= 0 && nact <= min(m, n), "0 <= NACT <= MIN(M, N)", srname);
                debug_obj.assert(numel(iact) == m, "SIZE(IACT) == M", srname);
                debug_obj.assert(all(iact(1:nact) >= 1 & iact(1:nact) <= m, 'all'), "1 <= IACT <= M", srname);
                debug_obj.assert(size(qfac, 1) == n && size(qfac, 2) == n, "SIZE(QFAC) == [N, N]", srname);
                orthtol = max(fortran.power(consts_obj.TEN, max(-10, -consts_obj.MAXPOW10)), min(0.1, fortran.power(consts_obj.TEN, min(8, consts_obj.MAXPOW10)) * consts_obj.EPS * double(n)));
                debug_obj.assert(linalg_obj.isorth(qfac, 'tol', orthtol), "QFAC is orthogonal", srname);
                debug_obj.assert(size(rfac, 1) == n && size(rfac, 2) == n, "SIZE(RFAC) == [N, N]", srname);
                debug_obj.assert(linalg_obj.istriu(rfac), "RFAC is upper triangular", srname);
                debug_obj.assert(size(qfac, 1) == n && size(qfac, 2) == n, "SIZE(QFAC) == [N, N]", srname);
            end

            %====================%
            % Calculation starts %
            %====================%

            % Scale the problem if GOPT contains large values. Otherwise, floating point exceptions may occur.
            % Note that the trust-region step is scale invariant.
            % N.B.: It is faster and safer to scale by multiplying a reciprocal than by division. See
            % https://fortran-lang.discourse.group/t/ifort-ifort-2021-8-0-1-0e-37-1-0e-38-0/
            if max(abs(gopt_in), [], 'all') > 1.0e12
                % The threshold is empirical.
                modscal = max(consts_obj.TWO * consts_obj.REALMIN, consts_obj.ONE / max(abs(gopt_in), [], 'all')); % MAX: precaution against underflow.
                gopt(:) = gopt_in * modscal;
                pq(:) = pq_in * modscal;
                hq(:, :) = hq_in * modscal;
            else
                gopt(:) = gopt_in;
                pq(:) = pq_in;
                hq(:, :) = hq_in;
            end

            % Return if G is not finite. Otherwise, GETACT will fail in the debugging mode.
            ipObj = inputParser();
            addParameter(ipObj, 'ngetact', NaN);
            parse(ipObj, varargin{:});
            ngetact = ipObj.Results.ngetact;
            if ~infnan_obj.is_finite(sum(abs(gopt), 'all'))
                s(:) = consts_obj.ZERO;
                if nargout >= 6
                    ngetact = 0;
                end
                return
            end

            % Set the initial elements of RESNEW, RESACT and S.

            % 1. RESNEW(J) < 0 indicates that the J-th constraint does not restrict the CG steps of the current
            % trust region calculation. In other words, RESCON >= DELTA.
            % 2. RESNEW(J) = 0 indicates that J is an entry of IACT(1:NACT).
            % 3. RESNEW(J) > 0 means that RESNEW(J) = max(B(J) - AMAT(:, J)^T*(XOPT+S), TINYCV), where S is the
            % step up to now, calculated by a sequence of (truncated) CG iterations.
            % N.B.: The order of the following lines is important, as the later ones override the earlier.
            resnew(:) = rescon;
            resnew(linalg_obj.trueloc(rescon >= 0)) = max(consts_obj.TINYCV, rescon(linalg_obj.trueloc(rescon >= 0)));
            resnew(linalg_obj.trueloc(rescon >= delta)) = -consts_obj.ONE;
            %%MATLAB:
            %%resnew = rescon; resnew(rescon >= 0) = max(TINYCV, rescon(rescon >= 0)); resnew(rescon >= delta) = -1;
            resnew(iact(1:nact)) = consts_obj.ZERO;

            % RESACT contains the constraint residuals of the constraints in IACT(1:NACT), namely the values
            % of B(J) - AMAT(:, J)^T*(XOPT+S) for the J in IACT(1:NACT). Here, IACT(1:NACT) is a set of
            % indicates such that the columns of AMAT(:, IACT(1:NACT)) form a basis of the constraint gradients
            % in the "active set". For the definition of the "active set", see (3.5) of Powell (2015) and the
            % comments at the beginning of the GETACT subroutine.
            % N.B.: Between two calls of GETACT, S is updated in the orthogonal complement of the "active"
            % gradients (i.e., null space of the "active" constraints). Therefore, RESACT remains unchanged.
            % RESACT is changed right after GETACT is called if the first search direction D is not PSD
            % but PSD + GAMMA * DPROJ.
            resact(1:nact) = rescon(iact(1:nact));

            g(:) = gopt;
            delsq = delta * delta;
            s(:) = consts_obj.ZERO;
            ss = consts_obj.ZERO;
            reduct = consts_obj.ZERO;
            ngetact_loc = 0;
            newact = true;

            % ITERCG is the number of CG iterations corresponding to the current "active set" obtained by
            % calling GETACT. These CG iterations are restricted in the orthogonal complement of the active
            % gradients (i.e., null space of the active constraints).
            % The following initial value of ITERCG is an artificial value that is not used. It is to entertain
            % Fortran compilers (can it be be removed?).
            itercg = -1;

            % What is the THEORETICAL upper bound of ITER? For the moment, we set the following MAXITER.
            % The formulation of MAXITER below contains a precaution against overflow. In MATLAB/Python/Julia/R,
            % we can write maxiter = min(10000, 10*(m + n))
            maxiter = fix(min(fortran.power(10, min(4, floor(log10(double(intmax('int64')))))), 10 * fix(m + n)));
            for iter = 1:maxiter                % Powell's code is essentially a DO WHILE loop. We impose an explicit MAXITER.
                if newact
                    % GETACT picks the active set for the current S. It also sets PSD to the vector closest to
                    % -G that is orthogonal to the normals of the active constraints. PSD is scaled to have
                    % length 0.2*DELTA. Then a move of PSD from S is allowed by the linear constraints: PSD
                    % reduces the values of the nearly active constraints; it changes the inactive constraints
                    % by at most 0.2*DELTA, but the residuals of these constraints at no less than 0.2*DELTA.
                    % N.B.: The magic number 0.2 appears also in GETACT (TDEL = 0.2_RP * DELTA). It works well.
                    ngetact_loc = ngetact_loc + 1;
                    [iact, nact, qfac, resact, resnew, rfac, psd] = getact_obj.getact(amat, delta, g, iact, nact, qfac, resact, resnew, rfac, psd);
                    dd = linalg_obj.inprod(psd, psd);
                    if dd <= consts_obj.EPS * delsq || infnan_obj.is_nan_sp(dd)
                        % Powell's code: IF (DD <= 0) THEN
                        break
                    end
                    psd(:) = (0.2 * delta / fortran.sqrt(dd)) * psd;

                    % If the modulus of the residual of an "active constraint" is substantial (i.e., more than
                    % 1.0E-4*DELTA), then modify the searching direction PSD by a projection step to the
                    % boundaries of the "active constraint". This modified step will reduce the constraint
                    % residuals of the "active constraints" (see the update of RESACT below). The motivation is
                    % that the constraints in the "active set" are presumed to be active, and hence should have
                    % zero residuals (no constraint is violated, as the current method is feasible). According
                    % to a test on 20220821, this modification is important for the performance of LINCOA.
                    % N.B.:
                    % 1. The residual of the constraint A*X <= B is defined as B - A*X. It is not the constraint
                    % violation. Indeed, the constraint violations of the iterates are 0 in the current method.
                    % 2. We prefer `ANY(X > Y)` to `MAXVAL(X) > Y`, as Fortran standards do not specify
                    % MAXVAL(X) when X contains NaN, and MATLAB/Python/R/Julia behave differently in this
                    % respect. Moreover, MATLAB defines max(X) = [] if X == [], differing from mathematics
                    % and other languages.
                    gamma = consts_obj.ZERO; % The steplength of the projection step to be taken.
                    if any(resact(1:nact) > 1.0e-4 * delta, 'all')
                        % Set DPROJ to the shortest move (projection step) from S to the boundaries of the
                        % active constraints. We will use DPROJ to modify PSD.
                        dproj(:) = linalg_obj.matprod21(qfac(:, 1:nact), linalg_obj.solve(rfac(1:nact, 1:nact)', resact(1:nact)));
                        %%MATLAB: dproj = qfac(:, 1:nact) * (rfac(1:nact, 1:nact)' \ resact(1:nact))

                        % The vector DPROJ is also the shortest move from S + PSD to the boundaries of the
                        % active constraints (this is because PSD is parallel to the boundaries of the active
                        % constraints). Set GAMMA to the greatest steplength of this move that satisfies both
                        % the trust region bound and the linear constraints.
                        ds = linalg_obj.inprod(dproj, s + psd);
                        dd = sum(fortran.power(dproj, 2), 'all');
                        resid = delsq - sum(fortran.power((s + psd), 2), 'all');
                        % Powell's condition for the following IF: RESID > 0.
                        if resid > 0 && dd > consts_obj.EPS * delsq && ~infnan_obj.is_nan_sp(ds)
                            % Set GAMMA to the greatest value so that S + PSD + GAMMA*DPROJ satisfies the trust
                            % region bound. SQRTD: square root of a discriminant. Powell's code for SQRTD is
                            % SQRT(DS * DS + DD * RESID), which may be below ABS(DS) due to underflow in DS*DS.
                            sqrtd = max([fortran.sqrt(ds * ds + dd * resid), abs(ds), fortran.sqrt(dd * resid)], [], 'all');
                            if ds <= 0
                                gamma = (sqrtd - ds) / dd;
                            else
                                gamma = resid / (sqrtd + ds);
                            end
                            % GAMMA < 0 should not happen. GAMMA can be 0 or NaN when, e.g., DS or DD becomes
                            % Inf. Powell's code does not handle this.
                            if gamma < 0 || ~infnan_obj.is_finite(gamma)
                                gamma = 0;
                            end

                            % Reduce GAMMA so that the move along DPROJ also satisfies the linear constraints.
                            ad(:) = -consts_obj.ONE;
                            ad(linalg_obj.trueloc(resnew > 0)) = linalg_obj.matprod12(dproj, amat(:, linalg_obj.trueloc(resnew > 0)));
                            frac(:) = consts_obj.ONE;
                            restmp(linalg_obj.trueloc(ad > 0)) = resnew(linalg_obj.trueloc(ad > 0)) - linalg_obj.matprod12(psd, amat(:, linalg_obj.trueloc(ad > 0)));
                            frac(linalg_obj.trueloc(ad > 0)) = restmp(linalg_obj.trueloc(ad > 0)) ./ ad(linalg_obj.trueloc(ad > 0));
                            gamma = min([gamma; consts_obj.ONE; reshape(frac, [], 1)], [], 'all'); % GAMMA = MINVAL([GAMMA, ONE, FRAC(TRUELOC(AD>0))])

                        end
                    end

                    % Set the next direction for seeking a reduction in the model function subject to the trust
                    % region bound and the linear constraints.
                    % Do NOT write D = PSD + GAMMA*DPROJ, as DPROJ may contain NaN/Inf, in which case GAMMA = 0.
                    if gamma > 0
                        d(:) = psd + gamma * dproj; % Modified searching direction.
                        itercg = -1;
                    else
                        d(:) = psd; % Original searching direction.
                        itercg = 0;
                    end
                end
                itercg = itercg + 1;
                % After the above line, ITERCG = 0 iff GETACT has been just called, and D is not PSD but a
                % modified step.

                % Set ALPHA to the steplength from S along D to the trust region boundary. Return if the first
                % derivative term of this step is sufficiently small or if no further progress is possible.
                resid = delsq - ss;
                dg = linalg_obj.inprod(d, g);
                ds = linalg_obj.inprod(d, s);
                dd = linalg_obj.inprod(d, d);
                % Powell's condition for the following IF: (RESID <= 0 .OR. DG >= 0). If DD is tiny (so is DS),
                % ALPHA may be mistakenly calculated as a huge value due to rounding errors, as observed on
                % 20221205. Therefore, we exit when DD is small. The test for DG is covered by the IF after the
                % calculation of ALPHA.
                if resid <= 0 || dd <= consts_obj.EPS * delsq || infnan_obj.is_nan_sp(ds)
                    break
                end
                % SQRTD: square root of a discriminant. Powell's code for SQRTD is SQRT(DS * DS + DD * RESID),
                % which may be below ABS(DS) due to underflow in DS*DS.
                sqrtd = max([fortran.sqrt(ds * ds + dd * resid), abs(ds), fortran.sqrt(dd * resid)], [], 'all');
                if ds <= 0
                    alpha = (sqrtd - ds) / dd;
                else
                    alpha = resid / (sqrtd + ds);
                end
                % ALPHA < 0 should not happen. ALPHA can be 0 or NaN when, e.g., DS or DD becomes Inf. Powell's
                % code does not handle this.
                if alpha <= 0 || ~infnan_obj.is_finite(alpha)
                    break
                end

                % Powell's condition for the following IF: -ALPHA * DG <= TOL * REDUCT. Note that the EXIT
                % will be triggered if DG >= 0, as ALPHA >= 0.
                if -alpha * dg <= tol * reduct || infnan_obj.is_nan_sp(alpha * dg)
                    break
                end

                % Set DHD to the curvature of the model along D. Then reduce ALPHA if necessary to the value
                % that minimizes the model.
                hd(:) = powalg_obj.hess_mul(d, xpt, pq, 'hq', hq);
                dhd = linalg_obj.inprod(d, hd);
                alpht = alpha;
                if dg + alpha * dhd > 0
                    alpha = -dg / dhd;
                end

                % Make a further reduction in ALPHA if necessary to preserve feasibility.
                alphm = alpha;
                ad(:) = -consts_obj.ONE;
                ad(linalg_obj.trueloc(resnew > 0)) = linalg_obj.matprod12(d, amat(:, linalg_obj.trueloc(resnew > 0)));
                frac(:) = alpha;
                frac(linalg_obj.trueloc(ad > 0)) = resnew(linalg_obj.trueloc(ad > 0)) ./ ad(linalg_obj.trueloc(ad > 0));
                frac(linalg_obj.trueloc(infnan_obj.is_nan(frac))) = alpha;
                jsav = 0;
                if any(frac < alpha, 'all')
                    jsav = fix(fortran.minloc(frac, 'dim', 1));
                    alpha = frac(jsav);
                end
                %----------------------------------------------------------------------------------------------%
                % Alternatively, JSAV and ALPHA can be calculated as below.
                % %JSAV = INT(MINLOC([ALPHA, FRAC], DIM=1), KIND(JSAV)) - 1_IK
                % %ALPHA = MINVAL([ALPHA, FRAC])  ! This line cannot be exchanged with the last.
                % We prefer our implementation as the code is more explicit; in addition, it is more flexible:
                % we can change the condition ANY(FRAC < ALPHA) to ANY(FRAC < (1 - EPS) * ALPHA) or
                % ANY(FRAC < (1 + EPS) * ALPHA), depending on whether we believe a false positive or a false
                % negative of JSAV > 0 is more harmful.
                %----------------------------------------------------------------------------------------------%

                % Post-process ALPHA according to some prior information.
                % N.B.:
                % 1. Since we set ALPHA=1 when ITERCG=0, the ALPHA calculated above is needed only if ITERCG>0.
                % 2. Zaikun 20220821: In theory, shouldn't this post-processing change nothing? According to
                % a test on 20220821, it does change ALPHA sometimes. Strange! Why?
                if itercg == 0
                    % Iff GETACT has been called, and D is not PSD but a modified step.
                    % By the definition of D, ALPHA = ONE is the largest ALPHA so that S + ALPHA*D satisfies the
                    % linear and trust region constraints.
                    alpha = consts_obj.ONE;
                elseif itercg == 1 && gamma <= 0
                    % Iff GETACT has been called, and D is not modified.
                    % Due to the scaling of PSD, S + D satisfies the linear and trust region constraints.
                    alpha = max(alpha, consts_obj.ONE);
                else
                    alpha = max(alpha, consts_obj.ZERO);
                end

                % Set ALPHA to the minimum between ALPHA and ALPHM, namely the steplength obtained by minimizing
                % the quadratic model along D.
                alpha = min(alpha, alphm);

                % Update S, G.
                sold(:) = s;
                s(:) = s + alpha * d;
                ss = sum(fortran.power(s, 2), 'all');
                if ~infnan_obj.is_finite(ss)
                    s(:) = sold;
                    break
                end
                g(:) = g + alpha * hd;
                if ~infnan_obj.is_finite(sum(abs(g), 'all'))
                    break
                end

                % Update RESNEW.
                restmp(:) = resnew - alpha * ad; % Only RESTMP(TRUELOC(RESNEW > 0)) is needed.
                resnew(linalg_obj.trueloc(resnew > 0)) = max(consts_obj.TINYCV, restmp(linalg_obj.trueloc(resnew > 0)));
                %%MATLAB: mask = (resnew > 0); resnew(mask) = max(TINYCV, resnew(mask) - alpha * ad(mask));

                % Update RESACT. This is done iff GETACT has been called, and D is not PSD but a modified step.
                %----------------------------------------------------------------------------------------------%
                % Zaikun 20220821: There seems be a typo here. Powell's original code does not take ALPHA into
                % account. Then RESACT seems to correspond to S + D, where D is defined as PSD + GAMMA*DPROJ
                % during the modification procedure after GETACT is called. Without this modification, RESACT
                % would remain unchanged because D = PSD, which is in the null space of the active constraints.
                % The GAMMA*DPROJ component in the modified step D reduces RESACT by GAMMA*RESACT. However,
                % since S is updated to S + ALPHA*D, shouldn't RESACT be reduced by ALPHA*GAMMA*RESACT?
                % Note that Powell chose to update RESACT after ALPHA is calculated (instead of right after
                % GAMMA is calculated), which might be an indication that he wanted to take ALPHA into account.
                % In the following code, we try correcting this apparent typo, but it has little impact on the
                % performance of LINCOA according to a test on 20220821.
                if itercg == 0
                    resact(1:nact) = (consts_obj.ONE - alpha * gamma) * resact(1:nact);
                    %resact(1:nact) = (ONE - gamma) * resact(1:nact)  ! Powell's code.

                end
                %----------------------------------------------------------------------------------------------%

                % Update REDUCT, the reduction up to now.
                reduct = reduct - alpha * (dg + consts_obj.HALF * alpha * dhd);
                if reduct <= 0 || infnan_obj.is_nan_sp(reduct)
                    s(:) = sold;
                    break
                end

                % Test for termination.
                if alpha >= alpht || -alphm * (dg + consts_obj.HALF * alphm * dhd) <= tol * reduct
                    break
                end

                % Branch to a new loop if there is a new active constraint.
                % When JSAV > 0, Powell's code branches back with NEWACT = .TRUE. only if ||S|| <= 0.8*DELTA,
                % and it exits if ||S|| > 0.8*DELTA, as mentioned at the end of Section 3 of Powell 2015. The
                % motivation seems to avoid small steps that changes the active set, because GETACT is expensive
                % in flops. However, according to a test on 20220820, removing this condition (essentially
                % replacing it with ||S|| < DELTA) improves the performance of LINCOA a bit. This may lead to
                % small steps, but tiny steps will lead to tiny reductions and trigger an exit.
                newact = (jsav > 0);
                if newact
                    continue
                end

                % If N-NACT CG iterations has been taken in the current null space (corresponding to the
                % current "active set"), then, in theory, a stationary point in this subspace has been found.
                % If the "active set" is the true active set, then a stationary point of the
                % linearly-constrained trust region subproblem is found. So a termination is reasonable.
                % However, the "active set" is not precisely the true active set, is it? See (3.5) of Powell
                % (2015) and the comments at the beginning of the GETACT subroutine. Also, we should take into
                % account the modification after GETACT is called.
                if itercg >= n - nact
                    % ITERCG > N - NACT is impossible.
                    break
                end

                % Calculate the next search direction, which is conjugate to the previous one if ITERCG /= NACT.
                % N.B.: NACT < 0 is impossible unless GETACT is buggy; NACT = 0 can happen, particularly if
                % there is no constraint. In theory, the code for the second case below covers the first as well.
                if nact <= 0
                    pg(:) = g;
                else
                    pg(:) = linalg_obj.matprod21(qfac(:, nact + 1:n), linalg_obj.matprod12(g, qfac(:, nact + 1:n)));
                    %%MATLAB: pg = qfac(:, nact+1:n) * (g' * qfac(:, nact+1:n))';
                end

                if itercg == 0
                    % Iff GETACT has been called, and D is not PSD but a modified step.
                    beta = consts_obj.ZERO;
                else
                    beta = linalg_obj.inprod(pg, hd) / dhd;
                end
                d(:) = -pg + beta * d;
            end

            if nargout >= 6
                ngetact = ngetact_loc;
            end

            %====================%
            %  Calculation ends  %
            %====================%

            % Postconditions
            if consts_obj.DEBUGGING
                debug_obj.assert(numel(s) == n && all(infnan_obj.is_finite(s), 'all'), "SIZE(S) == N, S is finite", srname);
                % Due to rounding, it may happen that ||S|| > DELTA, but ||S|| > 2*DELTA is highly improbable.
                debug_obj.assert(linalg_obj.p_norm(s) <= consts_obj.TWO * delta, "||S|| <= 2*DELTA", srname);
                debug_obj.assert(nact >= 0 && nact <= min(m, n), "0 <= NACT <= MIN(M, N)", srname);
                debug_obj.assert(size(qfac, 1) == n && size(qfac, 2) == n, "SIZE(QFAC) == [N, N]", srname);
                debug_obj.assert(linalg_obj.isorth(qfac, 'tol', orthtol), "QFAC is orthogonal", srname);
                debug_obj.assert(size(rfac, 1) == n && size(rfac, 2) == n, "SIZE(RFAC) == [N, N]", srname);
                debug_obj.assert(linalg_obj.istriu(rfac), "RFAC is upper triangular", srname);
                if nargout >= 6
                    debug_obj.assert(ngetact >= 1, "NGETACT >= 1", srname);
                end
            end

        end
        %--------------------------------------------------------------------------------------------------%
        % Zaikun 20220417:
        % For PG, the schemes below work evidently worse than the one above in a test on 20220417. Why?
        %-----------------------------------------------------------------------%
        % VERSION 1:
        % %pg = g - matprod(qfac(:, 1:nact), matprod(g, qfac(:, 1:nact)))
        %-----------------------------------------------------------------------%
        % VERSION 2:
        % %if (2 * nact < n) then
        % %    pg = g - matprod(qfac(:, 1:nact), matprod(g, qfac(:, 1:nact)))
        % %else
        % %    pg = matprod(qfac(:, nact + 1:n), matprod(g, qfac(:, nact + 1:n)))
        % %end if
        %-----------------------------------------------------------------------%
        %--------------------------------------------------------------------------------------------------%
        function delta = trrad(~, delta_in, dnorm, eta1, eta2, gamma1, gamma2, ratio)
            %--------------------------------------------------------------------------------------------------%
            % This function updates the trust region radius according to RATIO and DNORM.
            %--------------------------------------------------------------------------------------------------%

            % Generic module
            consts_obj = prima_mat.common.consts_mod();
            infnan_obj = prima_mat.common.infnan_mod();
            debug_obj = prima_mat.common.debug_mod();


            % Input
            % Current trust-region radius
            % Norm of current trust-region step
            % Ratio threshold for contraction
            % Ratio threshold for expansion
            % Contraction factor
            % Expansion factor
            % Reduction ratio

            % Outputs
            delta = NaN;

            % Local variables
            srname = "TRRAD";

            % Preconditions
            if consts_obj.DEBUGGING
                debug_obj.assert(delta_in >= dnorm && dnorm > 0, "DELTA_IN >= DNORM > 0", srname);
                debug_obj.assert(eta1 >= 0 && eta1 <= eta2 && eta2 < 1, "0 <= ETA1 <= ETA2 < 1", srname);
                debug_obj.assert(eta1 >= 0 && eta1 <= eta2 && eta2 < 1, "0 <= ETA1 <= ETA2 < 1", srname);
                debug_obj.assert(gamma1 > 0 && gamma1 < 1 && gamma2 > 1, "0 < GAMMA1 < 1 < GAMMA2", srname);
                % By the definition of RATIO in ratio.f90, RATIO cannot be NaN unless the actual reduction is
                % NaN, which should NOT happen due to the moderated extreme barrier.
                debug_obj.assert(~infnan_obj.is_nan_sp(ratio), "RATIO is not NaN", srname);
            end

            %====================%
            % Calculation starts %
            %====================%

            if ratio <= eta1
                delta = gamma1 * dnorm; % Powell's UOBYQA/NEWUOA.
                %delta = gamma1 * delta_in  ! Powell's COBYLA/LINCOA.
                %delta = min(gamma1 * delta_in, dnorm)  ! Powell's BOBYQA.

            elseif ratio <= eta2
                delta = max(gamma1 * delta_in, dnorm); % Powell's UOBYQA/NEWUOA/BOBYQA/LINCOA

            else
                delta = max(gamma1 * delta_in, gamma2 * dnorm); % Powell's NEWUOA/BOBYQA.
                %delta = max(delta_in, 1.25_RP * dnorm, dnorm + rho)  ! Powell's UOBYQA
                %delta = max(delta_in, gamma2 * dnorm)  ! Modified version. Works well for UOBYQA.
                % Powell's LINCOA code is as follows.
                %delta = min(max(gamma1 * delta_in, gamma2 * dnorm), sqrt(gamma2) * delta_in)
            end

            % For noisy problems, the following may work better.
            % %if (ratio <= eta1) then
            % %    delta = gamma1 * dnorm
            % %elseif (ratio <= eta2) then  ! Ensure DELTA >= DELTA_IN
            % %    delta = delta_in
            % %else  ! Ensure DELTA > DELTA_IN with a constant factor
            % %    delta = max(delta_in * (1.0_RP + gamma2) / 2.0_RP, gamma2 * dnorm)
            % %end if
            %

            %====================%
            %  Calculation ends  %
            %====================%

            % Postconditions
            if consts_obj.DEBUGGING
                debug_obj.assert(delta > 0, "DELTA > 0", srname);
            end

        end

    end
end