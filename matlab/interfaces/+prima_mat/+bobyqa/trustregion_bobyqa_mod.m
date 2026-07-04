classdef trustregion_bobyqa_mod
    %--------------------------------------------------------------------------------------------------%
    % This module provides subroutines concerning the trust-region calculations of BOBYQA.
    %
    % Coded by Zaikun ZHANG (www.zhangzk.net) based on Powell's code and the BOBYQA paper.
    %
    % Dedicated to the late Professor M. J. D. Powell FRS (1936--2015).
    %
    % Started: February 2022
    %
    % Last Modified: Thursday, April 04, 2024 PM09:26:23
    %--------------------------------------------------------------------------------------------------%

    methods
        function [crvmin, d] = trsbox(obj, delta, gopt_in, hq_in, pq_in, sl, su, tol, xopt, xpt, d)
            %--------------------------------------------------------------------------------------------------%
            % This subroutine approximately solves
            % minimize Q(XOPT + D) subject to ||D|| <= DELTA, SL <= XOPT + D <= SU.
            % See Section 3 of the BOBYQA paper.
            %
            % A version of the truncated conjugate gradient is applied. If a line search is restricted by
            % a constraint, then the procedure is restarted, the values of the variables that are at their
            % bounds being fixed. If the trust region boundary is reached, then further changes may be made to
            % D, each one being in the two dimensional space that is spanned by the current D and the gradient
            % of Q at XOPT+D, staying on the trust region boundary. Termination occurs when the reduction in
            % Q seems to be close to the greatest reduction that can be achieved.
            %
            % CRVMIN is set to zero if D reaches the trust region boundary. Otherwise it is set to the least
            % curvature of H that occurs in the conjugate gradient searches that are not restricted by any
            % constraints. In Powell's BOBYQA code, a negative value (-1) is assigned to CRVMIN if all of
            % the conjugate gradient searches are constrained. However, we set CRVMIN = 0 in this case, which
            % makes no difference to the algorithm, because CRVMIN is only used to tell whether the recent
            % models are sufficiently accurate, where both CRVMIN = 0 and CRVMIN < 0 provide a negative answer.
            %
            % XPT, XOPT, GOPT, HQ, PQ, SL and SU have the same meanings as the corresponding arguments of BOBYQB.
            % DELTA is the trust region radius for the present calculation
            % XNEW will be set to a new vector of variables that is approximately the one that minimizes the
            %   quadratic model within the trust region subject to the SL and SU constraints on the variables.
            %   It satisfies as equations the bounds that become active during the calculation.
            % GNEW holds the gradient of the quadratic model at XOPT+D. It is updated when D is updated.
            % XBDI is a working space vector. For I=1,2,...,N, the element XBDI(I) is set to -1.0, 0.0, or 1.0,
            %   the value being nonzero if and only if the I-th variable has become fixed at a bound, the bound
            %   being SL(I) or SU(I) in the case XBDI(I)=-1.0 or XBDI(I)=1.0, respectively. This information is
            %   accumulated during the construction of XNEW.
            % The arrays S and HS hold the current search direction and the change in the gradient of Q along S.
            %--------------------------------------------------------------------------------------------------%

            % Common modules
            consts_obj = prima_mat.common.consts_mod();
            debug_obj = prima_mat.common.debug_mod();
            infnan_obj = prima_mat.common.infnan_mod();
            linalg_obj = prima_mat.common.linalg_mod();
            powalg_obj = prima_mat.common.powalg_mod();
            univar_obj = prima_mat.common.univar_mod();


            % Inputs

            % GOPT_IN(N)
            % HQ_IN(N, N)
            % PQ_IN(NPT)
            % SL(N)
            % SU(N)

            % XOPT(N)
            % XPT(N, NPT)

            % Outputs

            % D(N)

            % Local variables
            srname = "TRSBOX";
            iact = NaN;
            n = NaN;
            npt = NaN;
            xbdi = NaN(numel(gopt_in), 1);
            grid_size = NaN;
            iter = NaN;
            itercg = NaN;
            maxiter = NaN;
            nact = NaN;
            nactsav = NaN;
            scaled = false;
            twod_search = false;
            beta = NaN;
            bstep = NaN;
            cth = NaN;
            delsq = NaN;
            dhd = NaN;
            dhs = NaN;
            dold = NaN(numel(d), 1);
            dredg = NaN;
            dredsq = NaN;
            ds = NaN;
            ggsav = NaN;
            gredsq = NaN;
            hangt = NaN;
            hangt_bd = NaN;
            hq = NaN(size(hq_in, 1), size(hq_in, 2));
            pq = NaN(numel(pq_in), 1);
            qred = NaN;
            rayleighq = NaN;
            resid = NaN;
            sbound = NaN(numel(gopt_in), 1);
            sdec = NaN;
            shs = NaN;
            sqrtd = NaN;
            sredg = NaN;
            stepsq = NaN;
            sth = NaN;
            stplen = NaN;
            temp = NaN;
            xtest = NaN(numel(xopt), 1);
            args = NaN(5, 1);
            dred = NaN(numel(gopt_in), 1);
            gnew = NaN(numel(gopt_in), 1);
            gopt = NaN(numel(gopt_in), 1);
            hdred = NaN(numel(gopt_in), 1);
            hs = NaN(numel(gopt_in), 1);
            modscal = NaN;
            s = NaN(numel(gopt_in), 1);
            sqdscr = NaN(numel(gopt_in), 1);
            ssq = NaN(numel(gopt_in), 1);
            tanbd = NaN(numel(gopt_in), 1);
            xnew = NaN(numel(gopt_in), 1);

            % Sizes
            n = fix(numel(gopt_in));
            npt = fix(numel(pq_in));

            % Preconditions
            if consts_obj.DEBUGGING
                debug_obj.assert(n >= 1 && npt >= n + 2, "N >= 1, NPT >= N + 2", srname);
                debug_obj.assert(delta > 0, "DELTA > 0", srname);
                debug_obj.assert(size(hq_in, 1) == n && linalg_obj.issymmetric(hq_in), "HQ is n-by-n and symmetric", srname);
                debug_obj.assert(numel(pq_in) == npt, "SIZE(PQ) == NPT", srname);
                debug_obj.assert(numel(sl) == n && all(sl <= 0, 'all'), "SIZE(SL) == N, SL <= 0", srname);
                debug_obj.assert(numel(su) == n && all(su >= 0, 'all'), "SIZE(SU) == N, SU >= 0", srname);
                debug_obj.assert(numel(xopt) == n && all(infnan_obj.is_finite(xopt), 'all'), "SIZE(XOPT) == N, XOPT is finite", srname);
                debug_obj.assert(all(xopt >= sl & xopt <= su, 'all'), "SL <= XOPT <= SU", srname);
                debug_obj.assert(size(xpt, 1) == n && size(xpt, 2) == npt, "SIZE(XPT) == [N, NPT]", srname);
                debug_obj.assert(all(infnan_obj.is_finite(xpt), 'all'), "XPT is finite", srname);
                debug_obj.assert(all(xpt >= fortran.spread(sl, 'dim', 2, 'ncopies', npt), 'all') && all(xpt <= fortran.spread(su, 'dim', 2, 'ncopies', npt), 'all'), "SL <= XPT <= SU", srname);
                debug_obj.assert(numel(d) == n, "SIZE(D) == N", srname);
                debug_obj.assert(numel(gnew) == n, "SIZE(GNEW) == N", srname);
                debug_obj.assert(numel(xnew) == n, "SIZE(XNEW) == N", srname);
            end

            %====================%
            % Calculation starts %
            %====================%

            % Scale the problem if GOPT contains large values. Otherwise, floating point exceptions may occur.
            % Note that CRVMIN must be scaled back if it is nonzero, but step is scale invariant.
            % N.B.: It is faster and safer to scale by multiplying a reciprocal than by division. See
            % https://fortran-lang.discourse.group/t/ifort-ifort-2021-8-0-1-0e-37-1-0e-38-0/
            if max(abs(gopt_in), [], 'all') > 1.0e12
                % The threshold is empirical.
                modscal = max(consts_obj.TWO * consts_obj.REALMIN, consts_obj.ONE / max(abs(gopt_in), [], 'all')); % MAX: precaution against underflow.
                gopt(:) = gopt_in * modscal;
                pq(:) = pq_in * modscal;
                hq(:, :) = hq_in * modscal;
                scaled = true;
            else
                modscal = consts_obj.ONE; % This value is not used, but Fortran compilers may complain without it.
                gopt(:) = gopt_in;
                pq(:) = pq_in;
                hq(:, :) = hq_in;
                scaled = false;
            end

            % The initial values of IACT, DREDSQ, and GGSAV are unused but to entertain Fortran compilers.
            % TODO: Check that GGSAV has been initialized before used.
            iact = 0;
            dredsq = consts_obj.ZERO;
            ggsav = consts_obj.ZERO;

            % The sign of GOPT(I) gives the sign of the change to the I-th variable that will reduce Q from its
            % value at XOPT. Thus XBDI(I) shows whether or not to fix the I-th variable at one of its bounds
            % initially, with NACT being set to the number of fixed variables.
            xbdi(:) = 0;
            xbdi(linalg_obj.trueloc(xopt >= su & gopt <= 0)) = 1;
            xbdi(linalg_obj.trueloc(xopt <= sl & gopt >= 0)) = -1;
            nact = fix(nnz(xbdi ~= 0));

            % Initialized D and CRVMIN.
            d(:) = consts_obj.ZERO;
            crvmin = -consts_obj.REALMAX;

            % GNEW is the gradient at the current iterate.
            gnew(:) = gopt;
            gredsq = sum(fortran.power(gnew(linalg_obj.trueloc(xbdi == 0)), 2), 'all');
            % DELSQ is the upper bound on the sum of squares of the free variables.
            delsq = delta * delta;
            % QRED is the reduction in Q so far.
            qred = consts_obj.ZERO;
            % BETA is the coefficient for the previous searching direction in the conjugate gradient method.
            beta = consts_obj.ZERO;

            % ITERCG is the number of CG iterations corresponding to the current set of active bounds.
            itercg = 0;

            % TWOD_SEARCH: whether to perform a 2-dimensional search after the truncated CG method.
            twod_search = false; % The default value of TWOD_SEARCH is FALSE!

            % Powell's code is essentially a DO WHILE loop. We impose an explicit MAXITER.
            % The formulation of MAXITER below contains a precaution against overflow. In MATLAB/Python/Julia/R,
            % we can write maxiter = min(10000, (n - nact)^2).
            % Powell commented in the BOBYQA paper (the paragraph above (3.7)) that "numerical experiments show
            % that it is very unusual for subroutine TRSBOX to make more than ten changes to d when seeking an
            % approximate solution to the subproblem (1.8), even if there are hundreds of variables."
            maxiter = fix(min(fortran.power(10, min(4, floor(log10(double(intmax('int64')))))), fortran.power(fix(n - nact), 2)));
            for iter = 1:maxiter
                resid = delsq - sum(fortran.power(d(linalg_obj.trueloc(xbdi == 0)), 2), 'all');
                if resid <= 0
                    twod_search = true;
                    break
                end

                % Set the next search direction of the conjugate gradient method. It is the steepest descent
                % direction initially and when the iterations are restarted because a variable has just been
                % fixed by a bound, and of course the components of the fixed variables are zero. MAXITER is an
                % upper bound on the indices of the conjugate gradient iterations.
                if itercg == 0
                    % TODO: If we are sure that S contain only finite values, we may merge this case into the next.
                    s(:) = -gnew;
                else
                    s(:) = beta * s - gnew;
                end
                s(linalg_obj.trueloc(xbdi ~= 0)) = consts_obj.ZERO;
                stepsq = sum(fortran.power(s, 2), 'all');
                ds = linalg_obj.inprod(d(linalg_obj.trueloc(xbdi == 0)), s(linalg_obj.trueloc(xbdi == 0)));

                if ~(stepsq > consts_obj.EPS * delsq && gredsq * delsq > fortran.power((tol * qred), 2) && ~infnan_obj.is_nan_sp(ds))
                    break
                end

                % Set BSTEP to the length of the step to the trust region boundary and STPLEN to the steplength,
                % ignoring the simple bounds.

                % SQRTD: square root of a discriminant. The MAXVAL avoids SQRTD < ABS(DS) due to underflow.
                sqrtd = max([fortran.sqrt(stepsq * resid + ds * ds), fortran.sqrt(stepsq * resid), abs(ds)], [], 'all');

                % Zaikun 20220210: For the IF ... ELSE ... END IF below, Powell's condition for the IF is DS>=0.
                % In theory, switching the condition to DS > 0 changes nothing; indeed, the two formulations
                % of BSTEP are equivalent. However, surprisingly, DS > 0 clearly worsens the performance of
                % BOBYQA in tests on 20220210, 20221206. Why? When DS = 0, what should be the best formulation?
                % What if we are at the first iteration? BSTEP = DELTA/||D||?
                % See TRSAPP.F90 of NEWUOA.
                %if (ds > 0) then  ! Zaikun 20210925
                if ds >= 0
                    bstep = resid / (sqrtd + ds);
                else
                    bstep = (sqrtd - ds) / stepsq;
                end
                % BSTEP < 0 should not happen. BSTEP can be 0 or NaN when, e.g., DS or STEPSQ becomes Inf.
                % Powell's code does not handle this.
                if bstep <= 0 || ~infnan_obj.is_finite(bstep)
                    break
                end

                hs(:) = powalg_obj.hess_mul(s, xpt, pq, 'hq', hq);
                shs = linalg_obj.inprod(s(linalg_obj.trueloc(xbdi == 0)), hs(linalg_obj.trueloc(xbdi == 0)));
                stplen = bstep;
                if shs > 0
                    stplen = min(bstep, gredsq / shs);
                end

                % Reduce STPLEN if necessary in order to preserve the simple bounds, letting IACT be the index
                % of the new constrained variable.
                % N.B. (Zaikun 20220422):
                % Theory and computation differ considerably in the calculation of STPLEN and IACT.
                % 1. Theoretically, the WHERE constructs can simplify (S > 0 .and. XTEST > SU) to (S > 0) and
                % (S < 0, XTEST < SL) to (S < 0), which will be equivalent to Powell's original code. However,
                % overflow will occur due to huge values in SU or SL that indicate the absence of bounds, and
                % Fortran compilers will complain. It is not an issue in MATLAB/Python/Julia/R.
                % 2. Theoretically, we can also simplify (S > 0 .and. XTEST > SU) to (XTEST > SU). This is
                % because the algorithm intends to ensure that SL <= XSUM <= SU, under which the inequality
                % XTEST(I) > SU(I) implies S(I) > 0. Numerically, however, XSUM may violate the bounds slightly
                % due to rounding. If we replace (S > 0 .and. XTEST > SU) with (XTEST > SU), then SBOUND(I) will
                % be -Inf when SU(I) - XSUM(I) is negative (although tiny) and S(I) is +0 (positively signed
                % zero), which will lead to STPLEN = -Inf and IACT = I > 0. This will trigger a restart of the
                % conjugate gradient method with DELSQ updated to DELSQ - D(IACT)**2; if D(IACT)**2 << DELSQ,
                % then DELSQ can remain unchanged due to rounding, leading to an infinite cycling.
                % 3. Theoretically, the WHERE construct corresponding to S > 0 can calculate SBOUND by
                % MIN(STPLEN * S, SU - XSUM) / S instead of (SU - XSUM) / S, since this quotient matters only if
                % it is less than STPLEN. The motivation is to avoid overflow even without checking XTEST > XU.
                % Yet such an implementation clearly worsens the performance of BOBYQA in our test on 20220422.
                % Why? Note that the conjugate gradient method restarts when IACT > 0. Due to rounding errors,
                % MIN(STPLEN * S, SU - XSUM) / S can frequently contain entries less than STPLEN, leading to a
                % positive IACT and hence a restart. This turns out harmful to the performance of the algorithm,
                % but WHY? It can be rectified in two ways: use MIN(STPLEN, (SU-XSUM) / S) instead of
                % MIN(STPLEN*S, SU-XSUM)/S, or set IACT to a positive value only if the minimum of SBOUND is
                % surely less STPLEN, e.g. ANY(SBOUND < (ONE-EPS) * STPLEN). The first method does not avoid
                % overflow and makes little sense.
                xnew(:) = xopt + d;
                xtest(:) = xnew + stplen * s;
                sbound(:) = stplen;
                sbound(s > 0 & xtest > su) = (su(s > 0 & xtest > su) - xnew(s > 0 & xtest > su)) ./ s(s > 0 & xtest > su);
                sbound(s < 0 & xtest < sl) = (sl(s < 0 & xtest < sl) - xnew(s < 0 & xtest < sl)) ./ s(s < 0 & xtest < sl);
                %%MATLAB:
                %%sbound(s > 0) = (su(s > 0) - xnew(s > 0)) / s(s > 0);
                %%sbound(s < 0) = (sl(s < 0) - xnew(s < 0)) / s(s < 0);
                %----------------------------------------------------------------------------------------------%
                % The code below is mathematically equivalent to the above but numerically inferior as explained.
                %where (s > 0) sbound = min(stplen * s, su - xnew) / s
                %where (s < 0) sbound = max(stplen * s, sl - xnew) / s
                %----------------------------------------------------------------------------------------------%
                sbound(linalg_obj.trueloc(infnan_obj.is_nan(sbound))) = stplen; % Needed? No if we are sure that D and S are finite.
                iact = 0;
                if any(sbound < stplen, 'all')
                    iact = fix(fortran.minloc(sbound, 'dim', 1));
                    stplen = sbound(iact);
                    %%MATLAB: [stplen, iact] = min(sbound);

                end
                %----------------------------------------------------------------------------------------------%
                % Alternatively, IACT and STPLEN can be calculated as below.
                % %IACT = INT(MINLOC([STPLEN, SBOUND], DIM=1), KIND(IACT)) - 1_IK
                % %STPLEN = MINVAL([STPLEN, SBOUND]) ! This line cannot be exchanged with the last
                % We prefer our implementation, as the code is more explicit; in addition, it is more flexible:
                % we can change the condition ANY(SBOUND < STPLEN) to ANY(SBOUND < (1 - EPS) * STPLEN) or
                % ANY(SBOUND < (1 + EPS) * STPLEN), depending on whether we believe a false positive or a false
                % negative of IACT > 0 is more harmful --- according to our test on 20220422, it is the former,
                % as mentioned above.
                %----------------------------------------------------------------------------------------------%

                % Update CRVMIN, GNEW, and D. Set SDEC to the decrease that occurs in Q.
                sdec = consts_obj.ZERO;
                if stplen > 0
                    itercg = itercg + 1;
                    rayleighq = shs / stepsq;
                    if iact == 0 && rayleighq > 0
                        if crvmin <= -consts_obj.REALMAX
                            % CRVMIN <= -REALMAX means CRVMIN has not been set.
                            crvmin = rayleighq;
                        else
                            crvmin = min(crvmin, rayleighq);
                        end
                    end
                    ggsav = gredsq;
                    gnew(:) = gnew + stplen * hs;
                    gredsq = sum(fortran.power(gnew(linalg_obj.trueloc(xbdi == 0)), 2), 'all');
                    dold(:) = d;
                    d(:) = d + stplen * s;

                    % Exit in case of Inf/NaN in D.
                    if ~infnan_obj.is_finite(sum(abs(d), 'all'))
                        d(:) = dold;
                        break
                    end

                    sdec = max(stplen * (ggsav - consts_obj.HALF * stplen * shs), consts_obj.ZERO);
                    qred = qred + sdec;
                end

                % Restart the conjugate gradient method if it has hit a new bound.
                if iact > 0
                    nact = nact + 1;
                    debug_obj.assert(abs(s(iact)) > 0, "S(IACT) /= 0", srname);
                    xbdi(iact) = round(fortran.sign(consts_obj.ONE, s(iact))); %%MATLAB: xbdi(iact) = sign(s(iact))
                    % Exit when NACT = N (NACT > N is impossible). We must update XBDI before exiting!
                    if nact >= n
                        break % This leads to a difference. Why?

                    end
                    delsq = delsq - fortran.power(d(iact), 2);
                    if delsq <= 0
                        twod_search = true;
                        % Why set TWOD_SEARCH to TRUE? Because DELSQ <= 0 just means that D reaches the trust
                        % region boundary.
                        break
                    end
                    beta = consts_obj.ZERO;
                    itercg = 0;
                    gredsq = sum(fortran.power(gnew(linalg_obj.trueloc(xbdi == 0)), 2), 'all');
                elseif stplen < bstep
                    % Either apply another conjugate gradient iteration or exit.
                    % N.B. ITERCG > N - NACT is impossible.
                    if itercg >= n - nact || sdec <= tol * qred || infnan_obj.is_nan_sp(sdec) || infnan_obj.is_nan_sp(qred)
                        break
                    end
                    beta = gredsq / ggsav; % Has GGSAV got the correct value yet?

                else
                    twod_search = true;
                    break
                end
            end

            % Set MAXITER for the 2-dimensional search on the trust region boundary. Powell's code essentially
            % sets MAXITER to infinity; the loop exits when NACT >= N-1 or the procedure cannot significantly
            % reduce the quadratic model. We set a finite but large MAXITER as a safeguard.
            if twod_search
                crvmin = consts_obj.ZERO;
                maxiter = 10 * (n - nact);
            else
                maxiter = 0;
            end

            % Improve D by a sequential 2-dimensional search on the boundary of the trust region for the
            % variables that have not reached a bound. See (3.6) of the BOBYQA paper and the elaborations nearby.
            % 1. At each iteration, the current D is improved by a search conducted on the circular arch
            % {D(THETA): D(THETA) = (I-P)*D + [COS(THETA)*P*D + SIN(THETA)*S], 0<=THETA<=PI/2, SL<=XOPT+D(THETA)<=SU},
            % where P is the orthogonal projection onto the space of the variables that have not reached their
            % bounds, and S is a linear combination of P*D and P*G(XOPT+D) with ||S|| = ||P*D|| and G(.) being
            % the gradient of the quadratic model. The iteration is performed only if P*D and P*G(XOPT+D) are
            % not nearly parallel. The arc lies in the hyperplane (I-P)*D + Span{P*D, P*G(XOPT+D)} and the trust
            % region boundary {D: ||D||=DELTA}; it is part of the circle (I-P)*D + {COS(THETA)*P*D + SIN(THETA)*S}
            % with THETA being in [0, PI/2] and restricted by the bounds on X.
            % 2. In (3.6) of the BOBYQA paper, Powell wrote that 0 <= THETA <= PI/4, which seems a typo.
            % 3. The search on the arch is done by calling INTERVAL_MAX, which maximizes INTERVAL_FUN_TRSBOX.
            % INTERVAL_FUN_TRSBOX is essentially Q(XOPT + D) - Q(XOPT + D(THETA)), but its independent variable
            % is not THETA but TAN(THETA/2), namely "tangent of the half angle" in Powell's code/comments. This
            % "half" may be the reason for the apparent typo mentioned above.
            % Question (Zaikun 20220424): Shouldn't we try something similar in GEOSTEP?

            nactsav = nact - 1;
            for iter = 1:maxiter
                xnew(:) = xopt + d;

                % Update XBDI. It indicates whether the lower (-1) or upper bound (+1) is reached or not (0).
                xbdi(linalg_obj.trueloc(xbdi == 0 & (xnew >= su))) = 1;
                xbdi(linalg_obj.trueloc(xbdi == 0 & (xnew <= sl))) = -1;
                nact = fix(nnz(xbdi ~= 0));
                if nact >= n - 1
                    break
                end

                % Update GREDSQ, DREDG, DREDSQ.
                gredsq = sum(fortran.power(gnew(linalg_obj.trueloc(xbdi == 0)), 2), 'all');
                dredg = linalg_obj.inprod(d(linalg_obj.trueloc(xbdi == 0)), gnew(linalg_obj.trueloc(xbdi == 0)));
                if iter == 1 || nact > nactsav
                    dredsq = sum(fortran.power(d(linalg_obj.trueloc(xbdi == 0)), 2), 'all'); % In theory, DREDSQ changes only when NACT increases.
                    dred(:) = d;
                    dred(linalg_obj.trueloc(xbdi ~= 0)) = consts_obj.ZERO;
                    hdred(:) = powalg_obj.hess_mul(dred, xpt, pq, 'hq', hq);
                    nactsav = nact;
                end

                % Let the search direction S be a linear combination of the reduced D and the reduced G that is
                % orthogonal to the reduced D.
                temp = gredsq * dredsq - dredg * dredg;
                if ~(temp > fortran.power(tol, 2) * max(gredsq * dredsq, fortran.power(qred, 2)))
                    % TEMP is tiny or NaN occurs
                    break
                end
                temp = fortran.sqrt(temp);
                s(:) = (dredg * d - dredsq * gnew) ./ temp;
                s(linalg_obj.trueloc(xbdi ~= 0)) = consts_obj.ZERO;
                sredg = -temp;

                % By considering the simple bounds on the free variables, calculate an upper bound on the
                % TANGENT of HALF the angle of the alternative iteration, namely ANGBD. The bounds are
                % SL - XOPT <= COS(THETA)*D + SIN(THETA)*S <= SU - XOPT for the free variables.
                % Defining HANGT = TAN(THETA/2), and using the tangent half-angle formula, we have
                % (1+HANGT^2)*(SL - XOPT) <= (1-HANGT^2)*D + 2*HANGT*S <= (1+HANGT^2)*(SU - XOPT),
                % which is required for all free variables. The indices of the free variable are those with
                % XBDI == 0. Solving this inequality system for HANGT in [0, PI/4], we get bounds for HANGT,
                % namely TANBD; the final bound for HANGT is the minimum of TANBD, which is HANGT_BD.
                % When solving the system, note that SL < XOPT < SU and SL < XOPT + D < SU if XBDI = 0.
                %
                % Note the following for the calculation of the first SQDSCR below (the second is similar).
                % 0. SQDSCR means "square root of discriminant".
                % 1. When calculating the first SQDSCR, Powell's code checks whether SSQ - (XOPT - SL)**2) is
                % positive. However, overflow will occur if SL contains large values that indicate absence of
                % bounds. It is not a problem in MATLAB/Python/Julia/R.
                % 2. Even if XOPT - SL < SQRT(SSQ), rounding errors may render SSQ - (XOPT - SL)**2) < 0.
                ssq(:) = fortran.power(d, 2) + fortran.power(s, 2); % Indeed, only SSQ(TRUELOC(XBDI == 0)) is needed.
                tanbd(:) = consts_obj.ONE;
                sqdscr(:) = -consts_obj.REALMAX;
                sqdscr(xbdi == 0 & xopt - sl < sqrt(ssq)) = fortran.sqrt(max(consts_obj.ZERO, ssq(xbdi == 0 & xopt - sl < fortran.sqrt(ssq)) - fortran.power((xopt(xbdi == 0 & xopt - sl < fortran.sqrt(ssq)) - sl(xbdi == 0 & xopt - sl < fortran.sqrt(ssq))), 2)));
                tanbd(sqdscr - s > 0) = min(tanbd(sqdscr - s > 0), (xnew(sqdscr - s > 0) - sl(sqdscr - s > 0)) ./ (sqdscr(sqdscr - s > 0) - s(sqdscr - s > 0)));
                sqdscr(:) = -consts_obj.REALMAX;
                sqdscr(xbdi == 0 & su - xopt < sqrt(ssq)) = fortran.sqrt(max(consts_obj.ZERO, ssq(xbdi == 0 & su - xopt < fortran.sqrt(ssq)) - fortran.power((su(xbdi == 0 & su - xopt < fortran.sqrt(ssq)) - xopt(xbdi == 0 & su - xopt < fortran.sqrt(ssq))), 2)));
                tanbd(sqdscr + s > 0) = min(tanbd(sqdscr + s > 0), (su(sqdscr + s > 0) - xnew(sqdscr + s > 0)) ./ (sqdscr(sqdscr + s > 0) + s(sqdscr + s > 0)));
                tanbd(linalg_obj.trueloc(infnan_obj.is_nan(tanbd))) = consts_obj.ZERO;
                %----------------------------------------------------------------------------------------------%
                %%MATLAB code for defining TANBD:
                %%xfree = (xbdi == 0);
                %%ssq = NaN(n, 1);
                %%ssq(xfree) = s(xfree).^2 + d(xfree).^2;
                %%discmn = NaN(n, 1);
                %%discmn(xfree) = ssq(xfree) - (xopt(xfree) - sl(xfree))**2;  % This is a discriminant.
                %%tanbd = 1;
                %%mask = (xfree & discmn > 0 & sqrt(discmn) - s > 0);
                %%tanbd(mask) = min(tanbd(mask), (xnew(mask) - sl(mask)) / (sqrt(discmn(mask)) - s(mask)));
                %%discmn(xfree) = ssq(xfree) - (su(xfree) - xopt(xfree))**2;  % This is a discriminant.
                %%mask = (xfree & discmn > 0 & sqrt(discmn) + s > 0);
                %%tanbd(mask) = min(tanbd(mask), (su(mask) - xnew(mask)) / (sqrt(discmn(mask)) + s(mask)));
                %%tanbd(isnan(tanbd)) = 0;
                %----------------------------------------------------------------------------------------------%

                iact = 0;
                hangt_bd = consts_obj.ONE;
                if any(tanbd < 1, 'all')
                    iact = fix(fortran.minloc(tanbd, 'dim', 1));
                    hangt_bd = tanbd(iact);
                    %%MATLAB: [hangt_bd, iact] = min(tanbd);

                end
                if hangt_bd <= 0
                    break
                end

                % Calculate HS and some curvatures for the alternative iteration.
                hs(:) = powalg_obj.hess_mul(s, xpt, pq, 'hq', hq);
                shs = linalg_obj.inprod(s(linalg_obj.trueloc(xbdi == 0)), hs(linalg_obj.trueloc(xbdi == 0)));
                dhs = linalg_obj.inprod(d(linalg_obj.trueloc(xbdi == 0)), hs(linalg_obj.trueloc(xbdi == 0)));
                dhd = linalg_obj.inprod(d(linalg_obj.trueloc(xbdi == 0)), hdred(linalg_obj.trueloc(xbdi == 0)));

                % Seek the greatest reduction in Q for a range of equally spaced values of HANGT in [0, ANGBD],
                % with HANGT being the TANGENT of HALF the angle of the alternative iteration.
                args(:) = [shs, dhd, dhs, dredg, sredg];
                if any(infnan_obj.is_nan(args), 'all')
                    break
                end
                % Define the grid size of the search for HANGT. Powell defined the size to be 4 if hangt_bd is
                % nearly zero and 20 if it is nearly one, with a linear interpolation in between. We double this
                % size, which improves the performance of BOBYQA in general according to a test on 20230827.
                %grid_size = nint(17.0_RP * hangt_bd + 4.1_RP, kind(grid_size))  ! Powell's version
                grid_size = 2 * round(17.0 * hangt_bd + 4.1);
                %%MATLAB: grid_size = 2 * round(17 * hangt_bd + 4.1_RP)
                hangt = univar_obj.interval_max(@(varargin) obj.interval_fun_trsbox(varargin{:}), consts_obj.ZERO, hangt_bd, args, grid_size);
                sdec = obj.interval_fun_trsbox(hangt, args);
                if ~(sdec > 0)
                    break
                end

                % Update GNEW, D and HDRED. If the angle of the alternative iteration is restricted by a bound
                % on a free variable, that variable is fixed at the bound. The MIN below is a precaution against
                % rounding errors.
                cth = min((consts_obj.ONE - fortran.power(hangt, 2)) / (consts_obj.ONE + fortran.power(hangt, 2)), consts_obj.ONE - fortran.power(hangt, 2));
                sth = min((hangt + hangt) / (consts_obj.ONE + fortran.power(hangt, 2)), hangt + hangt);
                gnew(:) = gnew + (cth - consts_obj.ONE) * hdred + sth * hs;
                dold(:) = d;
                d(linalg_obj.trueloc(xbdi == 0)) = cth * d(linalg_obj.trueloc(xbdi == 0)) + sth * s(linalg_obj.trueloc(xbdi == 0));

                % Exit in case of Inf/NaN in D.
                if ~infnan_obj.is_finite(sum(abs(d), 'all'))
                    d(:) = dold;
                    break
                end

                hdred(:) = cth * hdred + sth * hs;
                qred = qred + sdec;
                if iact >= 1 && iact <= n && hangt >= hangt_bd
                    % D(IACT) reaches lower/upper bound.
                    xbdi(iact) = round(fortran.sign(consts_obj.ONE, xopt(iact) + d(iact) - consts_obj.HALF * (sl(iact) + su(iact))));
                    %%MATLAB: xbdi(iact) = sign(xopt(iact)+d(iact) - 0.5*(sl+su));

                elseif ~(sdec > tol * qred)
                    % SDEC is small or NaN occurs
                    break
                end
            end

            % Set D, giving careful attention to the bounds.
            xnew(:) = max(sl, min(su, xopt + d));
            xnew(linalg_obj.trueloc(xbdi == -1)) = sl(linalg_obj.trueloc(xbdi == -1));
            xnew(linalg_obj.trueloc(xbdi == 1)) = su(linalg_obj.trueloc(xbdi == 1));
            d(:) = xnew - xopt;

            % Set CRVMIN to ZERO if it has never been set or becomes NaN due to ill conditioning.
            if crvmin <= -consts_obj.REALMAX || infnan_obj.is_nan_sp(crvmin)
                crvmin = consts_obj.ZERO;
            end

            % Scale CRVMIN back before return. Note that the trust-region step is scale invariant.
            if scaled && crvmin > 0
                crvmin = crvmin / modscal;
            end

            %====================%
            %  Calculation ends  %
            %====================%

            % Postconditions
            if consts_obj.DEBUGGING
                debug_obj.assert(numel(d) == n && all(infnan_obj.is_finite(d), 'all'), "SIZE(D) == N, D is finite", srname);
                % Due to rounding, it may happen that ||D|| > DELTA, but ||D|| > 2*DELTA is highly improbable.
                debug_obj.assert(linalg_obj.p_norm(d) <= consts_obj.TWO * delta, "||D|| <= 2*DELTA", srname);
                debug_obj.assert(crvmin >= 0, "CRVMIN >= 0", srname);
                % D is supposed to satisfy the bound constraints SL <= XOPT + D <= SU.
                debug_obj.assert(all(xopt + d >= sl - consts_obj.TEN * consts_obj.EPS * max(consts_obj.ONE, abs(sl)) & xopt + d <= su + consts_obj.TEN * consts_obj.EPS * max(consts_obj.ONE, abs(su)), 'all'), "SL <= XOPT + D <= SU", srname);
            end

        end
        function f = interval_fun_trsbox(~, hangt, args)
            %--------------------------------------------------------------------------------------------------%
            % This function defines the objective function of the search for HANGT in TRSBOX, with HANGT being
            % the TANGENT of HALF the angle of the "alternative iteration".
            %--------------------------------------------------------------------------------------------------%
            consts_obj = prima_mat.common.consts_mod();
            debug_obj = prima_mat.common.debug_mod();

            % Inputs



            % Outputs
            f = NaN;

            % Local variables
            srname = "INTERVAL_FUN_TRSBOX";
            sth = NaN;

            % Preconditions
            if consts_obj.DEBUGGING
                debug_obj.assert(numel(args) == 5, "SIZE(ARGS) == 5", srname);
            end

            %====================%
            % Calculation starts %
            %====================%

            f = consts_obj.ZERO;
            if abs(hangt) > 0
                sth = (hangt + hangt) / (consts_obj.ONE + hangt * hangt);
                f = args(1) + hangt * (hangt * args(2) - args(3) - args(3));
                f = sth * (hangt * args(4) - args(5) - consts_obj.HALF * sth * f);
                % N.B.: ARGS = [SHS, DHD, DHS, DREDG, SREDG]

            end

            %====================%
            %  Calculation ends  %
            %====================%
        end
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
                delta = min(gamma1 * delta_in, dnorm); % Powell's BOBYQA.
                %delta = gamma1 * dnorm  ! Powell's UOBYQA/NEWUOA.
                %delta = gamma1 * delta_in  ! Powell's COBYLA/LINCOA. Works poorly here.

            elseif ratio <= eta2
                delta = max(gamma1 * delta_in, dnorm); % Powell's UOBYQA/NEWUOA/BOBYQA/LINCOA

            else
                delta = max(gamma1 * delta_in, gamma2 * dnorm); % Powell's NEWUOA/BOBYQA.
                %delta = max(delta_in, gamma2 * dnorm)  ! Modified version. Works well for UOBYQA.
                %delta = max(delta_in, 1.25_RP * dnorm, dnorm + rho)  ! Powell's UOBYQA
                %delta = min(max(gamma1 * delta_in, gamma2 * dnorm), sqrt(gamma2) * delta_in)  ! Powell's LINCOA.
            end

            % For noisy problems, the following may work better.
            % %if (ratio <= eta1) then
            % %    delta = gamma1 * dnorm
            % %elseif (ratio <= eta2) then  ! Ensure DELTA >= DELTA_IN
            % %    delta = delta_in
            % %else  ! Ensure DELTA > DELTA_IN with a constant factor
            % %    delta = max(delta_in * (1.0_RP + gamma2) / 2.0_RP, gamma2 * dnorm)
            % %end if

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