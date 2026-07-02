classdef trustregion_newuoa_mod
    %--------------------------------------------------------------------------------------------------%
    % This module provides subroutines concerning the trust-region calculations of NEWUOA.
    %
    % Coded by Zaikun ZHANG (www.zhangzk.net) based on Powell's code and the NEWUOA paper.
    %
    % Dedicated to the late Professor M. J. D. Powell FRS (1936--2015).
    %
    % Started: July 2020
    %
    % Last Modified: Saturday, April 06, 2024 PM11:03:07
    %--------------------------------------------------------------------------------------------------%

    methods
        function [crvmin, s, info_loc] = trsapp(obj, delta, gopt_in, hq_in, pq_in, tol, xpt, s, varargin)
            %--------------------------------------------------------------------------------------------------%
            % This subroutine finds an approximate solution to the N-dimensional trust region subproblem
            %
            % min <S, GOPT> + 0.5*<S, HESSIAN*S> s.t. ||S|| <= DELTA
            %
            % Note that the HESSIAN here is the sum of an explicit part HQ and an implicit part (PQ, XPT):
            %
            % HESSIAN = HQ + sum_K=1^NPT PQ(K)*XPT(:, K)*XPT(:, K)' .
            %
            % The calculation of S begins with the truncated conjugate gradient method. If the boundary of the
            % trust region is reached, then further changes to S may be made, each one being in the 2-dimensional
            % space spanned by the current S and the corresponding gradient of Q. Thus S should provide a
            % substantial reduction to Q within the trust region. See Section 5 of the NEWUOA paper.
            %
            % At return, S will be the approximate solution. CRVMIN will be set to the least curvature of
            % HESSIAN along the conjugate directions that occur, except that it is set to ZERO if S goes all the
            % way to the trust-region boundary. INFO is an exit flag with the following possible values.
            % - INFO = 0: an approximate solution satisfying one of the following conditions is found:
            % 1. ||G+HS||/||G0|| <= TOL,
            % 2. ||S|| = DELTA and <S, -(G+HS)> >= (1 - TOL)*||S||*||G+HS||,
            % where TOL is set to 1e-2 in NEWUOA;
            % - INFO = 1: the last iteration reduces Q only insignificantly;
            % - INFO = 2: the maximal number of iterations is attained;
            % - INFO = -1: too much rounding error to continue.
            %--------------------------------------------------------------------------------------------------%

            % Common modules
            consts_obj = consts_mod();
            debug_obj = debug_mod();
            infnan_obj = infnan_mod();
            linalg_obj = linalg_mod();
            powalg_obj = powalg_mod();
            univar_obj = univar_mod();


            % Inputs

            % GOPT_IN(N)
            % HQ_IN(N, N)
            % PQ_IN(NPT)

            % XPT(N, NPT)

            % Outputs

            % S(N)


            % Local variables
            srname = "TRSAPP";
            info_loc = NaN;
            iter = NaN;
            maxiter = NaN;
            n = NaN;
            npt = NaN;
            scaled = false;
            twod_search = false;
            alpha = NaN;
            angle = NaN;
            args = NaN(4, 1);
            bstep = NaN;
            cth = NaN;
            d = NaN(numel(gopt_in), 1);
            dd = NaN;
            delsq = NaN;
            dg = NaN;
            dhd = NaN;
            dhs = NaN;
            ds = NaN;
            g = NaN(numel(gopt_in), 1);
            gg = NaN;
            gg0 = NaN;
            ggsav = NaN;
            gopt = NaN(numel(gopt_in), 1);
            hd = NaN(numel(gopt_in), 1);
            hq = NaN(size(hq_in, 1), size(hq_in, 2));
            hs = NaN(numel(gopt_in), 1);
            modscal = NaN;
            pq = NaN(numel(pq_in), 1);
            qadd = NaN;
            qred = NaN;
            reduc = NaN;
            resid = NaN;
            sg = NaN;
            shs = NaN;
            sold = NaN(numel(gopt_in), 1);
            sqrtd = NaN;
            ss = NaN;
            sth = NaN;

            % Sizes
            n = size(xpt, 1);
            npt = size(xpt, 2);

            % Preconditions
            if consts_obj.DEBUGGING
                debug_obj.assert(n >= 1 && npt >= n + 2, "N >= 1, NPT >= N + 2", srname);
                debug_obj.assert(delta > 0, "DELTA > 0", srname);
                debug_obj.assert(numel(gopt_in) == n, "SIZE(GOPT) = N", srname);
                debug_obj.assert(size(hq_in, 1) == n && linalg_obj.issymmetric(hq_in), "HQ is an NxN symmetric matrix", srname);
                debug_obj.assert(numel(pq_in) == npt, "SIZE(PQ) = NPT", srname);
                debug_obj.assert(all(infnan_obj.is_finite(xpt), 'all'), "XPT is finite", srname);
                debug_obj.assert(numel(s) == n, "SIZE(S) == N", srname);
            end

            %====================%
            % Calculation starts %
            %====================%

            % Scale the problem if GOPT contains large values. Otherwise, floating point exceptions may occur.
            % Note that CRVMIN must be scaled back if it is nonzero, but the step is scale invariant.
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

            s(:) = consts_obj.ZERO;
            crvmin = consts_obj.ZERO;
            qred = consts_obj.ZERO;
            info_loc = 2; % Default exit flag is 2, i.e., MAXITER is attained

            % Prepare for the first line search.
            %--------------------------------------------------------------------------------------------------%
            % N.B.: During the iterations, G is NOT updated, and it equals always GOPT, which is the gradient
            % of the trust-region model at the trust-region center X. However, GG is updated: GG = ||G + HS||^2,
            % which is the norm square of the gradient at the current iterate.
            g(:) = gopt;
            gg = linalg_obj.inprod(g, g);
            %--------------------------------------------------------------------------------------------------%
            gg0 = gg;
            d(:) = -g;
            dd = gg;
            ds = consts_obj.ZERO;
            ss = consts_obj.ZERO;
            hs(:) = consts_obj.ZERO;
            delsq = delta * delta;
            maxiter = n;

            twod_search = false;

            % The truncated-CG iterations.
            %
            % The iteration will be terminated in 4 possible cases:
            % 1. the maximal number of iterations is attained;
            % 2. QADD <= TOL*QRED or ||G|| <= TOL*||G0||, where QADD is the reduction of Q due to the latest
            % CG step, QRED is the reduction of Q since the beginning until the latest CG step, G is the
            % current gradient, and G0 is the initial gradient; see (5.13) of the NEWUOA paper;
            % 3. DS <= 0
            % 4. ||S|| = DELTA, i.e., CG path cuts the trust region boundary.
            %
            % In the 4th case, twod_search will be set to true, meaning that S will be improved by a sequence of
            % two-dimensional search, the two-dimensional subspace at each iteration being span(S, -G).
            for iter = 1:maxiter
                % Exit if G contains NaN.
                if infnan_obj.is_nan_sp(gg)
                    info_loc = -1;
                    break
                end
                % Exit if GG is small. This must be done first; otherwise, DD can be 0 and BSTEP is not well
                % defined. The inequality below must be non-strict so that GG = GG0 = 0 will trigger the exit.
                if gg <= (fortran.power(tol, 2)) * gg0
                    info_loc = 0;
                    break
                end

                % Set BSTEP to the step length such that ||S + BSTEP*D|| = DELTA.
                if iter == 1
                    bstep = delta / fortran.sqrt(dd);
                else
                    resid = delsq - ss;
                    if resid <= 0
                        twod_search = true;
                        break
                    end

                    % Powell's code does not have the following two IFs.
                    %--------------------------------------------------%
                    if dd <= consts_obj.EPS * delsq
                        info_loc = 0;
                        break
                    end
                    if infnan_obj.is_nan_sp(ds)
                        info_loc = -1;
                        break
                    end
                    %--------------------------------------------------%

                    % SQRTD: square root of a discriminant. The MAXVAL avoids SQRTD < ABS(DS) due to underflow.
                    sqrtd = max([fortran.sqrt(fortran.power(ds, 2) + dd * resid), abs(ds), fortran.sqrt(dd * resid)], [], 'all');
                    % Powell's code does not distinguish the following two cases, which have no difference in
                    % precise arithmetic. The following scheme stabilizes the calculation. Copied from LINCOA.
                    if ds <= 0
                        bstep = (sqrtd - ds) / dd;
                    else
                        bstep = resid / (sqrtd + ds);
                    end
                end

                % BSTEP < 0 should not happen. BSTEP may be 0 or NaN if, e.g., DS or DD becomes Inf.
                if bstep <= 0
                    break
                end
                if ~infnan_obj.is_finite(bstep)
                    info_loc = -1;
                    break
                end

                hd(:) = powalg_obj.hess_mul(d, xpt, pq, 'hq', hq);
                dhd = linalg_obj.inprod(d, hd);

                % Set the step-length ALPHA and update CRVMIN.
                if dhd <= 0
                    alpha = bstep;
                else
                    alpha = min(bstep, gg / dhd);
                    if iter == 1
                        crvmin = dhd / dd;
                    else
                        crvmin = min(crvmin, dhd / dd);
                    end
                end
                % QADD is the reduction of Q due to the new CG step.
                qadd = alpha * (gg - consts_obj.HALF * alpha * dhd);
                % QRED is the reduction of Q up to now.
                qred = qred + qadd;
                % QADD and QRED will be used in the 2-dimensional minimization if any.

                % Update S, HS, and GG.
                sold(:) = s;
                s(:) = s + alpha * d;
                ss = linalg_obj.inprod(s, s);
                hs(:) = hs + alpha * hd;
                ggsav = gg; % Gradient norm square before this iteration
                gg = linalg_obj.inprod(g + hs, g + hs); % Current gradient norm square
                % We may record g+hs for later usage:
                % gnew = g + hs
                % Note that we should NOT set g = g + hs, because g contains the gradient of Q at X.

                % Check whether to exit. This should be done after updating HS and GG, which will be used for
                % the 2-dimensional minimization if any.
                % Exit in case of Inf/NaN in S. This should come the first! Otherwise, we may return an S that
                % contains NaN and fulfills other exit conditions.
                if ~infnan_obj.is_finite(sum(abs(s), 'all'))
                    s(:) = sold;
                    info_loc = -1;
                    break
                end

                % Exit if CG path cuts the boundary. It is the only possibility that TWOD_SEARCH is true.
                if alpha >= bstep || ss >= delsq
                    crvmin = consts_obj.ZERO;
                    twod_search = (n >= 2 && gg > (fortran.power(tol, 2)) * gg0); % TWOD_SEARCH should be FALSE if N = 1.
                    break
                end

                % Exit due to small QADD.
                if qadd <= tol * qred
                    info_loc = 1;
                    break
                end

                % Prepare for the next CG iteration.
                d(:) = (gg / ggsav) * d - g - hs; % CG direction
                dd = linalg_obj.inprod(d, d);
                ds = linalg_obj.inprod(d, s);
                if ds <= 0
                    % DS is positive in theory.
                    info_loc = -1;
                    break
                end
            end

            if ss <= 0 || infnan_obj.is_nan_sp(ss)
                % This may occur for ill-conditioned problems due to rounding.
                info_loc = -1;
                twod_search = false;
            end

            if twod_search
                % At least 1 iteration of 2-dimensional minimization
                maxiter = max(1, maxiter - iter);
            else
                maxiter = 0;
            end

            % The 2-dimensional minimization
            % N.B.: During the iterations, G is NOT updated, and it equals always GOPT, which is the gradient
            % of the trust-region model at the trust-region center X. However, GG is updated: GG = ||G + HS||^2,
            % which is the norm square of the gradient at the current iterate.
            for iter = 1:maxiter
                % Exit if G contains NaN.
                if infnan_obj.is_nan_sp(gg)
                    info_loc = -1;
                    break
                end
                % Exit if GG is small. The inequality must be non-strict so that GG = GG0 = 0 triggers the exit.
                if gg <= (fortran.power(tol, 2)) * gg0
                    info_loc = 0;
                    break
                end
                sg = linalg_obj.inprod(s, g);
                shs = linalg_obj.inprod(s, hs);

                % Begin the 2-dimensional minimization by calculating D and HD and some scalar products.

                % Powell's code calculates D as follows. In precise arithmetic, INPROD(D, S) = 0, ||D|| = ||S||.
                % However, when DELSQ*GG - SGK**2 is tiny, the error in D can be large and hence damage these
                % equalities significantly. This did happen in tests, especially when using the single precision.
                % %sgk = sg + shs
                % %if (sgk / sqrt(gg * delsq) <= tol - ONE) then
                % %    info_loc = 0
                % %    exit
                % %end if
                % %t = sqrt(delsq * gg - sgk**2)
                % %d = (delsq / t) * (g + hs) - (sgk / t) * s

                % We calculate D as below. It did improve the performance of NEWUOA in our test.
                % PROJECT(X, V) returns the projection of X to SPAN(V): X'*(V/||V||)*(V/||V||).
                d(:) = (g + hs) - linalg_obj.project1(g + hs, s);
                % N.B.:
                % 1. The condition ||D||<=SQRT(TOL*GG) below is equivalent to |INPROD(G+HS,S)|<=SQRT((1-TOL)*GG*SS).
                % As given above, Powell's code triggers an exit if INPROD(G+HS,S)=SGK<=(TOL-1)*SQRT(GG*SS).
                % Since |SQRT(1-TOL) - (1-TOL)| <= TOL/2, our condition is close to |SGK| <= (TOL-1)*GG*SS.
                % When |SGK| is tiny, S and G+HS are nearly parallel and hence the 2-dimensional search cannot
                % continue. Note that SGK is unlikely positive if everything goes well.
                % 2. SQRT(TOL)*SQRT(GG) is less likely to encounter underflow than SQRT(TOL*GG).
                % 3. The condition below should be non-strict so that ||D|| = 0 can trigger the exit.
                if linalg_obj.p_norm(d) <= fortran.sqrt(tol) * fortran.sqrt(gg)
                    info_loc = 0;
                    break
                end
                d(:) = (linalg_obj.p_norm(s) / linalg_obj.p_norm(d)) * d;
                % In precise arithmetic, INPROD(D, S) = 0 and ||D|| = ||S|| = DELTA.
                if abs(linalg_obj.inprod(d, s)) >= consts_obj.TENTH * linalg_obj.p_norm(d) * linalg_obj.p_norm(s) || linalg_obj.p_norm(d) >= consts_obj.TWO * delta
                    info_loc = -1;
                    break
                end

                hd(:) = powalg_obj.hess_mul(d, xpt, pq, 'hq', hq);

                % Seek the value of the angle that minimizes Q.
                % First, calculate the coefficients of the objective function on the circle.
                dg = linalg_obj.inprod(d, g);
                dhd = linalg_obj.inprod(hd, d);
                dhs = linalg_obj.inprod(hd, s);
                args(:) = [sg, consts_obj.HALF * (shs - dhd), dg, dhs];
                % The 50 in the line below was chosen by Powell. It works the best in tests, MAGICALLY. Larger
                % (e.g., 60, 100) or smaller (e.g., 20, 40) values will worsen the performance of NEWUOA. Why??
                angle = univar_obj.circle_min(@(varargin) obj.circle_fun_trsapp(varargin{:}), args, 50);

                % Calculate the new S.
                cth = fortran.cos(angle);
                sth = fortran.sin(angle);
                sold(:) = s;
                s(:) = cth * s + sth * d;

                % Exit in case of Inf/NaN in S.
                if ~infnan_obj.is_finite(sum(abs(s), 'all'))
                    s(:) = sold;
                    info_loc = -1;
                    break
                end

                % Test for convergence.
                reduc = obj.circle_fun_trsapp(consts_obj.ZERO, args) - obj.circle_fun_trsapp(angle, args);
                qred = qred + reduc;
                if reduc / qred <= tol
                    info_loc = 1;
                    break
                end

                % Calculate HS.
                hs(:) = cth * hs + sth * hd;
                gg = linalg_obj.inprod(g + hs, g + hs);
            end

            % Set CRVMIN to zero if it is NaN, which may happen if the problem is ill-conditioned.
            if infnan_obj.is_nan_sp(crvmin)
                crvmin = consts_obj.ZERO;
            end

            % Scale CRVMIN back before return. Note that the trust-region step is scale invariant.
            if scaled && crvmin > 0
                crvmin = crvmin / modscal;
            end

            ipObj = inputParser();
            addParameter(ipObj, 'info', NaN);
            parse(ipObj, varargin{:});
            info = ipObj.Results.info;


            %====================%
            %  Calculation ends  %
            %====================%

            % Postconditions
            if consts_obj.DEBUGGING
                debug_obj.assert(numel(s) == n && all(infnan_obj.is_finite(s), 'all'), "SIZE(S) == N, S is finite", srname);
                % Due to rounding, it may happen that ||S|| > DELTA, but ||S|| > 2*DELTA is highly improbable.
                debug_obj.assert(linalg_obj.p_norm(s) <= consts_obj.TWO * delta, "||S|| <= 2*DELTA", srname);
                debug_obj.assert(crvmin >= 0, "CRVMIN >= 0", srname);
            end

        end
        function f = circle_fun_trsapp(~, theta, args)
            %--------------------------------------------------------------------------------------------------%
            % This function defines the objective function of the 2-dimensional search on a circle in TRSAPP.
            %--------------------------------------------------------------------------------------------------%
            consts_obj = consts_mod();
            debug_obj = debug_mod();
            % Inputs



            % Outputs
            f = NaN;

            % Local variables
            srname = "CIRCLE_FUN_TRSAPP";
            cth = NaN;
            sth = NaN;

            % Preconditions
            if consts_obj.DEBUGGING
                debug_obj.assert(numel(args) == 4, "SIZE(ARGS) == 4", srname);
            end

            %====================%
            % Calculation starts %
            %====================%

            cth = fortran.cos(theta);
            sth = fortran.sin(theta);
            f = (args(1) + args(2) * cth) * cth + (args(3) + args(4) * cth) * sth;

            %====================%
            %  Calculation ends  %
            %====================%
        end
        function delta = trrad(~, delta_in, dnorm, eta1, eta2, gamma1, gamma2, ratio)
            %--------------------------------------------------------------------------------------------------%
            % This function updates the trust region radius according to RATIO and DNORM.
            %--------------------------------------------------------------------------------------------------%

            % Generic module
            consts_obj = consts_mod();
            infnan_obj = infnan_mod();
            debug_obj = debug_mod();


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
                %delta = gamma1 * delta_in  ! Powell's COBYLA/LINCOA. Works poorly here.
                %delta = min(gamma1 * delta_in, dnorm)  ! Powell's BOBYQA.

            elseif ratio <= eta2
                delta = max(gamma1 * delta_in, dnorm); % Powell's UOBYQA/NEWUOA/BOBYQA/LINCOA

            else
                delta = max(gamma1 * delta_in, gamma2 * dnorm); % Powell's NEWUOA/BOBYQA.
                %delta = max(delta_in, gamma2 * dnorm)  ! Modified version. Works well for UOBYQA.
                % For noise-free CUTEst problems of <= 200 variables, Powell's version works slightly better
                % than the modified one.
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