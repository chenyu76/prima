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



            linalg_obj = prima_mat.common.linalg_mod();
            powalg_obj = prima_mat.common.powalg_mod();
            univar_obj = prima_mat.common.univar_mod();


            % GOPT_IN(N)
            % HQ_IN(N, N)
            % PQ_IN(NPT)

            % XPT(N, NPT)



            % S(N)



            info_loc = NaN;
            iter = NaN;


            alpha = NaN;
            angle = NaN;
            args = NaN(4, 1);
            bstep = NaN;
            cth = NaN;


            dg = NaN;
            dhd = NaN;
            dhs = NaN;


            ggsav = NaN;
            gopt = NaN(numel(gopt_in), 1);
            hd = NaN(numel(gopt_in), 1);
            hq = NaN(size(hq_in, 1), size(hq_in, 2));
            hs = NaN(numel(gopt_in), 1);

            pq = NaN(numel(pq_in), 1);
            qadd = NaN;

            reduc = NaN;
            resid = NaN;
            sg = NaN;
            shs = NaN;
            sold = NaN(numel(gopt_in), 1);
            sqrtd = NaN;

            sth = NaN;


            n = size(xpt, 1);


            %====================%
            % Calculation starts %
            %====================%

            % Scale the problem if GOPT contains large values. Otherwise, floating point exceptions may occur.
            % Note that CRVMIN must be scaled back if it is nonzero, but the step is scale invariant.
            % N.B.: It is faster and safer to scale by multiplying a reciprocal than by division. See
            % https://fortran-lang.discourse.group/t/ifort-ifort-2021-8-0-1-0e-37-1-0e-38-0/
            if max(abs(gopt_in), [], 'all') > 1.0e12
                % The threshold is empirical.
                modscal = max(2.0 * realmin, 1.0 / max(abs(gopt_in), [], 'all')); % MAX: precaution against underflow.
                gopt(:) = gopt_in * modscal;
                pq(:) = pq_in * modscal;
                hq(:, :) = hq_in * modscal;
                scaled = true;
            else
                modscal = 1.0; % This value is not used, but Fortran compilers may complain without it.
                gopt(:) = gopt_in;
                pq(:) = pq_in;
                hq(:, :) = hq_in;
                scaled = false;
            end

            s(:) = 0.0;
            crvmin = 0.0;
            qred = 0.0;
            info_loc = 2; % Default exit flag is 2, i.e., MAXITER is attained

            % Prepare for the first line search.
            %--------------------------------------------------------------------------------------------------%
            % N.B.: During the iterations, G is NOT updated, and it equals always GOPT, which is the gradient
            % of the trust-region model at the trust-region center X. However, GG is updated: GG = ||G + HS||^2,
            % which is the norm square of the gradient at the current iterate.
            g = gopt;
            gg = sum(g .* g, 'all');
            %--------------------------------------------------------------------------------------------------%
            gg0 = gg;
            d = -g;
            dd = gg;
            ds = 0.0;
            ss = 0.0;
            hs(:) = 0.0;
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
                if isnan(gg)
                    info_loc = -1;
                    break
                end
                % Exit if GG is small. This must be done first; otherwise, DD can be 0 and BSTEP is not well
                % defined. The inequality below must be non-strict so that GG = GG0 = 0 will trigger the exit.
                if gg <= (tol ^ 2) * gg0
                    info_loc = 0;
                    break
                end

                % Set BSTEP to the step length such that ||S + BSTEP*D|| = DELTA.
                if iter == 1
                    bstep = delta / sqrt(dd);
                else
                    resid = delsq - ss;
                    if resid <= 0
                        twod_search = true;
                        break
                    end

                    % Powell's code does not have the following two IFs.
                    %--------------------------------------------------%
                    if dd <= eps(1.0) * delsq
                        info_loc = 0;
                        break
                    end
                    if isnan(ds)
                        info_loc = -1;
                        break
                    end
                    %--------------------------------------------------%

                    % SQRTD: square root of a discriminant. The MAXVAL avoids SQRTD < ABS(DS) due to underflow.
                    sqrtd = max([sqrt(ds ^ 2 + dd * resid), abs(ds), sqrt(dd * resid)], [], 'all');
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
                if ~isfinite(bstep)
                    info_loc = -1;
                    break
                end

                hd(:) = powalg_obj.hess_mul(d, xpt, pq, 'hq', hq);
                dhd = sum(d .* hd, 'all');

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
                qadd = alpha * (gg - 0.5 * alpha * dhd);
                % QRED is the reduction of Q up to now.
                qred = qred + qadd;
                % QADD and QRED will be used in the 2-dimensional minimization if any.

                % Update S, HS, and GG.
                sold(:) = s;
                s(:) = s + alpha * d;
                ss = sum(s .* s, 'all');
                hs = hs + alpha * hd;
                ggsav = gg; % Gradient norm square before this iteration
                gg = sum((g + hs) .* (g + hs), 'all'); % Current gradient norm square
                % We may record g+hs for later usage:
                % gnew = g + hs
                % Note that we should NOT set g = g + hs, because g contains the gradient of Q at X.

                % Check whether to exit. This should be done after updating HS and GG, which will be used for
                % the 2-dimensional minimization if any.
                % Exit in case of Inf/NaN in S. This should come the first! Otherwise, we may return an S that
                % contains NaN and fulfills other exit conditions.
                if ~isfinite(sum(abs(s), 'all'))
                    s(:) = sold;
                    info_loc = -1;
                    break
                end

                % Exit if CG path cuts the boundary. It is the only possibility that TWOD_SEARCH is true.
                if alpha >= bstep || ss >= delsq
                    crvmin = 0.0;
                    twod_search = (n >= 2 && gg > (tol ^ 2) * gg0); % TWOD_SEARCH should be FALSE if N = 1.
                    break
                end

                % Exit due to small QADD.
                if qadd <= tol * qred
                    info_loc = 1;
                    break
                end

                % Prepare for the next CG iteration.
                d = (gg / ggsav) * d - g - hs; % CG direction
                dd = sum(d .* d, 'all');
                ds = sum(d .* s, 'all');
                if ds <= 0
                    % DS is positive in theory.
                    info_loc = -1;
                    break
                end
            end

            if ss <= 0 || isnan(ss)
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
                if isnan(gg)
                    info_loc = -1;
                    break
                end
                % Exit if GG is small. The inequality must be non-strict so that GG = GG0 = 0 triggers the exit.
                if gg <= (tol ^ 2) * gg0
                    info_loc = 0;
                    break
                end
                sg = sum(s .* g, 'all');
                shs = sum(s .* hs, 'all');

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
                d = (g + hs) - linalg_obj.project1(g + hs, s);
                % N.B.:
                % 1. The condition ||D||<=SQRT(TOL*GG) below is equivalent to |INPROD(G+HS,S)|<=SQRT((1-TOL)*GG*SS).
                % As given above, Powell's code triggers an exit if INPROD(G+HS,S)=SGK<=(TOL-1)*SQRT(GG*SS).
                % Since |SQRT(1-TOL) - (1-TOL)| <= TOL/2, our condition is close to |SGK| <= (TOL-1)*GG*SS.
                % When |SGK| is tiny, S and G+HS are nearly parallel and hence the 2-dimensional search cannot
                % continue. Note that SGK is unlikely positive if everything goes well.
                % 2. SQRT(TOL)*SQRT(GG) is less likely to encounter underflow than SQRT(TOL*GG).
                % 3. The condition below should be non-strict so that ||D|| = 0 can trigger the exit.
                if norm(d) <= sqrt(tol) * sqrt(gg)
                    info_loc = 0;
                    break
                end
                d = (norm(s) / norm(d)) * d;
                % In precise arithmetic, INPROD(D, S) = 0 and ||D|| = ||S|| = DELTA.
                if abs(sum(d .* s, 'all')) >= 0.1 * norm(d) * norm(s) || norm(d) >= 2.0 * delta
                    info_loc = -1;
                    break
                end

                hd(:) = powalg_obj.hess_mul(d, xpt, pq, 'hq', hq);

                % Seek the value of the angle that minimizes Q.
                % First, calculate the coefficients of the objective function on the circle.
                dg = sum(d .* g, 'all');
                dhd = sum(hd .* d, 'all');
                dhs = sum(hd .* s, 'all');
                args(:) = [sg, 0.5 * (shs - dhd), dg, dhs];
                % The 50 in the line below was chosen by Powell. It works the best in tests, MAGICALLY. Larger
                % (e.g., 60, 100) or smaller (e.g., 20, 40) values will worsen the performance of NEWUOA. Why??
                angle = univar_obj.circle_min(@(varargin) obj.circle_fun_trsapp(varargin{:}), args, 50);

                % Calculate the new S.
                cth = cos(angle);
                sth = sin(angle);
                sold(:) = s;
                s(:) = cth * s + sth * d;

                % Exit in case of Inf/NaN in S.
                if ~isfinite(sum(abs(s), 'all'))
                    s(:) = sold;
                    info_loc = -1;
                    break
                end

                % Test for convergence.
                reduc = obj.circle_fun_trsapp(0.0, args) - obj.circle_fun_trsapp(angle, args);
                qred = qred + reduc;
                if reduc / qred <= tol
                    info_loc = 1;
                    break
                end

                % Calculate HS.
                hs = cth * hs + sth * hd;
                gg = sum((g + hs) .* (g + hs), 'all');
            end

            % Set CRVMIN to zero if it is NaN, which may happen if the problem is ill-conditioned.
            if isnan(crvmin)
                crvmin = 0.0;
            end

            % Scale CRVMIN back before return. Note that the trust-region step is scale invariant.
            if scaled && crvmin > 0
                crvmin = crvmin / modscal;
            end

            ipObj = inputParser();
            addParameter(ipObj, 'info', NaN);
            parse(ipObj, varargin{:});


            %====================%
            %  Calculation ends  %
            %====================%



        end
        function f = circle_fun_trsapp(~, theta, args)
            %--------------------------------------------------------------------------------------------------%
            % This function defines the objective function of the 2-dimensional search on a circle in TRSAPP.
            %--------------------------------------------------------------------------------------------------%



            f = NaN;


            %====================%
            % Calculation starts %
            %====================%

            cth = cos(theta);
            sth = sin(theta);
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



            % Input
            % Current trust-region radius
            % Norm of current trust-region step
            % Ratio threshold for contraction
            % Ratio threshold for expansion
            % Contraction factor
            % Expansion factor
            % Reduction ratio


            delta = NaN;


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



        end

    end
end