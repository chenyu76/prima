classdef trustregion_uobyqa_mod
    %--------------------------------------------------------------------------------------------------%
    % This module provides subroutines concerning the trust-region calculations of UOBYQA.
    %
    % Coded by Zaikun ZHANG (www.zhangzk.net) based on Powell's code and the UOBYQA paper.
    %
    % Dedicated to the late Professor M. J. D. Powell FRS (1936--2015).
    %
    % Started: February 2022
    %
    % Last Modified: Monday, August 07, 2023 AM03:56:41
    %--------------------------------------------------------------------------------------------------%

    methods
        function [d, crvmin] = trstep(~, delta, g, h, tol, d)
            %--------------------------------------------------------------------------------------------------%
            % This subroutine solves the trust-region subproblem
            %
            %     minimize <G, D> + 0.5 * <D, H*D> subject to ||D|| <= DELTA.
            %
            % D will be set to the calculated vector of variables. CRVMIN will be the least eigenvalue of H iff
            % D is a Newton-Raphson step. Then CRVMIN will be positive, but otherwise it will be set to zero.
            % TOL is the value of a tolerance from the open interval (0,1). Let MAXRED be the maximum of
            %   Q(0)-Q(D) subject to ||D|| <= DELTA, and let ACTRED be the value of Q(0)-Q(D) that is actually
            %   calculated. We take the view that any D is acceptable if it has the properties
            %
            %             ||D|| <= DELTA  and  ACTRED <= (1-TOL)*MAXRED.
            %
            % The algorithm first tridiagonalizes H and then applies the More-Sorensen method in
            % More and Sorensen, "Computing a trust region step", SIAM J. Sci. Stat. Comput. 4: 553-572, 1983.
            %
            % The major calculations of the More-Sorensen method lie in the Cholesky factorization or LDL
            % factorization of H + PAR*I with the iteratively selected values of PAR (in More-Sorensen (1983),
            % the parameter is named LAMBDA; in Powell's UOBYQA paper, it is THETA). Powell's method in this
            % code simplifies the calculations by first tridiagonalizing H with an orthogonal transformation.
            % If a matrix T is tridiagonal, its LDL factorization, if exits, can be obtained easily:
            %
            %       T = L*diag(PIV)*L^T,
            %
            % where diag(PIV) is the diagonal matrix with the diagonal entries being PIV(1:N), i.e., "the pivots
            % of the Cholesky factorization" in Powell's comments on his code, and L is the lower triangular
            % matrix with all the diagonal entries being 1, the subdiagonal being the subdiagonal of T divided
            % by PIV(1:N-1), and all the other entries being 0. PIV can be obtained by a simple recursion.
            %
            % For more information, see Section 2 of the UOBYQA paper and
            % Powell, M. J. D., "Trust region calculations revisited", Numerical Analysis 1997: Proceedings of
            % the 17th Dundee Biennial Numerical Analysis Conference, 1997, 193--211,
            % Powell, M. J. D., "The use of band matrices for second derivative approximations in trust region
            % algorithms", Advances in Nonlinear Programming: Proceedings of the 96 International Conference on
            % Nonlinear Programming, 1998, 3--28.
            %--------------------------------------------------------------------------------------------------%



            linalg_obj = prima_mat.common.linalg_mod();


            % G(N)
            % H(N, N)



            % D(N)



            i = NaN;

            k = NaN;


            negcrv = false;


            dhd = NaN;
            dnewton = NaN(numel(g), 1); % Newton-Raphson step; only calculated when N = 1.
            dnorm = NaN;
            dold = NaN(numel(g), 1);

            dtg = NaN;
            dtz = NaN;
            gam = NaN;
            gg = NaN(numel(g), 1);


            hh = NaN(numel(g));


            partmp = NaN;


            phi = NaN;


            piv = NaN(numel(g), 1);
            slope = NaN;
            td = NaN(numel(g), 1);
            tempa = NaN;
            tempb = NaN;
            tn = NaN(numel(g) + -1, 1);
            tnz = NaN;
            wsq = NaN;
            wwsq = NaN;
            z = NaN(numel(g), 1);
            zsq = NaN;


            n = numel(g);


            %====================%
            % Calculation starts %
            %====================%

            % The initial values of DSQ, PHIU, and PHIL are unused but to entertain Fortran compilers.
            % TODO: Check that DSQ, PHIU, PHIL have been initialized before used.
            dsq = 0.0;
            phiu = 0.0;
            phil = 0.0;

            % Scale the problem if G contains large values. Otherwise, floating point exceptions may occur. In
            % the sequel, GG and HH are used instead of G and H, which are INTENT(IN) and hence cannot be
            % changed. Note that CRVMIN must be scaled back if it is nonzero, but the step is scale invariant.
            % N.B.: It is faster and safer to scale by multiplying a reciprocal than by division. See
            % https://fortran-lang.discourse.group/t/ifort-ifort-2021-8-0-1-0e-37-1-0e-38-0/
            if max(abs(g), [], 'all') > 1.0e8
                % The threshold is empirical.
                modscal = max(2.0 * realmin, 1.0 / max(abs(g), [], 'all')); % MAX: precaution against underflow.
                gg(:) = g * modscal;
                hh(:, :) = h * modscal;
                scaled = true;
            else
                modscal = 1.0; % This value is not used, but Fortran compilers may complain without it.
                gg(:) = g;
                hh(:, :) = h;
                scaled = false;
            end

            % Initialize D and CRVMIN.
            d(:) = 0.0;
            crvmin = 0.0;

            gsq = sum(gg .^ 2, 'all');
            gnorm = sqrt(gsq);

            if isnan(gsq)
                return
            end
            if ~any(abs(hh) > 0, 'all')
                if gnorm > 0
                    d(:) = -(delta / gnorm) * gg;
                end
                return
            end

            % Handle the case with N = 1. This should be done after the case where GSQ is NaN.
            % Powell's original code requires that N >= 2.  When N = 1, the code does not work (sometimes even
            % encounters memory errors). This is indeed why the original UOBYQA code constantly terminates with
            % "a trust region step has failed to reduce the quadratic model" when applied to univariate problems.
            if n == 1
                d(:) = delta .* ((-g > 0) .* 2 - 1); %%MATLAB: d = -delta * sign(g)
                if h(1, 1) > 0
                    dnewton(:) = -g ./ h(1, 1);
                    if abs(dnewton(1)) <= delta
                        d(:) = dnewton;
                        crvmin = h(1, 1); % If we use HH(1, 1) here, then we need to scale it back!

                    end
                end
                return
            end

            % Apply Householder transformations to get a tridiagonal matrix similar to H (i.e., the Hessenberg
            % form of H), and put the elements of the Householder vectors in the lower triangular part of HH.
            % Further, TD and TN will contain the diagonal and other nonzero elements of the tridiagonal matrix.
            % In the comments hereafter, H indeed means this tridiagonal matrix.
            [hh, td, tn] = linalg_obj.hessenberg_hhd_trid(hh, td, tn); %%MATLAB: [P, hh] = hess(hh); td = diag(hh); tn = diag(hh, 1)

            % Form GG by applying the similarity transformation.
            for k = 1:n - 1
                gg(k + 1:n) = gg(k + 1:n) - sum(gg(k + 1:n) .* hh(k + 1:n, k), 'all') * hh(k + 1:n, k);
            end
            %%MATLAB: gg = (gg'*P)';  % gg = P'*gg;

            %--------------------------------------------------------------------------------------------------%
            % Zaikun 20220303: Exit if GG, HH, TD, or TN is not finite. Otherwise, the behavior of this
            % subroutine is not predictable. For example, if HNORM = GNORM = Inf, it is observed that the
            % initial value of PARL defined below will change when we add code that should not affect PARL
            % (e.g., print it, or add TD = 0, TN = 0, PIV = 0 at the beginning of this subroutine).
            % This is probably because the behavior of MAX is undefined if it receives NaN (if GNORM and HNORM
            % are both Inf, then GNORM/DELTA - HNORM = NaN).
            %--------------------------------------------------------------------------------------------------%
            if ~isfinite(sum(abs(gg), 'all') + sum(abs(hh), 'all') + sum(abs(td), 'all') + sum(abs(tn), 'all'))
                return
            end

            % Begin the trust region calculation with a tridiagonal matrix by calculating the L_1-norm of the
            % Hessenberg form of H, which is an upper bound for the spectral norm of H.
            hnorm = max(abs([0.0; tn]) + abs(td) + abs([tn; 0.0]), [], 'all');
            delsq = delta * delta;

            % Set the initial values of PAR and its bounds.
            % N.B.: PAR is the parameter LAMBDA in More-Sorensen (1983) and Powell (1997), as well as the THETA
            % in Section 2 of the UOBYQA paper. The algorithm looks for the optimal PAR characterized in Lemmas
            % 2.1--2.3 of More-Sorensen (1983).
            parl = max([0.0, -min(td, [], 'all'), gnorm / delta - hnorm], [], 'all'); % Lower bound for the optimal PAR
            parlest = parl; % Estimation for PARL
            par = parl;
            paru = 0.0; % Upper bound for the optimal PAR ??? The initial value is less than PARL. Why?
            paruest = 0.0; % Estimation for PARU
            posdef = false;
            dold(:) = 0.0;

            maxiter = min(1000, 100 * n); % Unlikely to be reached.
            % Zaikun 26-06-2019: Powell's original code can encounter infinite cycling, which did happen when
            % testing the CUTEst problems GAUSS1LS, GAUSS2LS, and GAUSS3LS. Indeed, in all these cases, Inf
            % and NaN appear in D due to extremely large values in the Hessian matrix (up to 10^219).

            for iter = 1:maxiter
                if isfinite(sum(abs(d), 'all'))
                    dold(:) = d;
                else
                    d(:) = dold;
                    break
                end
                if iter > maxiter
                    break
                end

                % Calculate the pivots of the Cholesky factorization of (H + PAR*I), which correspond to the
                % squares of the diagonal entries of L in the Cholesky factorization LL^T, or the diagonal
                % matrix in the LDL factorization. After getting PIV, we can get the LDL factorization of
                % H + PAR*I easily: it is L*diag(PIV)*L^T, where diag(PIV) is the diagonal matrix with PIV being
                % the diagonal, and L is the lower triangular matrix with all the diagonal entries being 1, the
                % subdiagonal being the vector TN/PIV(1:N-1) (entrywise), and all the other entries being 0.
                piv(:) = 0.0; % Initialize PIV, so that we know that any NaN in PIV is due to the loop below.
                piv(1) = td(1) + par;
                % Powell implemented the loop by a GOTO, and K = N when the loop exits. It may not be true here.
                for k = 1:n - 1
                    if piv(k) > 0
                        piv(k + 1) = td(k + 1) + par - tn(k) ^ 2 / piv(k);
                    elseif abs(piv(k)) + abs(tn(k)) <= 0
                        % PIV(K) == 0 == TN(K)
                        piv(k + 1) = td(k + 1) + par;
                    else                        % PIV(K) < 0 .OR. (PIV(K) == 0 .AND. TN(K) /= 0)
                        break
                    end
                end

                % Zaikun 20220509
                if any(isnan(piv), 'all')
                    break % Better action to take???

                end

                % NEGCRV is TRUE iff H + PAR*I has at least one negative eigenvalue (CRV means curvature).
                negcrv = any(piv < 0 | (piv <= 0 & abs([tn; 0.0]) > 0), 'all');

                % Handle the case where H + PAR*I is positive semidefinite and the gradient at the trust region
                % center is zero.
                if gsq <= 0 && ~negcrv
                    paru = par;
                    paruest = par;
                    if par <= 0
                        % PAR == 0. A rare case: the trust region center is optimal.
                        break
                    end
                end

                if negcrv
                    % Set K to the first index corresponding to a negative curvature.
                    % N.B.: In theory, we need not prepend N to TRUELOC(...), because TRUELOC must return
                    % a nonempty array when NEGCRV is TRUE, and hence K <= N; however, the Fortran code may not
                    % behave in this way when compiled with aggressive optimization options; on 20221220, it is
                    % observed that K = HUGE(K) = 32767 with Flang -Ofast.
                    k = min([n; find(piv < 0 | (piv <= 0 & abs([tn; 0.0]) > 0))], [], 'all');
                else
                    % Set K to the last index corresponding to a zero curvature; K = 0 if no such curvature exits.
                    k = max([0; find(abs(piv) + abs([tn; 0.0]) <= 0)], [], 'all');
                end

                % At this point, K == 0 iff H + PAR*I is positive definite.
                % Handle the case where H + PAR*I has at least one nonpositive eigenvalue.
                if k >= 1
                    % Set D to a direction of nonpositive curvature of the tridiagonal matrix, and revise PARLEST.

                    %------------------------------------------------------------------------------------------%
                    % Zaikun 20220512: Powell's code does not include the following initialization. Consequently,
                    % D(KSAV+1:N) or D(KSAV+2:N) will not be initialized but inherit values from the previous
                    % iteration. Is this intended?
                    d(:) = 0.0;
                    %------------------------------------------------------------------------------------------%

                    d(k) = 1.0; % Zaikun 20220512: D(K+1:N) = ?

                    %------------------------------------------------------------------------------------------%
                    % The code until "Terminate with D set to a multiple of the current D ..." sets only D(1:KSAV)
                    % or D(1:KSAV+1), with the KSAV defined later. D_INITIALIZED indicates whether D(1:N) is
                    % fully initialized in this process (TRUE) or not (FALSE). See the comments above for details.
                    %------------------------------------------------------------------------------------------%

                    dhd = piv(k);

                    % In Fortran, the following two IFs CANNOT be merged into
                    % IF(K < N .AND. ABS(TN(K)) > ABS(PIV(K))).
                    % This is because Fortran may not perform a short-circuit evaluation of this logic expression,
                    % and hence TN(K) may be accessed even if K >= N, leading to an out-of-boundary index since
                    % SIZE(TN) is only N-1. This is not a problem in C, MATLAB, Python, Julia, or R, where short
                    % circuit is ensured.
                    if k < n
                        if abs(tn(k)) > abs(piv(k))
                            % PIV(K+1) was named as "TEMP" in Powell's code. Is PIV(K+1) consistent with the meaning of PIV?
                            piv(k + 1) = td(k + 1) + par;
                            if piv(k + 1) <= abs(piv(k))
                                d(k + 1) = 1.0 .* ((-tn(k) > 0) .* 2 - 1); %%MATLAB: d(k + 1) = -sing(tn(k))
                                dhd = piv(k) + piv(k + 1) - 2.0 * abs(tn(k));
                            else
                                d(k + 1) = -tn(k) / piv(k + 1);
                                dhd = piv(k) + tn(k) * d(k + 1);
                            end
                        end
                    end

                    for i = k - 1:-1:1
                        % It may happen that TN(I) == 0 == PIV(I). Without checking TN(I), we will get D(I)=NaN.
                        % Once we encounter a zero TN(I), D(I) is set to zero, and D(1:I-1) will consequently be
                        % zero as well, because D(J) is a multiple of D(J+1) for each J.
                        if abs(tn(i)) > 0
                            d(i) = -tn(i) * d(i + 1) / piv(i);
                        else
                            d(1:i) = 0.0;
                            break
                        end
                    end

                    dsq = sum(d .^ 2, 'all');
                    parl = par;
                    parlest = par - dhd / dsq;
                end

                if gsq <= 0 || k >= 1
                    % Handle the case where the gradient at the trust region center is zero or H + PAR*I is not
                    % positive definite.

                    % Terminate with D set to a multiple of the current D if the following test suggests so.
                    if gsq <= 0
                        partmp = paruest * (1.0 - tol);
                    else
                        partmp = paruest;
                    end
                    if paruest > 0 && parlest >= partmp
                        %--------------------------------------------------------------------------------------%
                        % Zaikun 20220512:
                        % The definition of D below requires that D is initialized. In Powell's code, it may
                        % happen that only D(1:KSAV) or D(1:KSAV+1) is initialized during the current iteration,
                        % but the other entries are inherited from the previous iteration OR from the initial
                        % value before the iterations start, which is 0. If such inheriting happens,
                        % D_INITIALIZED will be FALSE. In tests on 20220514, both cases did occur.
                        % Interestingly, in both cases, the inherited values were all zero or close to zero
                        % (1E-16), and hence not very different from the initial value zero that we set above.
                        % Is this intended?
                        %--------------------------------------------------------------------------------------%

                        dtg = sum(d .* gg, 'all');
                        if dtg > 0
                            % Has DSQ got the correct value?
                            d(:) = -(delta / sqrt(dsq)) * d;
                        else                            % This ELSE covers the unlikely yet possible case where DTG is zero or even NaN.
                            d(:) = (delta / sqrt(dsq)) * d;
                        end
                        % N.B.: As per Powell's code, the lines above would be D = -SIGN(DELTA/SQRT(DSQ), DTG)*D.
                        % However, our version here seems more reasonable in case DTG == 0, which is unlikely
                        % but did happen numerically. Note that SIGN(A, 0) = |A| /= -SIGN(A, 0).
                        break
                    end
                else

                    % Handle the case where the gradient at the trust region center is nonzero and H + PAR*I
                    % is positive definite.
                    % Calculate D = -(H + PAR*I)^{-1}*G for the current PAR. The loops below find D using the
                    % LDL factorization of the (tridiagonalized) H + PAR*I = L*diag(PIV)*L^T.
                    d(1) = -gg(1) / piv(1);
                    % The loop sets D = -PIV^{-1}L^{-1}*GG
                    for k = 1:n - 1
                        d(k + 1) = -(gg(k + 1) + tn(k) * d(k)) / piv(k + 1);
                    end
                    wsq = sum(piv .* d .^ 2, 'all'); % GG^T*(H+PAR*I)^{-1}*GG. Needed in the convergence test.
                    % The loop sets D = L^{-T}*D = -L^{-T}*PIV^{-1}*L^{-1}*GG = -(H+PAR*I)^{-1}*GG.
                    for k = n - 1:-1:1
                        d(k) = d(k) - tn(k) * d(k + 1) / piv(k);
                    end

                    if ~isfinite(sum(abs(d), 'all'))
                        d(:) = dold;
                        break
                    end

                    dsq = sum(d .^ 2, 'all');

                    % Return if the Newton-Raphson step is feasible, setting CRVMIN to the least eigenvalue of H.
                    if par <= 0 && dsq <= delsq
                        % PAR <= 0 indeed means PAR == 0.
                        crvmin = linalg_obj.eigmin_sym_trid(td, tn, 'tol', 1.0e-2);
                        %%MATLAB:
                        %%% It is critical for the efficiency to use `spdiags` to construct `tridh` sparsely.
                        %%tridh = spdiags([[tn; 0], td, [0; tn]], -1:1, n, n);
                        %%crvmin = eigs(tridh, 1, 'smallestreal');
                        break
                    end

                    % Make the usual test for acceptability of a full trust region step.
                    dnorm = sqrt(dsq);

                    phi = 1.0 / dnorm - 1.0 / delta;
                    if tol * (1.0 + par * dsq / wsq) - dsq * phi * phi >= 0
                        d(:) = (delta / dnorm) * d;
                        break
                    end
                    if iter >= 2 && par <= parl
                        break
                    end
                    if paru > 0 && par >= paru
                        break
                    end

                    % Complete the iteration when PHI is negative.
                    if phi < 0
                        parlest = par;
                        if posdef
                            if phi <= phil
                                break % Has PHIL got the correct value

                            end
                            slope = (phi - phil) / (par - parl);
                            parlest = par - phi / slope;
                        end
                        if paru > 0
                            slope = (phiu - phi) / (paru - par); % Has PHIU got the correct value?

                        else
                            slope = 1.0 / gnorm;
                        end
                        partmp = par - phi / slope;
                        if paruest > 0
                            paruest = min(partmp, paruest);
                        else
                            paruest = partmp;
                        end
                        posdef = true;
                        parl = par;
                        phil = phi;
                    else

                        % If required, calculate Z for the alternative test for convergence.
                        % For Z, see the discussions below (16) in Section 2 of the UOBYQA paper (the 2002 version
                        % in Math. Program.; in the DAMTP 2000/NA14 report, it is below (2.8) in Section 2). The two
                        % loops below find Z using the LDL factorization of the (tridiagonalized) H + PAR*I.
                        if ~posdef
                            z(1) = 1.0 / piv(1);
                            for k = 1:n - 1
                                tnz = tn(k) * z(k);
                                if tnz > 0
                                    z(k + 1) = -(1.0 + tnz) / piv(k + 1);
                                else
                                    z(k + 1) = (1.0 - tnz) / piv(k + 1);
                                end
                            end
                            wwsq = sum(piv .* z .^ 2, 'all'); % Needed in the convergence test.
                            for k = n - 1:-1:1
                                z(k) = z(k) - tn(k) * z(k + 1) / piv(k);
                            end

                            zsq = sum(z .^ 2, 'all');
                            dtz = sum(d .* z, 'all');

                            % Apply the alternative test for convergence.
                            tempa = abs(delsq - dsq);
                            tempb = sqrt(dtz * dtz + tempa * zsq);
                            if abs(dtz) > 0
                                gam = tempa / (tempb .* ((dtz > 0) .* 2 - 1) + dtz); %%MATLAB: gam = tempa / (sign(dtz)*tempb + dtz)

                            else                                % This ELSE covers the unlikely yet possible case where DTZ is zero or even NaN.
                                gam = sqrt(tempa / zsq);
                            end
                            if tol * (wsq + par * delsq) - gam * gam * wwsq >= 0
                                d(:) = d + gam * z;
                                break
                            end
                            parlest = max(parlest, par - wwsq / zsq);
                        end

                        % Complete the iteration when PHI is positive.
                        slope = 1.0 / gnorm;
                        if paru > 0
                            if phi >= phiu
                                break % Has PHIU got the correct value?

                            end
                            slope = (phiu - phi) / (paru - par);
                        end
                        parlest = max(parlest, par - phi / slope);
                        paruest = par;
                        if posdef
                            slope = (phi - phil) / (par - parl); % Has PHIL got the correct value?
                            paruest = par - phi / slope;
                        end
                        paru = par;
                        phiu = phi;
                    end
                end

                % Pick the value of PAR for the next iteration.
                if paru <= 0
                    % PARU == 0
                    par = 2.0 * parlest + gnorm / delta;
                else
                    par = 0.5 * (parl + paru);
                    par = max(par, parlest);
                end
                if paruest > 0
                    par = min(par, paruest);
                end
            end

            % Apply the inverse Householder transformations to recover D.
            for k = n - 1:-1:1
                d(k + 1:n) = d(k + 1:n) - sum(d(k + 1:n) .* hh(k + 1:n, k), 'all') * hh(k + 1:n, k);
            end
            %%MATLAB: d = P*d;

            % If the More-Sorensen algorithm breaks down abnormally (e.g., NaN in the computation), then ||D||
            % may be (much) more than DELTA. This is handled in the following naive way.
            if norm(d) > delta
                d(:) = (delta / norm(d)) * d;
            end

            % Set CRVMIN to zero if it is NaN, which may happen if the problem is ill-conditioned.
            if isnan(crvmin)
                crvmin = 0.0;
            end

            % Scale CRVMIN back before return. Note that the trust-region step is scale invariant.
            if scaled && crvmin > 0
                crvmin = crvmin / modscal;
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
                delta = gamma1 * dnorm; % Powell's UOBYQA/NEWUOA
                %delta = gamma1 * delta_in  ! Powell's COBYLA/LINCOA.
                %delta = min(gamma1 * delta_in, dnorm)  ! Powell's BOBYQA.

            elseif ratio <= eta2
                delta = max(gamma1 * delta_in, dnorm); % Powell's UOBYQA/NEWUOA/BOBYQA/LINCOA

            else
                delta = max(gamma1 * delta_in, gamma2 * dnorm); % Powell's NEWUOA/BOBYQA. Works well for UOBYQA.
                %delta = max(delta_in, 1.25_RP * dnorm, dnorm + rho)  ! Powell's original UOBYQA code.
                %delta = max(delta_in, gamma2 * dnorm)  ! This works evidently better than Powell's version.
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