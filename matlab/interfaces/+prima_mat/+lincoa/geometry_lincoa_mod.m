classdef geometry_lincoa_mod
    %--------------------------------------------------------------------------------------------------%
    % This module contains subroutines concerning the geometry-improving of the interpolation set XPT.
    %
    % Coded by Zaikun ZHANG (www.zhangzk.net) based on Powell's code.
    %
    % Dedicated to the late Professor M. J. D. Powell FRS (1936--2015).
    %
    % Started: February 2022
    %
    % Last Modified: Sunday, April 21, 2024 PM03:15:36
    %--------------------------------------------------------------------------------------------------%

    methods
        function knew = setdrop_tr(~, idz, kopt, ximproved, bmat, d, delta, rho, xpt, zmat)
            %--------------------------------------------------------------------------------------------------%
            % This subroutine sets KNEW to the index of the interpolation point to be deleted AFTER A TRUST
            % REGION STEP. KNEW will be set in a way ensuring that the geometry of XPT is "optimal" after
            % XPT(:, KNEW) is replaced with XNEW = XOPT + D, where D is the trust-region step.
            % N.B.:
            % It is tempting to take the function value into consideration when defining KNEW, for example,
            % set KNEW so that FVAL(KNEW) = MAX(FVAL) as long as F(XNEW) < MAX(FVAL), unless there is a better
            % choice. However, this is not a good idea, because the definition of KNEW should benefit the
            % quality of the model that interpolates f at XPT. A set of points with low function values is not
            % necessarily a good interpolation set. In contrast, a good interpolation set needs to include
            % points with relatively high function values; otherwise, the interpolant will unlikely reflect the
            % landscape of the function sufficiently.
            %--------------------------------------------------------------------------------------------------%


            powalg_obj = prima_mat.common.powalg_mod();

            % BMAT(N, NPT + N)
            % D(N)


            % XPT(N, NPT)
            % ZMAT(NPT, NPT - N - 1)


            knew = NaN;

            distsq = NaN(size(xpt, 2), 1);

            %====================%
            % Calculation starts %
            %====================%

            % Calculate the distance squares between the interpolation points and the "optimal point". When
            % identifying the optimal point, it is reasonable to take into account the new trust-region trial
            % point XPT(:, KOPT) + D, which will become the optimal point in the next iteration if XIMPROVED
            % is TRUE. Powell suggested this in
            % - (56) of the UOBYQA paper, lines 276--297 of uobyqb.f,
            % - (7.5) and Box 5 of the NEWUOA paper, lines 383--409 of newuob.f,
            % - the last paragraph of page 26 of the BOBYQA paper, lines 435--465 of bobyqb.f.
            % However, Powell's LINCOA code is different. In his code, the KNEW after a trust-region step is
            % picked in lines 72--96 of the update.f for LINCOA, where DISTSQ is calculated as the square of the
            % distance to XPT(KOPT, :) (Powell recorded the interpolation points in rows). However, note that
            % the trust-region trial point has not been included into XPT yet --- it cannot be included without
            % knowing KNEW (see lines 332-344 and 404--431 of lincob.f). Hence Powell's LINCOA code picks KNEW
            % based on the distance to the un-updated "optimal point", which is unreasonable. This has been
            % corrected in our implementation of LINCOA, yet it does not boost the performance.
            if ximproved
                distsq(:) = sum((xpt - (xpt(:, kopt) + d)) .^ 2, 1);
                %%MATLAB: distsq = sum((xpt - (xpt(:, kopt) + d)).^2)  % d should be a column!! Implicit expansion

            else
                distsq(:) = sum((xpt - xpt(:, kopt)) .^ 2, 1);
                %%MATLAB: distsq = sum((xpt - xpt(:, kopt)).^2)  % Implicit expansion
            end
            %distsq = sum((xpt - spread(xpt(:, kopt), dim=2, ncopies=npt))**2, dim=1)  ! Powell's code

            weight = max(1.0, distsq ./ max(0.1 * delta, rho) ^ 2) .^ 3; % Powell's NEWUOA code
            % Other possible definitions of WEIGHT.
            % %weight = distsq**2  ! Powell's code. WRONG.
            % %weight = max(ONE, distsq / max(TENTH * delta, rho)**2)**2.5  ! Worse than power 3
            % %weight = max(ONE, distsq / max(TENTH * delta, rho)**2)**3.5  ! Worse than power 3
            % %weight = (distsq / delta**2)**2   ! Works the same as DISTSQ**2 (as it should be).
            % %weight = (distsq / delta**2)**3  ! Not bad
            % %weight = max(1.0_RP, 10.0_RP * distsq / rho**2)**3
            % %weight = max(1.0_RP, 1.0E2 * distsq / rho**2)**3
            % %weight = max(1.0_RP, 10.0_RP * distsq / delta**2)**3
            % %weight = max(1.0_RP, 1.0E2_RP * distsq / delta**2)**3
            %--------------------------------------------------------------------------------------------------%
            % N.B.: If DISTSQ is the square of distances to the updated XOPT, then it is WRONG to set WEIGHT to
            % DISTSQ**2 or any power of DISTSQ. Why?
            %
            % Consider a scenario where XIMPROVED is TRUE and the new interpolation point XNEW is quite close to
            % one of the existing points in the old XPT, e.g., XPT(:, J). In this case, the only appropriate
            % value of KNEW is J; otherwise, the new interpolation problem will be close to degenerate and its
            % KKT system will be close to singular due to the two close points in the updated interpolation set.
            %
            % What KNEW will be generated if KNEW = MAXLOC(DISTSQ**p * ABS(DEN))? If ||XNEW - XPT(:, J)|| = E,
            % we have the following.
            % 1. DISTSQ(J) = O(E**2) and DEN(J) = O(1);
            % 2. for any K /= J, DIST(K) = O(1) and DEN(K) = O(E).
            % Therefore, for any p > 1/2, KNEW /= J when E is small. As analyzed above, this is inappropriate.
            % In addition, small values of p (e.g., p <= 1/2) always performs poorly for all Powell's methods
            % in our numerical experiments.
            %
            % For the order of DEN, note the DEN(K) is the denominator in the Sherman-Morrison-Woodbury update
            % of the KKT matrix, and DEN(K) = det(new KKT matrix) / det(old KKT matrix), where "KKT matrix"
            % refers to the coefficient matrix of the KKT system for the interpolation problem. See equations
            % (3.10)--(3.12) of the NEWUOA paper for this matrix.
            %
            % Similar arguments can be made if the interpolation is fully determined. Indeed, the usage of
            % DISTSQ as the weight led to a problem in COBYLA during a test on 20230501, which is the very
            % motivation for the current comment.
            %
            % Note that Powell's LINCOA code sets DISTSQ to the square of the distance to the old XOPT, which
            % avoids this problem. However, such a DISTSQ itself seems not ideal, as mentioned above.
            %--------------------------------------------------------------------------------------------------%

            den = powalg_obj.calden(kopt, bmat, d, xpt, zmat, 'idz', idz);
            score = weight .* abs(den);

            % If the new F is not better than FVAL(KOPT), we set SCORE(KOPT) = -1 to avoid KNEW = KOPT.
            if ~ximproved
                score(kopt) = -1.0;
            end

            % SCORE(K) is NaN implies ABS(DEN(K)) is NaN, but we want ABS(DEN) to be big. So we exclude such K.
            score(isnan(score)) = -1.0;

            knew = 0;
            % The following IF works a bit better than `IF (ANY(SCORE > 1) .OR. ANY(SCORE > 0) .AND. XIMPROVED)`
            % from Powell's UOBYQA and NEWUOA code.
            if any(score > 0, 'all')
                % Powell's BOBYQA and LINCOA code
                [~, knew] = max(score);
                %%MATLAB: [~, knew] = max(score);

            end

            % Powell's code does not include the following instructions. With Powell's code, if DEN consists of
            % only NaN, then KNEW can be 0 even when XIMPROVED is TRUE. Here, we set KNEW to the following value,
            % to make sure that the new trial point is included in the interpolation set. However, the updating
            % subroutine will likely need to skip the update of the Lagrange polynomials (i.e., H), or they
            % would be destroyed by the NaNs.
            if (ximproved && knew == 0) || knew < 0
                % KNEW < 0 is impossible in theory.
                [~, knew] = max(distsq);
            end

            %====================%
            %  Calculation ends  %
            %====================%


        end
        function [feasible, s] = geostep(~, iact, idz, knew, kopt, nact, amat, bmat, delbar, qfac, rescon, xpt, zmat, s)
            %--------------------------------------------------------------------------------------------------%
            % This subroutine finds a step S hat intends to improve the geometry of the interpolation set when
            % XPT(:, KNEW) is changed to XOPT + S, where XOPT = XPT(:, KOPT).
            %
            % S is chosen to provide a relatively large value of the modulus of the denominator SIGMA in the
            % updating formula (4.11) of the NEWUOA paper (in theory, SIGMA is positive, yet this may not hold
            % numerically due to rounding errors; in the code, SIGMA is represented by DEN and its modulus by
            % DENABS). This is done by solving
            %
            %   |LFUNC(XOPT + S)|, subject to ||S|| <= DELBAR,
            %
            % because SIGMA >= |LFUNC(XOPT + S)|^2 according to (4.12) of the NEWUOA paper. We do not solve this
            % problem exactly, but calculate three approximate solutions as follows and then choose S from them.
            % 1. The step that maximizes |LFUNC| within the trust region on the lines through XOPT and other
            % interpolation points.
            % 2. The gradient step that maximizes |LFUNC| within the trust region.
            % 3. A projected gradient step that maximizes |LFUNC| within the trust region, the projection being
            % made onto the orthogonal complement of the space spanned by the active gradients.
            %
            % We select S from these three steps by the following criteria.
            % 1. First, set S to either the first or second step, whichever renders a larger value of |SIGMA|.
            % 2. Second, override S by the third step if the latter provides a |SIGMA| that is not small
            % compared with the above one, while leading to a good feasibility.
            %
            % N.B.:
            % 1. The linear constraints are NOT considered in the calculation of the first two steps.
            % 2. For the selection of S, Powell adopted a different set of criteria as follows.
            % 2.1. First, set S to either the first or second step, whichever renders a larger value of |LFUNC|.
            % 2.2. Second, override S by the third step if the latter provides a |LFUNC| that is not small
            % compared with the above one, while being either feasible or with a constraint violation that is at
            % least 0.2*DELBAR.
            % 2.3. If S is not feasible and its constraint violation is less than 0.2*DELBAR, then it is
            % perturbed so that the constraint violation becomes 0.2*DELBAR or more.
            % 3. Powell required the positive constraint violation to be at least 0.2*DELBAR in order to keep
            % the interpolation points apart. In our criteria specified above, we have removed this requirement
            % as we believe that it is implied by the maximization of |LFUNC|. Our criteria works well in tests.
            % 4. In terms of flops, our criteria are more expensive than Powell's due to the evaluation of SIGMA.
            % In derivative-free optimization, we are willing to save function evaluations at the cost of flops.
            % 5. The geometry step of BOBYQA is calculated in a fashion similar to this subroutine: first obtain
            % a step by maximizing |LFUNC| along the lines through XOPT and other interpolation points, then
            % find the Cauchy step, which maximizes |LFUNC| in the 1D space spanned by the gradient of LFUNC,
            % and finally select the geometry step from the aforesaid steps according to the value of |LFUNC| or
            % SIGMA. Yet there still exist three major differences between geometry steps of BOBYQA and LINCOA.
            % 5.1. The geometry step of BOBYQA is calculated subject to the bound constraints, and the computed
            % step is always feasible; the geometry step of LINCOA may violate the linear constraints.
            % 5.2. BOBYQA uses an estimated value of SIGMA (see (3.11) of the BOBYQA paper) to select the step
            % along the lines through XOPT and other interpolation points; LINCOA uses |LFUNC|. According to
            % a test on 20220529, using the estimated SIGMA does not improve the performance of LINCOA.
            % 5.3. LINCOA tries a projected gradient step and prefers this step when its quality is reasonable.
            % BOBYQA does not compute such a step.
            % Additionally, we observe in BOBYQA that it is beneficial for bound constrained problems (but NOT
            % unconstrained ones) to take the Cauchy step only in the late stage of the algorithm, e.g., when
            % DELBAR <= 1.0E-2. Similar phenomenon is not observed in LINCOA, where it worsens the performance
            % of the algorithm to skip the gradient step or the projected gradient step even in the early stage.
            %
            % AMAT, XPT, NACT, IACT, RESCON, QFAC, KOPT are the same as the terms with these names in subroutine
            % LINCOB. KNEW is the index of the interpolation point that is going to be moved. DELBAR is the
            % restriction on the length of S, which is never greater than the current trust region radius DELTA.
            %--------------------------------------------------------------------------------------------------%


            powalg_obj = prima_mat.common.powalg_mod();

            % IACT(M)


            % AMAT(N, M)
            % BMAT(N, NPT+N)

            % QFAC(N, N)
            % RESCON(M)
            % XPT(N, NPT)
            % ZMAT(NPT, NPT-N-1)


            % S(N)


            rstat = NaN(size(amat, 2), 1);

            dderiv = NaN(size(xpt, 2), 1);

            distsq = NaN(size(xpt, 2), 1);
            glag = NaN(size(xpt, 1), 1);

            pglag = NaN(size(xpt, 1), 1);

            pqlag = NaN(size(xpt, 2), 1);

            xopt = NaN(size(xpt, 1), 1);

            n = size(xpt, 1);

            %====================%
            % Calculation starts %
            %====================%

            % Read XOPT.
            xopt(:) = xpt(:, kopt);

            % PQLAG contains the leading NPT elements of the KNEW-th column of H, and it provides the second
            % derivative parameters of LFUNC. Set GLAG to the gradient of LFUNC at the trust region centre.
            pqlag(:) = powalg_obj.omega_col(idz, zmat, knew);
            glag(:) = bmat(:, knew) + powalg_obj.hess_mul(xopt, xpt, pqlag);

            % Maximize |LFUNC| within the trust region on the lines through XOPT and other interpolation points,
            % without considering the linear constraints. In the following, VLAGABS(K) is set to the maximum of
            % |PHI_K(t)| subject to the trust-region constraint with PHI_K(t) = LFUNC((1-t)*XOPT + t*XPT(:, K)).
            dderiv(:) = xpt.' * glag - sum(glag .* xopt, 'all'); % The derivatives PHI_K'(0).
            distsq(:) = sum((xpt - xopt) .^ 2, 1);
            % Set DISTSQ(KOPT) to a positive artificial value. Otherwise, the calculation of STPLEN will raise a
            % floating point exception. This artificial value will NOT be used.
            distsq(kopt) = 1.0;
            % For each K /= KNEW, |PHI_K(t)| is maximized by STPLEN(K), the maximum being VLAGABS(K). Note that
            % PHI_K(t) is a quadratic function with PHI_K'(0) = DDERIV(K) and PHI_K(0) = 0 = PHI_K(1).
            stplen = -delbar ./ sqrt(distsq);
            vlagabs = abs(stplen .* (1.0 - stplen) .* dderiv);
            % The maximization of |PHI_K(t)| is as follows. Note that PHI_K(t) is a quadratic function with
            % PHI_K'(0) = DDERIV(K), PHI_K(0) = 0, and PHI_K(1) = 1.
            if dderiv(knew) * (dderiv(knew) - 1.0) < 0
                stplen(knew) = -stplen(knew);
            end
            vlagabs(knew) = abs(stplen(knew) * dderiv(knew)) + stplen(knew) ^ 2 * abs(dderiv(knew) - 1.0);
            % It does not make sense to consider "the straight line through XOPT and XPT(:, KOPT)". Thus we set
            % VLAGABS(KOPT) to -1 so that KOPT is skipped when we maximize VLAGABS.
            vlagabs(kopt) = -1.0;
            % Find K so that VLAGABS(K) is maximized. We define K in a way slightly different from Powell's
            % code, which sets K to MAXLOC(VLAGABS) by comparing the entries of VLAGABS sequentially.
            % 1. If VLAGABS contains only NaN, which can happen, Powell's code leaves K uninitialized.
            % 2. If VLAGABS(KNEW) = MAXVAL(VLAGABS) = VLAGABS(K) and K < KNEW, Powell's code does not set K=KNEW.
            k = knew;
            if any(vlagabs > vlagabs(knew), 'all')
                [~, k] = max(vlagabs, [], 'omitnan');
                %%MATLAB: [~, k] = max(vlagabs, [], 'omitnan');

            end
            % Set S to the step corresponding to VLAGABS(K), and calculate DENABS for it.
            s(:) = stplen(k) * (xpt(:, k) - xopt);
            den = powalg_obj.calden(kopt, bmat, s, xpt, zmat, 'idz', idz); % Indeed, only DEN(KNEW) is needed.
            denabs = abs(den(knew));

            % Replace S with a steepest ascent step from XOPT if the latter provides a larger value of DENABS.
            gnorm = norm(glag);
            if gnorm > eps(1.0) && isfinite(gnorm)
                gstp = (delbar / gnorm) * glag;
                if sum(gstp .* powalg_obj.hess_mul(gstp, xpt, pqlag), 'all') < 0
                    % <GSTP, HESS_LAG*GSTP> is negative
                    gstp = -gstp;
                end
                den = powalg_obj.calden(kopt, bmat, gstp, xpt, zmat, 'idz', idz); % Indeed, only DEN(KNEW) is needed.
                if abs(den(knew)) > denabs || isnan(denabs)
                    denabs = abs(den(knew));
                    s(:) = gstp;
                end
            end

            % RSTAT identifies the constraints that need evaluation. RSTAT(J) is -1, 0, or 1 respectively means
            % constraint J is irrelevant, active, or inactive and relevant. Do NOT change the order of the lines
            % that set RSTAT, as the later lines override the earlier.
            rstat(:) = 1; % Inactive and relevant
            rstat(abs(rescon) >= delbar) = -1; % Irrelevant
            rstat(iact(1:nact)) = 0; % Active

            % Set FEASIBLE for the calculated S.
            cstrv = max([0.0; amat(:, find(rstat >= 0)).' * s - rescon(find(rstat >= 0))], [], 'all');
            feasible = (cstrv <= 0);

            % If NACT <= 0 or NACT >= N, the calculation has finished. Otherwise, define PGSTP by maximizing
            % |LFUNC| within the trust region from XOPT along the projection of GLAG onto the column space of
            % QFAC(:, NACT+1:N), i.e., the orthogonal complement of the space spanned by the active gradients.
            % In precise arithmetic, moving along PGSTP does not change the values of the active constraints.
            % This projected gradient step is preferred and will override S if it renders a denominator not too
            % small and leads to good feasibility. *** This is critical for the performance of LINCOA. ***
            % In the following, NORMG > EPS prevents floating point exception, and it implies NACT < N.
            pglag(:) = qfac(:, nact + 1:n) * (qfac(:, nact + 1:n).' * glag);
            %%MATLAB: pglag = qfac(:, nact+1:n) * (glag' * qfac(:, nact+1:n))';
            gnorm = norm(pglag);
            if nact > 0 && gnorm > eps(1.0) && isfinite(gnorm)
                pgstp = (delbar / gnorm) * pglag;
                if sum(pgstp .* powalg_obj.hess_mul(pgstp, xpt, pqlag), 'all') < 0
                    % <PGSTP, HESS_LAG*PGSTP> is negative.
                    pgstp = -pgstp;
                end

                % Decide whether to replace S with PGSTP and set FEASIBLE accordingly. CSTRV is the constraint
                % violation of XOPT+PGSTP. Note that we only need to check the constraints that are inactive and
                % relevant, as the value of the active constraints is not changed by moving along PGSTP.
                cstrv = max([0.0; amat(:, find(rstat == 1)).' * pgstp - rescon(find(rstat == 1))], [], 'all');
                % The purpose of CVTOL below is to provide a check on feasibility that includes a tolerance for
                % contributions from computer rounding errors.
                % Powell's code is as follows. Note that MATPROD(PGSTP, AMAT(:, IACT(1:NACT))) is 0 in theory.
                % %cvtol = min(0.01_RP * norm(pgstp), TEN * norm(matprod(pgstp, amat(:, iact(1:nact))), 'inf'))
                % The following code works essentially the same as Powell's code.
                cvtol = max(eps(1.0) * norm(pgstp), 10.0 * norm(amat(:, iact(1:nact)).' * pgstp, "inf"));
                take_pgstp = false;
                if cstrv <= cvtol
                    den = powalg_obj.calden(kopt, bmat, pgstp, xpt, zmat, 'idz', idz); % Indeed, only DEN(KNEW) is needed.
                    take_pgstp = (abs(den(knew)) > 0.1 * denabs);
                end
                if take_pgstp || isnan(denabs)
                    s(:) = pgstp;
                    feasible = (cstrv <= cvtol);
                end
            end

            % In case S is zero or contains Inf/NaN, replace it with a displacement from XPT(:, KNEW) to
            % XOPT. Powell's code does not have this.
            if sum(abs(s), 'all') <= 0 || ~isfinite(sum(abs(s), 'all'))
                s(:) = xpt(:, knew) - xopt;
                scaling = delbar / norm(s);
                s(:) = max(0.6 * scaling, min(0.5, scaling)) * s; % 0.6: ensure |D| > DELBAR/2
                cstrv = max([0.0; amat(:, find(rstat >= 0)).' * s - rescon(find(rstat >= 0))], [], 'all');
                feasible = (cstrv <= 0);
            end

            %====================%
            %  Calculation ends  %
            %====================%


        end

    end
end