classdef geometry_bobyqa_mod
    %--------------------------------------------------------------------------------------------------%
    % This module contains subroutines concerning the geometry-improving of the interpolation set XPT.
    %
    % Coded by Zaikun ZHANG (www.zhangzk.net) based on Powell's code and the BOBYQA paper.
    %
    % Dedicated to the late Professor M. J. D. Powell FRS (1936--2015).
    %
    % Started: February 2022
    %
    % Last Modified: Tue 10 Feb 2026 02:07:34 PM CET
    %--------------------------------------------------------------------------------------------------%

    methods
        function knew = setdrop_tr(~, kopt, ximproved, bmat, d, delta, rho, xpt, zmat)
            %--------------------------------------------------------------------------------------------------%
            % This subroutine sets KNEW to the index of the interpolation point to be deleted AFTER A TRUST
            % REGION STEP. KNEW will be set in a way ensuring that the geometry of XPT is "optimal" after
            % XPT(:, KNEW) is replaced with XNEW = XOPT + D, where D is the trust-region step. See discussions
            % around (6.1) of the BOBYQA paper.
            % N.B.:
            % 1. If XIMPROVED = TRUE, then KNEW > 0 so that XNEW is included into XPT. Otherwise, it is a bug.
            % 2. If XIMPROVED = FALSE, then KNEW /= KOPT so that XPT(:, KOPT) stays. Otherwise, it is a bug.
            % 3. It is tempting to take the function value into consideration when defining KNEW, for example,
            % set KNEW so that FVAL(KNEW) = MAX(FVAL) as long as F(XNEW) < MAX(FVAL), unless there is a better
            % choice. However, this is not a good idea, because the definition of KNEW should benefit the
            % quality of the model that interpolates f at XPT. A set of points with low function values is not
            % necessarily a good interpolation set. In contrast, a good interpolation set needs to include
            % points with relatively high function values; otherwise, the interpolant will unlikely reflect the
            % landscape of the function sufficiently.
            %--------------------------------------------------------------------------------------------------%

            % Common modules



            powalg_obj = prima_mat.common.powalg_mod();


            % Inputs


            % BMAT(N, NPT + N)
            % D(N)


            % XPT(N, NPT)
            % ZMAT(NPT, NPT - N - 1)

            % Outputs
            knew = NaN;

            % Local variables



            distsq = NaN(size(xpt, 2), 1);


            % Sizes



            % Preconditions


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
                %%MATLAB: distsq = sum((xpt - (xpt(:, kopt) + d)).^2)  % d should be a column! Implicit expansion

            else
                distsq(:) = sum((xpt - xpt(:, kopt)) .^ 2, 1);
                %%MATLAB: distsq = sum((xpt - xpt(:, kopt)).^2)  % Implicit expansion
            end

            weight = max(1.0, distsq ./ rho ^ 2) .^ 4;
            % Other possible definitions of WEIGHT.
            % %weight = max(ONE, distsq / rho**2)**3.5  ! Quite similar to power 4
            % %weight = max(ONE, distsq / rho**2)**3  ! Not bad
            % %weight = max(ONE, distsq / delta**2)**2  ! Powell's code. Does not works as well as the above.
            % %weight = max(ONE, distsq / rho**2)**2  ! Similar to Powell's code, not better.
            % %weight = max(ONE, distsq / delta**2)  ! Defined in (6.1) of the BOBYQA paper. It works poorly!
            % %weight = max(ONE, distsq / max(TENTH * delta, rho)**2)**3.5  ! The same as DISTSQ/RHO**2.
            % The following WEIGHT all perform a bit worse than the above one.
            % %weight = max(ONE, distsq / delta**2)**3.5
            % %weight = max(ONE, distsq / delta**2)**2.5
            % %weight = max(ONE, distsq / delta**2)**3
            % %weight = max(ONE, distsq / delta**2)**4
            % %weight = max(ONE, distsq / delta**2)**4.5
            % %weight = max(ONE, distsq / rho**2)**2.5
            % %weight = max(ONE, distsq / rho**2)**3
            % %weight = max(ONE, distsq / rho**2)**4.5

            % Different from NEWUOA/LINCOA, the possibility that entries in DEN become negative is handled by
            % RESCUE. Hence the SCORE here uses DEN in contrast to ABS(DEN) in NEWUOA/LINCOA.
            den = powalg_obj.calden(kopt, bmat, d, xpt, zmat);
            score = weight .* den;

            % If the new F is not better than FVAL(KOPT), we set SCORE(KOPT) = -1 to avoid KNEW = KOPT.
            if ~ximproved
                score(kopt) = -1.0;
            end

            % SCORE(K) = NaN implies DEN(K) = NaN. We exclude such K as we want DEN to be big.
            score(isnan(score)) = -1.0;

            knew = 0;
            % The following IF works slightly better than `IF (ANY(SCORE > 0))` from Powell's BOBYQA/LINCOA code.
            if any(score > 1, 'all') || (ximproved && any(score > 0, 'all'))
                % Powell's UOBYQA and NEWUOA code.
                % See (6.1) of the BOBYQA paper for the definition of KNEW in this case.
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

            % Postconditions


        end
        function d = geostep(~, knew, kopt, bmat, delbar, sl, su, xpt, zmat)
            %--------------------------------------------------------------------------------------------------%
            % This subroutine finds a step D that intends to improve the geometry of the interpolation set
            % when XPT(:, KNEW) is changed to XOPT + D, where XOPT = XPT(:, KOPT). See Section 3 of the BOBYQA
            % paper, particularly the discussions starting from (3.7).
            %
            % The arguments XPT, BMAT, ZMAT, SL and SU all have the same meanings as in BOBYQB.
            % KOPT is the index of the optimal interpolation point.
            % KNEW is the index of the interpolation point that is going to be moved.
            % DELBAR is the trust region bound for the geometry step.
            % XLINE will be a suitable new position for the interpolation point XPT(:, KNEW). Specifically, it
            %   satisfies the SL, SU and trust region bounds and it should provide a large denominator in the
            %   next call of UPDATE. The step XLINE-XOPT from XOPT is restricted to moves along the straight
            %   lines through XOPT and another interpolation point.
            % XCAUCHY provides a large value of the modulus of the KNEW-th Lagrange function subject to the
            %   constraints that have been mentioned, its main difference from XLINE being that XCAUCHY-XOPT
            %   is a bound-constrained version of the Cauchy step within the trust region.
            %--------------------------------------------------------------------------------------------------%

            % Common modules



            powalg_obj = prima_mat.common.powalg_mod();


            % Inputs


            % BMAT(N, NPT + N)

            % SL(N)
            % SU(N)
            % XPT(N, NPT)
            % ZMAT(NPT, NPT-N-1)

            % Outputs
            d = NaN(size(xpt, 1), 1); % D(N)

            % Local variables


            ilbd = NaN;
            isbd = NaN(3, size(xpt, 2));

            iubd = NaN;
            k = NaN;


            mask_fixl = false(size(xpt, 1), 1);
            mask_fixu = false(size(xpt, 1), 1);


            curv = NaN;
            dderiv = NaN(size(xpt, 2), 1);


            distsq = NaN(size(xpt, 2), 1);

            glag = NaN(size(xpt, 1), 1);
            grdstp = NaN;
            gs = NaN;
            lfrac = NaN(size(xpt, 1), 1);
            pqlag = NaN(size(xpt, 2), 1);

            resis = NaN;
            s = NaN(size(xpt, 1), 1);
            scaling = NaN;
            sfixsq = NaN;
            slbd = NaN;
            slbd_test = NaN(size(xpt, 1), 1);
            ssqsav = NaN;
            stplen = NaN(3, size(xpt, 2));
            stpm = NaN;
            subd = NaN;
            subd_test = NaN(size(xpt, 1), 1);
            sumin = NaN;
            sxpt = NaN(size(xpt, 2), 1);
            ufrac = NaN(size(xpt, 1), 1);

            vlagsq = NaN;

            x = NaN(size(xpt, 1), 1);

            xdiff = NaN(size(xpt, 1), 1);
            xline = NaN(size(xpt, 1), 1);
            xopt = NaN(size(xpt, 1), 1);
            xtemp = NaN(size(xpt, 1), 1);

            % Sizes.
            n = size(xpt, 1);
            npt = size(xpt, 2);

            % Preconditions


            %====================%
            % Calculation starts %
            %====================%

            % PQLAG contains the leading NPT elements of the KNEW-th column of H, and it provides the second
            % derivative parameters of LFUNC, which is the KNEW-th Lagrange function. ALPHA will is the KNEW-th
            % diagonal element of the H matrix.
            pqlag(:) = zmat * zmat(knew, :).';
            alpha = pqlag(knew);

            % Read XOPT.
            xopt(:) = xpt(:, kopt);

            % Calculate the gradient GLAG of the KNEW-th Lagrange function at XOPT.
            glag(:) = bmat(:, knew) + powalg_obj.hess_mul(xopt, xpt, pqlag);

            % In case GLAG contains NaN, set D to a displacement from XOPT to XPT(:, KNEW) and return. Powell's
            % code does not have this, and D may be NaN in the end. Note that it is crucial to ensure that a
            % geometry step is nonzero.
            if ~isfinite(sum(abs(glag), 'all'))
                d(:) = xpt(:, knew) - xopt;
                d = min(0.5, delbar / norm(d)) * d; % Since XPT respects the bounds, so does XOPT + D.
                return
            end

            % Search for a large denominator along the straight lines through XOPT and another interpolation
            % point, subject to the bound constraints and the trust region. According to these constraints, SLBD
            % and SUBD will be lower and upper bounds on the step along each of these lines in turn. On each
            % line, we will evaluate the value of the KNEW-th Lagrange function at 3 trial points, and estimate
            % the denominator accordingly. The three points take the form (1-t)*XOPT + t*XPT(:, K) with step
            % lengths t = SLBD, SUBD, and STPM, corresponding to the upper (U) bound of t, the lower (L) bound
            % of t, and a medium (M) step. In total, 3*(NPT-1) trial points will be considered. On the K-th line,
            % we intend to maximize the modulus of PHI_K(t) = LFUNC((1-t)*XOPT + t*XPT(:,K)); overall, we intend
            % to find a trial point rendering a large value of the PREDSQ defined in (3.11) of the BOBYQA paper.
            %
            % We start with the following DO loop, the purpose of which is to define two 3-by-NPT arrays STPLEN
            % and ISBD. For each K, STPLEN(1:3, K) and ISBD(1:3, K) corresponds to the straight line through
            % XOPT and XPT(:, K). STPLEN(1:3, K) contains SLBD, SUBD, and STPM in this order, which are the step
            % lengths for the three trial points on this line. The three entries of SBDI(1:3, K) indicate
            % whether the corresponding trial points lie on bounds; SBDI(I, K) = J > 0 means that the I-th trail
            % point on the K-th line attains the J-th upper bound, SBDI(I, K) = -J < 0 indicates reaching the
            % J-th lower bound, and SBDI(I, K) = 0 means not touching any bound.
            dderiv(:) = xpt.' * glag - sum(glag .* xopt, 'all'); % The derivatives PHI_K'(0).
            distsq(:) = sum((xpt - xopt) .^ 2, 1);
            for k = 1:npt
                % It does not make sense to consider "straight line through XOPT and XPT(:, KOPT)". Hence set
                % STPLEN(:, KOPT) = 0 and ISBD(:, KOPT) = 0 so that VLAG(:, K) and PREDSQ(:, K) obtained after
                % this loop will be both zero and the search will skip K = KOPT. To avoid undesired/unpredictable
                % behavior due to possible NaN, set DDERIV(K) = 0 if K = KOPT or if DDERIV(K) is originally NaN.
                if k == kopt || isnan(dderiv(k))
                    dderiv(k) = 0.0;
                    stplen(:, k) = 0.0;
                    isbd(:, k) = 0;
                    continue
                end

                subd = delbar / sqrt(distsq(k)); % DISTSQ(K) > 0 unless K == KOPT or the input is incorrect.
                slbd = -subd;
                ilbd = 0;
                iubd = 0;
                sumin = min(1.0, subd);

                % Revise SLBD and SUBD if necessary because of the bounds in SL and SU according to LFRAC, UFRAC.
                % N.B.: We calculate LFRAC only at the positions where SL - XOPT > -ABS(XDIFF) * SUBD, because
                % the values of LFRAC are relevant only at the positions where ABS(LFRAC) < SUBD. Powell's code
                % does not check this inequality before evaluating LFRAC, and overflow may occur due to large
                % entries of SL. Note that SL - XOPT > -ABS(XDIFF) * SUBD implies that XDIFF /= 0, as long as
                % SL <= XOPT is ensured. In addition, when initializing LFRAC to SIGN(SUBD, -XDIFF), we do not
                % need to worry about the case where XDIFF = 0, because we only use LFRAC when XDIFF /= 0.
                % Similar things can be said about UFRAC.
                xdiff(:) = xpt(:, k) - xopt;
                lfrac = subd .* ((-xdiff > 0) .* 2 - 1);
                lfrac(sl - xopt > -abs(xdiff) * subd) = (sl(sl - xopt > -abs(xdiff) * subd) - xopt(sl - xopt > -abs(xdiff) * subd)) ./ xdiff(sl - xopt > -abs(xdiff) * subd);
                ufrac = subd .* ((xdiff > 0) .* 2 - 1);
                ufrac(su - xopt < abs(xdiff) * subd) = (su(su - xopt < abs(xdiff) * subd) - xopt(su - xopt < abs(xdiff) * subd)) ./ xdiff(su - xopt < abs(xdiff) * subd);
                %%MATLAB code for LFRAC and UFRAC (the code is simpler as we are not concerned about overflow):
                %%xdiff = xpt(:, k) - xopt;
                %%lfrac = (sl - xopt) / xdiff;
                %%ufrac = (su - xopt) / xdiff;

                % First, revise SLBD. Note that SLBD_TEST <= 0 unless the input violates XOPT >= SL.
                slbd_test(:) = slbd;
                slbd_test(xdiff > 0) = lfrac(xdiff > 0);
                slbd_test(xdiff < 0) = ufrac(xdiff < 0);
                if any(slbd_test > slbd, 'all')
                    [~, ilbd] = max(slbd_test, [], 'omitnan');
                    slbd = slbd_test(ilbd);
                    ilbd = -ilbd * round(1.0 .* ((xdiff(ilbd) > 0) .* 2 - 1));
                    %%MATLAB:
                    %%[slbd, ilbd] = max(slbd_test, [], 'omitnan');
                    %%ilbd = -ilbd * sign(xdiff(ilbd));

                end

                % Second, revise SUBD. Note that SUBD_TEST >= 0 unless the input violates XOPT <= SU.
                subd_test(:) = subd;
                subd_test(xdiff > 0) = ufrac(xdiff > 0);
                subd_test(xdiff < 0) = lfrac(xdiff < 0);
                if any(subd_test < subd, 'all')
                    [~, iubd] = min(subd_test, [], 'omitnan');
                    subd = max(sumin, subd_test(iubd));
                    iubd = iubd * round(1.0 .* ((xdiff(iubd) > 0) .* 2 - 1));
                    %%MATLAB:
                    %%[subd, iubd] = min(subd_test, [], 'omitnan');
                    %%subd = max(sumin, subd);
                    %%iubd = iubd * sign(xdiff(iubd));

                end


                % Now, define the step length STPM between SLBD and SUBD by finding the critical point of the
                % function PHI_K(t) = LFUNC((1-t)*XOPT + t*XPT(:,K)) mentioned above. It is a quadratic since
                % LFUNC is the KNEW-th Lagrange function. For K /= KNEW, the critical point is 0.5, as
                % PHI_K(0) = 1 = PHI_K(1); when K = KNEW, it is -0.5*PHI_K'(0) / (1 - PHI_K'(0)), because
                % PHI_K(0) = 0 and PHI_K(1) = 1.
                stpm = 0.5;
                if k == knew
                    stpm = slbd;
                    if abs(1.0 - dderiv(k)) > 0
                        stpm = -0.5 * dderiv(k) / (1.0 - dderiv(k));
                    end
                end
                stpm = max(slbd, min(subd, stpm));

                stplen(:, k) = [slbd, subd, stpm];
                isbd(:, k) = [ilbd, iubd, 0];
            end

            % The following lines calculate PREDSQ for all the 3*(NPT-1) trial points.
            % First, compute VLAG = PHI(STPLEN). Using the fact that PHI_K(0) = 0, PHI_K(1) = delta_{K, KNEW}
            % (Kronecker delta), and recalling the PHI_K is quadratic, we can find that
            % PHI_K(t) = t*(1-t)*PHI_K'(0) for K /= KNEW, and PHI_KNEW = t*[t*(1-PHI_K'(0)) + PHI_K'(0)].
            vlag = stplen .* (1.0 - stplen) .* dderiv.';
            %%MATLAB: vlag = stplen .* (1 - stplen) .* dderiv; % Implicit expansion; dderiv is a row!
            vlag(:, knew) = stplen(:, knew) .* (stplen(:, knew) * (1.0 - dderiv(knew)) + dderiv(knew));
            % Set NaNs in VLAG to 0 so that the behavior of MAXVAL(ABS(VLAG)) is predictable. VLAG does not have
            % NaN unless XPT does, which would be a bug. MAXVAL(ABS(VLAG)) appears in Powell's code, not here.
            vlag(isnan(vlag)) = 0.0; %%MATLAB: vlag(isnan(vlag)) = 0;
            %
            % Second, BETABD is the upper bound of BETA given in (3.10) of the BOBYQA paper.
            betabd = 0.5 * (stplen .* (1.0 - stplen) .* distsq.') .^ 2;
            %%MATLAB: betabd = 0.5 * (stplen .* (1-stplen) .* distsq).^2 % Implicit expansion; distsq is a row!
            %
            % Finally, PREDSQ is the quantity defined in (3.11) of the BOBYQA paper.
            predsq = vlag .* vlag .* (vlag .* vlag + alpha * betabd);
            % Set NaNs in PREDSQ to 0 so that the behavior of MAXLOC(PREDSQ) is predictable. PREDSQ does not
            % have NaN unless XPT does, which would be a bug.
            predsq(isnan(predsq)) = 0.0; %%MATLAB: predsq(isnan(predsq)) = 0

            % Locate the trial point the renders the maximum of PREDSQ. It is the ISQ-th trial point on the
            % straight line through XOPT and XPT(:, KSQ).
            % N.B.: 1. The strategy is a bit different from Powell's original code. In Powell's code and the
            % BOBYQA paper, we first select the trial point that gives the largest value of ABS(VLAG) on each
            % straight line, and then maximize PREDSQ among the (NPT-1) selected points. Here we maximize PREDSQ
            % among all the trial points. It works slightly better than Powell's version in a test on 20220428.
            % Powell's version is as follows.
            %---------------------------------------------------------------------%
            %isqs = int(maxloc(abs(vlag), dim=1), kind(isqs))  ! SIZE(ISQS) = NPT
            %ksq = int(maxloc([(predsq(isqs(k), k), k=1, npt)], dim=1), kind(ksq))
            %isq = isqs(ksq)
            %---------------------------------------------------------------------%
            % 2. Recall that we have set the NaN entries of PREDSQ to zero, if there is any. Thus the KSQS below
            % is a well defined integer array, all the three entries lying between 1 and NPT.
            ksqs = fortran.maxloc(predsq, 'dim', 2);
            isq = fortran.maxloc([predsq(1, ksqs(1)), predsq(2, ksqs(2)), predsq(3, ksqs(3))], 'dim', 1);
            ksq = ksqs(isq);
            %%MATLAB:
            %%[~, ksqs] = max(predsq, [], 'omitnan');
            %%[~, isq] = max([predsq(1, ksqs(1)), predsq(2, ksqs(2)), predsq(3, ksqs(3))]);
            %%ksq = ksqs(isq);

            % Construct XLINE in a way that satisfies the bound constraints exactly.
            stpsiz = stplen(isq, ksq);
            ibd = isbd(isq, ksq);

            xline(:) = max(sl, min(su, xopt + stpsiz * (xpt(:, ksq) - xopt)));
            if ibd < 0
                xline(-ibd) = sl(-ibd);
            end
            if ibd > 0
                xline(ibd) = su(ibd);
            end

            % Calculate DENOM for the current choice of D. Indeed, only DEN_LINE(KNEW) is needed.
            % Zaikun 20250907: It was observed numerically that D could be ZERO here (i.e., XLINE = XOPT).
            % Should this be impossible in theory?
            d = xline - xopt;
            den_line = powalg_obj.calden(kopt, bmat, d, xpt, zmat);

            %--------------------------------------------------------------------------------------------------%
            % The following IF ... END IF does not exist in Powell's code. SURPRISINGLY, the performance of
            % BOBYQA on bound constrained problems (but NOT unconstrained ones) is evidently improved by this IF
            % ... END IF, which means to try the Cauchy step only in the late stage of the algorithm, e.g., when
            % DELBAR is relatively small. WHY? In the following condition, 1.0E-2 works well if we use
            % DEN_CAUCHY to decide whether to take the Cauchy step; 1.0E-3 works well if we use VLAGSQ instead.
            % How to make this condition adaptive? A naive idea is to replace the thresholds to,
            % e.g.,1.0E-2*RHOBEG. However, in a test on 20220517, this adaptation worsened the performance. In
            % such a test, RHOBEG must take a value that is quite different from one. We tried RHOBEG = 0.9E-2.
            %if (delbar > 1.0E-3_RP) then
            %if (delbar > 1.0E-1_RP) then
            if delbar > 1.0e-2
                return
            end
            %--------------------------------------------------------------------------------------------------%

            % Prepare for the method that assembles the constrained Cauchy step in S. The sum of squares of the
            % fixed components of S is formed in SFIXSQ, and the free components of S are set to BIGSTP. When
            % UPHILL = 0, the method calculates the downhill version of XCAUCHY, which intends to minimize the
            % KNEW-th Lagrange function; when UPHILL = 1, it calculates the uphill version that intends to
            % maximize the Lagrange function.
            bigstp = delbar + delbar; % N.B.: In the sequel, S <= BIGSTP.
            xcauchy = xopt;
            vlagsq_cauchy = 0.0;
            for uphill = 0:1
                if uphill == 1
                    glag = -glag;
                end
                s(:) = 0.0;
                mask_free = (min(xopt - sl, glag) > 0 | max(xopt - su, glag) < 0);
                s(mask_free) = bigstp;
                ggfree = sum(glag(find(mask_free)) .^ 2, 'all');
                % In Powell's code, the subroutine returns immediately if GGFREE is 0. However, GGFREE depends
                % on GLAG, which in turn depends on UPHILL. It can happen that GGFREE is 0 when UPHILL = 0 but
                % not so when UPHILL= 1. Thus we skip the iteration for the current UPHILL but do not return.
                if ggfree <= 0 || isnan(ggfree)
                    continue
                end

                % Investigate whether more components of S can be fixed. Note that the loop counter K does not
                % appear in the loop body. The purpose of K is only to impose an explicit bound on the number of
                % loops. Powell's code does not have such a bound. The bound is not a true restriction, because
                % we can check that (SFIXSQ > SSQSAV .AND. GGFREE > 0) must fail within N loops.
                sfixsq = 0.0;
                grdstp = 0.0;
                for k = 1:n
                    resis = delbar ^ 2 - sfixsq;
                    if resis <= 0
                        break
                    end
                    ssqsav = sfixsq;
                    grdstp = sqrt(resis / ggfree);
                    xtemp = xopt - grdstp * glag;
                    mask_fixl = (s >= bigstp & xtemp <= sl); % S == BIGSTP & XTEMP == SL
                    mask_fixu = (s >= bigstp & xtemp >= su); % S == BIGSTP & XTEMP == SU
                    mask_free = (s >= bigstp & ~(mask_fixl | mask_fixu));
                    s(mask_fixl) = sl(mask_fixl) - xopt(mask_fixl);
                    s(mask_fixu) = su(mask_fixu) - xopt(mask_fixu);
                    sfixsq = sfixsq + sum(s(find(mask_fixl | mask_fixu)) .^ 2, 'all');
                    ggfree = sum(glag(find(mask_free)) .^ 2, 'all');
                    if ~(sfixsq > ssqsav && ggfree > 0)
                        break
                    end
                end

                % Set the remaining free components of S and all components of XCAUCHY. S may be scaled later.
                x(glag > 0) = sl(glag > 0);
                x(glag <= 0) = su(glag <= 0);
                x(abs(s) <= 0) = xopt(abs(s) <= 0);
                xtemp(:) = max(sl, min(su, xopt - grdstp * glag));
                x(s >= bigstp) = xtemp(s >= bigstp); % S == BIGSTP
                s(s >= bigstp) = -grdstp * glag(s >= bigstp); % S == BIGSTP
                gs = sum(glag .* s, 'all');

                % Set CURV to the curvature of the KNEW-th Lagrange function along S. Scale S by a factor less
                % than ONE if that can reduce the modulus of the Lagrange function at XOPT+S. Set CAUCHY to the
                % final value of the square of this function.
                sxpt(:) = xpt.' * s;
                curv = sum(sxpt .* (pqlag .* sxpt), 'all'); % CURV = INPROD(S, HESS_MUL(S, XPT, PQLAG))
                if uphill == 1
                    curv = -curv;
                end
                if curv > -gs && curv < -(1.0 + sqrt(2.0)) * gs
                    scaling = -gs / curv;
                    x(:) = max(sl, min(su, xopt + scaling * s));
                    vlagsq = (0.5 * gs * scaling) ^ 2;
                else
                    vlagsq = (gs + 0.5 * curv) ^ 2;
                end

                if vlagsq > vlagsq_cauchy
                    xcauchy = x;
                    vlagsq_cauchy = vlagsq;
                end
            end

            % Calculate the denominator rendered by the Cauchy step. Indeed, only DEN_CAUCHY(KNEW) is needed.
            s = xcauchy - xopt;
            den_cauchy = powalg_obj.calden(kopt, bmat, s, xpt, zmat);

            % Take the Cauchy step if it is likely to render a larger denominator.
            %IF (VLAGSQ_CAUCHY > MAX(DEN_LINE(KNEW), ZERO) .OR. IS_NAN(DEN_LINE(KNEW))) THEN  ! Powell's version
            if den_cauchy(knew) > max(den_line(knew), 0.0) || isnan(den_line(knew))
                % Works better
                d = s;
            end

            % In case D is zero or contains Inf/NaN, replace it with a displacement from XPT(:, KNEW) to XOPT.
            % Powell's code does not have this. Note that it is crucial to ensure that a geometry step is nonzero.
            if sum(abs(d), 'all') <= 0 || ~isfinite(sum(abs(d), 'all'))
                d(:) = xpt(:, knew) - xopt;
                d = min(0.5, delbar / norm(d)) * d; % Since XPT respects the bounds, so does XOPT + D.

            end

            %====================%
            %  Calculation ends  %
            %====================%

            % Postconditions


        end

    end
end