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
            consts_obj = consts_mod();
            debug_obj = debug_mod();
            infnan_obj = infnan_mod();
            linalg_obj = linalg_mod();
            powalg_obj = powalg_mod();


            % Inputs


            % BMAT(N, NPT + N)
            % D(N)


            % XPT(N, NPT)
            % ZMAT(NPT, NPT - N - 1)

            % Outputs
            knew = NaN;

            % Local variables
            srname = "SETDROP_TR";
            n = NaN;
            npt = NaN;
            den = NaN(size(xpt, 2), 1);
            distsq = NaN(size(xpt, 2), 1);
            score = NaN(size(xpt, 2), 1);
            weight = NaN(size(xpt, 2), 1);

            % Sizes
            n = size(xpt, 1);
            npt = size(xpt, 2);

            % Preconditions
            if consts_obj.DEBUGGING
                debug_obj.assert(n >= 1 && npt >= n + 2, "N >= 1, NPT >= N + 2", srname);
                debug_obj.assert(kopt >= 1 && kopt <= npt, "1 <= KOPT <= NPT", srname);
                debug_obj.assert(delta >= rho && rho > 0, "DELTA >= RHO > 0", srname);
                debug_obj.assert(numel(d) == n && all(infnan_obj.is_finite(d), 'all'), "SIZE(D) == N, D is finite", srname);
                debug_obj.assert(all(infnan_obj.is_finite(xpt), 'all'), "XPT is finite", srname);
                debug_obj.assert(size(bmat, 1) == n && size(bmat, 2) == npt + n, "SIZE(BMAT)==[N, NPT+N]", srname);
                debug_obj.assert(linalg_obj.issymmetric(bmat(:, npt + 1:npt + n)), "BMAT(:, NPT+1:NPT+N) is symmetric", srname);
                debug_obj.assert(size(zmat, 1) == npt && size(zmat, 2) == npt - n - 1, "SIZE(ZMAT) == [NPT, NPT - N - 1]", srname);
            end

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
                distsq(:) = sum(fortran.power((xpt - fortran.spread(xpt(:, kopt) + d, 'dim', 2, 'ncopies', npt)), 2), 1);
                %%MATLAB: distsq = sum((xpt - (xpt(:, kopt) + d)).^2)  % d should be a column! Implicit expansion

            else
                distsq(:) = sum(fortran.power((xpt - fortran.spread(xpt(:, kopt), 'dim', 2, 'ncopies', npt)), 2), 1);
                %%MATLAB: distsq = sum((xpt - xpt(:, kopt)).^2)  % Implicit expansion
            end

            weight(:) = fortran.power(max(consts_obj.ONE, distsq ./ fortran.power(rho, 2)), 4);
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
            den(:) = powalg_obj.calden(kopt, bmat, d, xpt, zmat);
            score(:) = weight .* den;

            % If the new F is not better than FVAL(KOPT), we set SCORE(KOPT) = -1 to avoid KNEW = KOPT.
            if ~ximproved
                score(kopt) = -consts_obj.ONE;
            end

            % SCORE(K) = NaN implies DEN(K) = NaN. We exclude such K as we want DEN to be big.
            score(linalg_obj.trueloc(infnan_obj.is_nan(score))) = -consts_obj.ONE;

            knew = 0;
            % The following IF works slightly better than `IF (ANY(SCORE > 0))` from Powell's BOBYQA/LINCOA code.
            if any(score > 1, 'all') || (ximproved && any(score > 0, 'all'))
                % Powell's UOBYQA and NEWUOA code.
                % See (6.1) of the BOBYQA paper for the definition of KNEW in this case.
                knew = fix(fortran.maxloc(score, 'dim', 1));
                %%MATLAB: [~, knew] = max(score);

            end

            % Powell's code does not include the following instructions. With Powell's code, if DEN consists of
            % only NaN, then KNEW can be 0 even when XIMPROVED is TRUE. Here, we set KNEW to the following value,
            % to make sure that the new trial point is included in the interpolation set. However, the updating
            % subroutine will likely need to skip the update of the Lagrange polynomials (i.e., H), or they
            % would be destroyed by the NaNs.
            if (ximproved && knew == 0) || knew < 0
                % KNEW < 0 is impossible in theory.
                knew = fix(fortran.maxloc(distsq, 'dim', 1));
            end

            %====================%
            %  Calculation ends  %
            %====================%

            % Postconditions
            if consts_obj.DEBUGGING
                debug_obj.assert(knew >= 0 && knew <= npt, "0 <= KNEW <= NPT", srname);
                debug_obj.assert(knew ~= kopt || ximproved, "KNEW /= KOPT unless XIMPROVED = TRUE", srname);
                debug_obj.assert(knew >= 1 || ~ximproved, "KNEW >= 1 unless XIMPROVED = FALSE", srname);
                % KNEW >= 1 when XIMPROVED = TRUE unless NaN occurs in DISTSQ, which should not happen if the
                % starting point does not contain NaN and the trust-region/geometry steps never contain NaN.

            end

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
            consts_obj = consts_mod();
            debug_obj = debug_mod();
            infnan_obj = infnan_mod();
            linalg_obj = linalg_mod();
            powalg_obj = powalg_mod();


            % Inputs


            % BMAT(N, NPT + N)

            % SL(N)
            % SU(N)
            % XPT(N, NPT)
            % ZMAT(NPT, NPT-N-1)

            % Outputs
            d = NaN(size(xpt, 1), 1); % D(N)

            % Local variables
            srname = "GEOSTEP";
            ibd = NaN;
            ilbd = NaN;
            isbd = NaN(3, size(xpt, 2));
            isq = NaN;
            iubd = NaN;
            k = NaN;
            ksq = NaN;
            ksqs = NaN(3, 1);
            n = NaN;
            npt = NaN;
            uphill = NaN;
            mask_fixl = false(size(xpt, 1), 1);
            mask_fixu = false(size(xpt, 1), 1);
            mask_free = false(size(xpt, 1), 1);
            alpha = NaN; stpsiz = NaN;
            betabd = NaN(3, size(xpt, 2));
            bigstp = NaN;
            curv = NaN;
            dderiv = NaN(size(xpt, 2), 1);
            den_cauchy = NaN(size(xpt, 2), 1);
            den_line = NaN(size(xpt, 2), 1);
            distsq = NaN(size(xpt, 2), 1);
            ggfree = NaN;
            glag = NaN(size(xpt, 1), 1);
            grdstp = NaN;
            gs = NaN;
            lfrac = NaN(size(xpt, 1), 1);
            pqlag = NaN(size(xpt, 2), 1);
            predsq = NaN(3, size(xpt, 2));
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
            vlag = NaN(3, size(xpt, 2));
            vlagsq = NaN;
            vlagsq_cauchy = NaN;
            x = NaN(size(xpt, 1), 1);
            xcauchy = NaN(size(xpt, 1), 1);
            xdiff = NaN(size(xpt, 1), 1);
            xline = NaN(size(xpt, 1), 1);
            xopt = NaN(size(xpt, 1), 1);
            xtemp = NaN(size(xpt, 1), 1);

            % Sizes.
            n = size(xpt, 1);
            npt = size(xpt, 2);

            % Preconditions
            if consts_obj.DEBUGGING
                debug_obj.assert(n >= 1 && npt >= n + 2, "N >= 1, NPT >= N + 2", srname);
                debug_obj.assert(knew >= 1 && knew <= npt, "1 <= KNEW <= NPT", srname);
                debug_obj.assert(kopt >= 1 && kopt <= npt, "1 <= KOPT <= NPT", srname);
                debug_obj.assert(knew ~= kopt, "KNEW /= KOPT", srname);
                debug_obj.assert(delbar > 0, "DELBAR > 0", srname);
                debug_obj.assert(numel(sl) == n && all(sl <= 0, 'all'), "SIZE(SL) == N, SL <= 0", srname);
                debug_obj.assert(numel(su) == n && all(su >= 0, 'all'), "SIZE(SU) == N, SU >= 0", srname);
                debug_obj.assert(all(infnan_obj.is_finite(xpt), 'all'), "XPT is finite", srname);
                debug_obj.assert(all(xpt >= fortran.spread(sl, 'dim', 2, 'ncopies', npt), 'all') && all(xpt <= fortran.spread(su, 'dim', 2, 'ncopies', npt), 'all'), "SL <= XPT <= SU", srname);
                debug_obj.assert(size(bmat, 1) == n && size(bmat, 2) == npt + n, "SIZE(BMAT) == [N, NPT+N]", srname);
                debug_obj.assert(linalg_obj.issymmetric(bmat(:, npt + 1:npt + n)), "BMAT(:, NPT+1:NPT+N) is symmetric", srname);
                debug_obj.assert(size(zmat, 1) == npt && size(zmat, 2) == npt - n - 1, "SIZE(ZMAT) == [NPT, NPT-N-1]", srname);
            end

            %====================%
            % Calculation starts %
            %====================%

            % PQLAG contains the leading NPT elements of the KNEW-th column of H, and it provides the second
            % derivative parameters of LFUNC, which is the KNEW-th Lagrange function. ALPHA will is the KNEW-th
            % diagonal element of the H matrix.
            pqlag(:) = linalg_obj.matprod21(zmat, zmat(knew, :));
            alpha = pqlag(knew);

            % Read XOPT.
            xopt(:) = xpt(:, kopt);

            % Calculate the gradient GLAG of the KNEW-th Lagrange function at XOPT.
            glag(:) = bmat(:, knew) + powalg_obj.hess_mul(xopt, xpt, pqlag);

            % In case GLAG contains NaN, set D to a displacement from XOPT to XPT(:, KNEW) and return. Powell's
            % code does not have this, and D may be NaN in the end. Note that it is crucial to ensure that a
            % geometry step is nonzero.
            if ~infnan_obj.is_finite(sum(abs(glag), 'all'))
                d(:) = xpt(:, knew) - xopt;
                d(:) = min(consts_obj.HALF, delbar / linalg_obj.p_norm(d)) * d; % Since XPT respects the bounds, so does XOPT + D.
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
            dderiv(:) = linalg_obj.matprod12(glag, xpt) - linalg_obj.inprod(glag, xopt); % The derivatives PHI_K'(0).
            distsq(:) = sum(fortran.power((xpt - fortran.spread(xopt, 'dim', 2, 'ncopies', npt)), 2), 1);
            for k = 1:npt
                % It does not make sense to consider "straight line through XOPT and XPT(:, KOPT)". Hence set
                % STPLEN(:, KOPT) = 0 and ISBD(:, KOPT) = 0 so that VLAG(:, K) and PREDSQ(:, K) obtained after
                % this loop will be both zero and the search will skip K = KOPT. To avoid undesired/unpredictable
                % behavior due to possible NaN, set DDERIV(K) = 0 if K = KOPT or if DDERIV(K) is originally NaN.
                if k == kopt || infnan_obj.is_nan_sp(dderiv(k))
                    dderiv(k) = consts_obj.ZERO;
                    stplen(:, k) = consts_obj.ZERO;
                    isbd(:, k) = 0;
                    continue
                end

                subd = delbar / fortran.sqrt(distsq(k)); % DISTSQ(K) > 0 unless K == KOPT or the input is incorrect.
                slbd = -subd;
                ilbd = 0;
                iubd = 0;
                sumin = min(consts_obj.ONE, subd);

                % Revise SLBD and SUBD if necessary because of the bounds in SL and SU according to LFRAC, UFRAC.
                % N.B.: We calculate LFRAC only at the positions where SL - XOPT > -ABS(XDIFF) * SUBD, because
                % the values of LFRAC are relevant only at the positions where ABS(LFRAC) < SUBD. Powell's code
                % does not check this inequality before evaluating LFRAC, and overflow may occur due to large
                % entries of SL. Note that SL - XOPT > -ABS(XDIFF) * SUBD implies that XDIFF /= 0, as long as
                % SL <= XOPT is ensured. In addition, when initializing LFRAC to SIGN(SUBD, -XDIFF), we do not
                % need to worry about the case where XDIFF = 0, because we only use LFRAC when XDIFF /= 0.
                % Similar things can be said about UFRAC.
                xdiff(:) = xpt(:, k) - xopt;
                lfrac(:) = fortran.sign(subd, -xdiff);
                lfrac(sl - xopt > -abs(xdiff) * subd) = (sl(sl - xopt > -abs(xdiff) * subd) - xopt(sl - xopt > -abs(xdiff) * subd)) ./ xdiff(sl - xopt > -abs(xdiff) * subd);
                ufrac(:) = fortran.sign(subd, xdiff);
                ufrac(su - xopt < abs(xdiff) * subd) = (su(su - xopt < abs(xdiff) * subd) - xopt(su - xopt < abs(xdiff) * subd)) ./ xdiff(su - xopt < abs(xdiff) * subd);
                %%MATLAB code for LFRAC and UFRAC (the code is simpler as we are not concerned about overflow):
                %%xdiff = xpt(:, k) - xopt;
                %%lfrac = (sl - xopt) / xdiff;
                %%ufrac = (su - xopt) / xdiff;

                % First, revise SLBD. Note that SLBD_TEST <= 0 unless the input violates XOPT >= SL.
                slbd_test = repmat(slbd, size(slbd_test));
                slbd_test(linalg_obj.trueloc(xdiff > 0)) = lfrac(linalg_obj.trueloc(xdiff > 0));
                slbd_test(linalg_obj.trueloc(xdiff < 0)) = ufrac(linalg_obj.trueloc(xdiff < 0));
                if any(slbd_test > slbd, 'all')
                    ilbd = fix(fortran.maxloc(slbd_test, 'mask', (~infnan_obj.is_nan(slbd_test)), 'dim', 1));
                    slbd = slbd_test(ilbd);
                    ilbd = -ilbd * round(fortran.sign(consts_obj.ONE, xdiff(ilbd)));
                    %%MATLAB:
                    %%[slbd, ilbd] = max(slbd_test, [], 'omitnan');
                    %%ilbd = -ilbd * sign(xdiff(ilbd));

                end

                % Second, revise SUBD. Note that SUBD_TEST >= 0 unless the input violates XOPT <= SU.
                subd_test = repmat(subd, size(subd_test));
                subd_test(linalg_obj.trueloc(xdiff > 0)) = ufrac(linalg_obj.trueloc(xdiff > 0));
                subd_test(linalg_obj.trueloc(xdiff < 0)) = lfrac(linalg_obj.trueloc(xdiff < 0));
                if any(subd_test < subd, 'all')
                    iubd = fix(fortran.minloc(subd_test, 'mask', (~infnan_obj.is_nan(subd_test)), 'dim', 1));
                    subd = max(sumin, subd_test(iubd));
                    iubd = iubd * round(fortran.sign(consts_obj.ONE, xdiff(iubd)));
                    %%MATLAB:
                    %%[subd, iubd] = min(subd_test, [], 'omitnan');
                    %%subd = max(sumin, subd);
                    %%iubd = iubd * sign(xdiff(iubd));

                end

                if consts_obj.DEBUGGING
                    debug_obj.assert(slbd <= 0 && subd >= 0, "SLBD <= 0 <= SUBD", srname);
                end

                % Now, define the step length STPM between SLBD and SUBD by finding the critical point of the
                % function PHI_K(t) = LFUNC((1-t)*XOPT + t*XPT(:,K)) mentioned above. It is a quadratic since
                % LFUNC is the KNEW-th Lagrange function. For K /= KNEW, the critical point is 0.5, as
                % PHI_K(0) = 1 = PHI_K(1); when K = KNEW, it is -0.5*PHI_K'(0) / (1 - PHI_K'(0)), because
                % PHI_K(0) = 0 and PHI_K(1) = 1.
                stpm = consts_obj.HALF;
                if k == knew
                    stpm = slbd;
                    if abs(consts_obj.ONE - dderiv(k)) > 0
                        stpm = -consts_obj.HALF * dderiv(k) / (consts_obj.ONE - dderiv(k));
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
            vlag(:, :) = stplen .* (consts_obj.ONE - stplen) .* fortran.spread(dderiv, 'dim', 1, 'ncopies', 3);
            %%MATLAB: vlag = stplen .* (1 - stplen) .* dderiv; % Implicit expansion; dderiv is a row!
            vlag(:, knew) = stplen(:, knew) .* (stplen(:, knew) * (consts_obj.ONE - dderiv(knew)) + dderiv(knew));
            % Set NaNs in VLAG to 0 so that the behavior of MAXVAL(ABS(VLAG)) is predictable. VLAG does not have
            % NaN unless XPT does, which would be a bug. MAXVAL(ABS(VLAG)) appears in Powell's code, not here.
            vlag(infnan_obj.is_nan(vlag)) = consts_obj.ZERO; %%MATLAB: vlag(isnan(vlag)) = 0;
            %
            % Second, BETABD is the upper bound of BETA given in (3.10) of the BOBYQA paper.
            betabd(:, :) = consts_obj.HALF * fortran.power((stplen .* (consts_obj.ONE - stplen) .* fortran.spread(distsq, 'dim', 1, 'ncopies', 3)), 2);
            %%MATLAB: betabd = 0.5 * (stplen .* (1-stplen) .* distsq).^2 % Implicit expansion; distsq is a row!
            %
            % Finally, PREDSQ is the quantity defined in (3.11) of the BOBYQA paper.
            predsq(:, :) = vlag .* vlag .* (vlag .* vlag + alpha * betabd);
            % Set NaNs in PREDSQ to 0 so that the behavior of MAXLOC(PREDSQ) is predictable. PREDSQ does not
            % have NaN unless XPT does, which would be a bug.
            predsq(infnan_obj.is_nan(predsq)) = consts_obj.ZERO; %%MATLAB: predsq(isnan(predsq)) = 0

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
            ksqs = repmat(fix(fortran.maxloc(predsq, 'dim', 2)), size(ksqs));
            isq = fix(fortran.maxloc([predsq(1, ksqs(1)), predsq(2, ksqs(2)), predsq(3, ksqs(3))], 'dim', 1));
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
            d(:) = xline - xopt;
            den_line(:) = powalg_obj.calden(kopt, bmat, d, xpt, zmat);

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
            xcauchy(:) = xopt;
            vlagsq_cauchy = consts_obj.ZERO;
            for uphill = 0:1
                if uphill == 1
                    glag(:) = -glag;
                end
                s = repmat(consts_obj.ZERO, size(s));
                mask_free(:) = (min(xopt - sl, glag) > 0 | max(xopt - su, glag) < 0);
                s(linalg_obj.trueloc(mask_free)) = bigstp;
                ggfree = sum(fortran.power(glag(linalg_obj.trueloc(mask_free)), 2), 'all');
                % In Powell's code, the subroutine returns immediately if GGFREE is 0. However, GGFREE depends
                % on GLAG, which in turn depends on UPHILL. It can happen that GGFREE is 0 when UPHILL = 0 but
                % not so when UPHILL= 1. Thus we skip the iteration for the current UPHILL but do not return.
                if ggfree <= 0 || infnan_obj.is_nan_sp(ggfree)
                    continue
                end

                % Investigate whether more components of S can be fixed. Note that the loop counter K does not
                % appear in the loop body. The purpose of K is only to impose an explicit bound on the number of
                % loops. Powell's code does not have such a bound. The bound is not a true restriction, because
                % we can check that (SFIXSQ > SSQSAV .AND. GGFREE > 0) must fail within N loops.
                sfixsq = consts_obj.ZERO;
                grdstp = consts_obj.ZERO;
                for k = 1:n
                    resis = fortran.power(delbar, 2) - sfixsq;
                    if resis <= 0
                        break
                    end
                    ssqsav = sfixsq;
                    grdstp = fortran.sqrt(resis / ggfree);
                    xtemp(:) = xopt - grdstp * glag;
                    mask_fixl(:) = (s >= bigstp & xtemp <= sl); % S == BIGSTP & XTEMP == SL
                    mask_fixu(:) = (s >= bigstp & xtemp >= su); % S == BIGSTP & XTEMP == SU
                    mask_free(:) = (s >= bigstp & ~(mask_fixl | mask_fixu));
                    s(linalg_obj.trueloc(mask_fixl)) = sl(linalg_obj.trueloc(mask_fixl)) - xopt(linalg_obj.trueloc(mask_fixl));
                    s(linalg_obj.trueloc(mask_fixu)) = su(linalg_obj.trueloc(mask_fixu)) - xopt(linalg_obj.trueloc(mask_fixu));
                    sfixsq = sfixsq + sum(fortran.power(s(linalg_obj.trueloc(mask_fixl | mask_fixu)), 2), 'all');
                    ggfree = sum(fortran.power(glag(linalg_obj.trueloc(mask_free)), 2), 'all');
                    if ~(sfixsq > ssqsav && ggfree > 0)
                        break
                    end
                end

                % Set the remaining free components of S and all components of XCAUCHY. S may be scaled later.
                x(linalg_obj.trueloc(glag > 0)) = sl(linalg_obj.trueloc(glag > 0));
                x(linalg_obj.trueloc(glag <= 0)) = su(linalg_obj.trueloc(glag <= 0));
                x(linalg_obj.trueloc(abs(s) <= 0)) = xopt(linalg_obj.trueloc(abs(s) <= 0));
                xtemp(:) = max(sl, min(su, xopt - grdstp * glag));
                x(linalg_obj.trueloc(s >= bigstp)) = xtemp(linalg_obj.trueloc(s >= bigstp)); % S == BIGSTP
                s(linalg_obj.trueloc(s >= bigstp)) = -grdstp * glag(linalg_obj.trueloc(s >= bigstp)); % S == BIGSTP
                gs = linalg_obj.inprod(glag, s);

                % Set CURV to the curvature of the KNEW-th Lagrange function along S. Scale S by a factor less
                % than ONE if that can reduce the modulus of the Lagrange function at XOPT+S. Set CAUCHY to the
                % final value of the square of this function.
                sxpt(:) = linalg_obj.matprod12(s, xpt);
                curv = linalg_obj.inprod(sxpt, pqlag .* sxpt); % CURV = INPROD(S, HESS_MUL(S, XPT, PQLAG))
                if uphill == 1
                    curv = -curv;
                end
                if curv > -gs && curv < -(consts_obj.ONE + fortran.sqrt(consts_obj.TWO)) * gs
                    scaling = -gs / curv;
                    x(:) = max(sl, min(su, xopt + scaling * s));
                    vlagsq = fortran.power((consts_obj.HALF * gs * scaling), 2);
                else
                    vlagsq = fortran.power((gs + consts_obj.HALF * curv), 2);
                end

                if vlagsq > vlagsq_cauchy
                    xcauchy(:) = x;
                    vlagsq_cauchy = vlagsq;
                end
            end

            % Calculate the denominator rendered by the Cauchy step. Indeed, only DEN_CAUCHY(KNEW) is needed.
            s(:) = xcauchy - xopt;
            den_cauchy(:) = powalg_obj.calden(kopt, bmat, s, xpt, zmat);

            % Take the Cauchy step if it is likely to render a larger denominator.
            %IF (VLAGSQ_CAUCHY > MAX(DEN_LINE(KNEW), ZERO) .OR. IS_NAN(DEN_LINE(KNEW))) THEN  ! Powell's version
            if den_cauchy(knew) > max(den_line(knew), consts_obj.ZERO) || infnan_obj.is_nan_sp(den_line(knew))
                % Works better
                d(:) = s;
            end

            % In case D is zero or contains Inf/NaN, replace it with a displacement from XPT(:, KNEW) to XOPT.
            % Powell's code does not have this. Note that it is crucial to ensure that a geometry step is nonzero.
            if sum(abs(d), 'all') <= 0 || ~infnan_obj.is_finite(sum(abs(d), 'all'))
                d(:) = xpt(:, knew) - xopt;
                d(:) = min(consts_obj.HALF, delbar / linalg_obj.p_norm(d)) * d; % Since XPT respects the bounds, so does XOPT + D.

            end

            %====================%
            %  Calculation ends  %
            %====================%

            % Postconditions
            if consts_obj.DEBUGGING
                debug_obj.assert(numel(d) == n, "SIZE(D) == N", srname);
                debug_obj.assert(all(infnan_obj.is_finite(d), 'all'), "D is finite", srname);
                % In theory, ||D|| <= DELBAR, which may be false due to rounding, but ||D|| >= 2*DELBAR is unlikely.
                % It is crucial to ensure that the geometry step is nonzero, which holds in theory. However, due
                % to the bound constraints, ||D|| may be much smaller than DELBAR.
                debug_obj.assert(linalg_obj.p_norm(d) > 0 && linalg_obj.p_norm(d) < consts_obj.TWO * delbar, "0 < ||D|| < 2*DELBAR", srname);
                % D is supposed to satisfy the bound constraints SL <= XOPT + D <= SU.
                debug_obj.assert(all(xopt + d >= sl - consts_obj.TEN * consts_obj.EPS_custom * max(consts_obj.ONE, abs(sl)) & xopt + d <= su + consts_obj.TEN * consts_obj.EPS_custom * max(consts_obj.ONE, abs(su)), 'all'), "SL <= XOPT + D <= SU", srname);
            end

        end

    end
end