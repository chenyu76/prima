classdef geometry_uobyqa_mod
    %--------------------------------------------------------------------------------------------------%
    % This module contains subroutines concerning the geometry-improving of the interpolation set XPT.
    %
    % Coded by Zaikun ZHANG (www.zhangzk.net) based on Powell's code and the UOBYQA paper.
    %
    % Dedicated to the late Professor M. J. D. Powell FRS (1936--2015).
    %
    % Started: February 2022
    %
    % Last Modified: Tue 10 Feb 2026 02:43:25 PM CET
    %--------------------------------------------------------------------------------------------------%

    methods
        function knew = setdrop_tr(~, kopt, ximproved, d, pl, rho, xpt)
            %--------------------------------------------------------------------------------------------------%
            % This subroutine sets KNEW to the index of the interpolation point to be deleted AFTER A TRUST
            % REGION STEP. KNEW will be set in a way ensuring that the geometry of XPT is "optimal" after
            % XPT(:, KNEW) is replaced with XNEW = XOPT + D, where D is the trust-region step. See the
            % discussions around (56) of the UOBYQA paper.
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


            powalg_obj = prima_mat.common.powalg_mod();

            % D(N)
            % PL(NPT-1, NPT)

            % XPT(N, NPT)


            knew = NaN;

            distsq = NaN(size(xpt, 2), 1);

            vlag = NaN(size(xpt, 2), 1);

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
            % %weight = max(ONE, distsq / rho**2)**3.5_RP ! ! No better than power 4.
            % %weight = max(ONE, distsq / delta**2)**3.5_RP  ! Not better than DISTSQ/RHO**2.
            % %weight = max(ONE, distsq / rho**2)**1.5_RP  ! Powell's origin code: power 1.5.
            % %weight = max(ONE, distsq / rho**2)**2  ! Better than power 1.5.
            % %weight = max(ONE, distsq / delta**2)**2  ! Not better than DISTSQ/RHO**2.
            % %weight = max(ONE, distsq / max(TENTH * delta, rho)**2)**2  ! The same as DISTSQ/RHO**2.
            % %weight = distsq**2  ! Not better than MAX(ONE, DISTSQ/RHO**2)**2
            % %weight = max(ONE, distsq / rho**2)**3  ! Better than power 2.
            % %weight = max(ONE, distsq / delta**2)**3  ! Similar to DISTSQ/RHO**2; not better than it.
            % %weight = max(ONE, distsq / max(TENTH * delta, rho)**2)**3  ! The same as DISTSQ/RHO**2.
            % %weight = distsq**3  ! Not better than MAX(ONE, DISTSQ/RHO**2)**3
            % %weight = max(ONE, distsq / delta**2)**4  ! Not better than DISTSQ/RHO**2.
            % %weight = max(ONE, distsq / max(TENTH * delta, rho)**2)**4  ! The same as DISTSQ/RHO**2.
            % %weight = distsq**4  ! Not better than MAX(ONE, DISTSQ/RHO**2)**4

            % Here, VLAG is the counterpart of DEN in NEWUOA/BOBYQA/LINCOA, representing the denominator in the
            % update of the Lagrange functions (or, inverse of the coefficient or KKT matrix of the
            % interpolation system).
            vlag(:) = powalg_obj.calvlag_qint(pl, d, xpt(:, kopt), kopt);
            score = weight .* abs(vlag);

            % If the new F is not better than FVAL(KOPT), we set SCORE(KOPT) = -1 to avoid KNEW = KOPT.
            if ~ximproved
                score(kopt) = -1.0;
            end

            % SCORE(K) is NaN implies VLAG(K) is NaN, but we want ABS(VLAG) to be big. So we exclude such K.
            score(isnan(score)) = -1.0;

            knew = 0;
            % It makes almost no difference if we change the IF below to `IF (ANY(SCORE>0))`, which is used
            % in Powell's BOBYQA and LINCOA code.
            if any(score > 1, 'all') || (ximproved && any(score > 0, 'all'))
                % Powell's UOBYQA and NEWUOA code
                [~, knew] = max(score);
                %%MATLAB: [~, knew] = max(score);

            end

            % Powell's code does not include the following instructions. With Powell's code, if VLAG consists of
            % only NaN, then KNEW can be 0 even when XIMPROVED is TRUE. Here, we set KNEW to the following value,
            % to make sure that the new trial point is included in the interpolation set. However, the updating
            % subroutine will likely need to skip the update of the Lagrange polynomial, or they would be
            % destroyed by the NaNs.
            if (ximproved && knew == 0) || knew < 0
                % KNEW < 0 is impossible in theory.
                [~, knew] = max(distsq);
            end

            %====================%
            %  Calculation ends  %
            %====================%


        end
        function d = geostep(~, knew, kopt, delbar, pl, xpt)
            %--------------------------------------------------------------------------------------------------%
            % This function calculates a step D that approximately solves
            %
            % maximize |LFUNC(XOPT + D)| subject to ||D|| <= DELBAR,
            %
            % so that the geometry of the interpolation step is improved when XPT(:, KNEW) becomes XOPT + D.
            % Here, LFUNC is the Lagrange polynomial at the KNEW-th interpolation point. See (8) and Section 2
            % of the UOBYQA paper.
            %
            % PL contains the parameters of the Lagrange functions. PL(KNEW) defines LFUNC.
            %
            % Powell's comments on the cost of linear algebra is as follows.
            % Calculating the D that maximizes |LFUNC(XOPT + D)| subject to ||D|| <= DELBAR requires of order
            % N^3 operations, but sometimes it is adequate if |LFUNC(XOPT + D)| is within about 0.9 of its
            % greatest possible value. This subroutine provides such a solution in only of order N^2 operations,
            % where the claim of accuracy has been tested by numerical experiments.
            %
            % N.B.: In Powell's UOBYQA code, DELBAR = RHO. We take the DELBAR of NEWUOA, which works better.
            %--------------------------------------------------------------------------------------------------%


            linalg_obj = prima_mat.common.linalg_mod();
            powalg_obj = prima_mat.common.powalg_mod();

            % PL(NPT-1, NPT)
            % XPT(N, NPT)


            d = NaN(size(xpt, 1), 1); % D(N)


            dcauchy = NaN(size(xpt, 1), 1);

            g = NaN(size(xpt, 1), 1);

            h = NaN(size(xpt, 1));
            hv = NaN(size(xpt, 1), 1);

            vlag = NaN(size(xpt, 2), 1);
            vlagc = NaN(size(xpt, 2), 1);

            xopt = NaN(size(xpt, 1), 1);

            n = size(xpt, 1);
            npt = size(xpt, 2);

            %====================%
            % Calculation starts %
            %====================%

            % Read XOPT.
            xopt(:) = xpt(:, kopt);

            % For the KNEW-th Lagrange function, evaluate the gradient at XOPT and the Hessian.
            g(:) = pl(1:n, knew) + linalg_obj.smat_mul_vec(pl(n + 1:npt - 1, knew), xopt);
            h(:, :) = linalg_obj.vec2smat(pl(n + 1:npt - 1, knew));

            % Evaluate GG = G^T*G and GHG = G^T*H*G. They will be used later.
            gg = sum(g .^ 2, 'all');
            ghg = sum(g .* (h * g), 'all');

            % Calculate the Cauchy step as a backup. Powell's code does not have this, and D may be 0 or NaN.
            if gg > 0 && isfinite(gg)
                dcauchy = (delbar / sqrt(gg)) * g;
                if ghg < 0
                    dcauchy = -dcauchy;
                end
            else                % GG is 0 or NaN due to rounding errors. Set DCAUCHY to a displacement from XOPT to XPT(:, KNEW).
                dcauchy(:) = xpt(:, knew) - xopt;
                scaling = delbar / norm(dcauchy);
                dcauchy = max(0.6 * scaling, min(0.5, scaling)) * dcauchy; % 0.6: ensure |D| > DELBAR/2
                if sum(g .* dcauchy, 'all') * sum(dcauchy .* (h * dcauchy), 'all') < 0
                    dcauchy = -dcauchy;
                end
            end

            % Return if H or G contains NaN or H is zero. Powell's code does not do this.
            if any(isnan(h), 'all') || any(isnan(g), 'all') || all(abs(h) <= 0, 'all')
                d = dcauchy;
                return
            end

            % Handle the case with N = 1. This should be done after the case where G or H contains NaN. Powell's
            % code does not contain this part.
            if n == 1
                if g(1) * h(1, 1) > 0
                    d(:) = delbar;
                else
                    d(:) = -delbar;
                end
                return
            end

            % Pick V such that ||HV|| / ||V|| is large.
            v = h(:, fortran.maxloc(sum(h .^ 2, 1), 'dim', 1));
            % Normalize V. Powell's code does not do this. It does not change the algorithm as only its
            % direction matters. It slightly improves the performance in the noiseless case.
            v = v ./ norm(v);

            % Set D to a vector in the subspace span{V, HV} that maximizes |(D, HD)|/(D, D), except that we set
            % D = HV if V and HV are nearly parallel.
            vv = sum(v .^ 2, 'all');
            d(:) = h * v;
            vhv = sum(v .* d, 'all');
            if vhv * vhv <= 0.9999 * sum(d .^ 2, 'all') * vv
                d = d - (vhv / vv) * v;
                dd = sum(d .^ 2, 'all');
                scaling = sqrt(dd / vv);
                dhd = sum(d .* (h * d), 'all');
                v = scaling * v;
                vhv = scaling * scaling * vhv;
                vhd = scaling * dd;
                temp = 0.5 * (dhd - vhv);
                if dhd + vhv < 0
                    d = vhd * v + (temp - sqrt(temp ^ 2 + vhd ^ 2)) * d;
                else
                    d = vhd * v + (temp + sqrt(temp ^ 2 + vhd ^ 2)) * d;
                end
            end

            % We now turn our attention to the subspace span{G, D}. A multiple of the current D is returned if
            % that choice seems to be adequate.
            dd = sum(d .^ 2, 'all');
            gd = sum(g .* d, 'all');
            dhd = sum(d .* (h * d), 'all');

            % Zaikun 20220504: GG and DD can become 0 at this point due to rounding. Detected by IFORT.
            if ~(gg > 0 && dd > 0)
                d = dcauchy;
                return
            end

            v = d - (gd / gg) * g;
            vv = sum(v .^ 2, 'all');
            if gd * dhd < 0
                scaling = -delbar / sqrt(dd);
            else
                scaling = delbar / sqrt(dd);
            end
            d = scaling * d;
            gnorm = sqrt(gg);

            if ~(gnorm * dd > 5.0e-3 * delbar * abs(dhd) && vv > 1.0e-4 * dd)
                % It may happen that D = 0 due to overflow in DD, which is used to define SCALING.
                if sum(abs(d), 'all') <= 0 || ~isfinite(sum(abs(d), 'all'))
                    d = dcauchy;
                end
                return
            end

            % G and V are now orthogonal in the subspace span{G, D}. Hence we generate an orthonormal basis of
            % this subspace such that (D, HV) is negligible or 0, where D and V will be the basis vectors.
            hv(:) = h * v;
            vhg = sum(g .* hv, 'all');
            vhv = sum(v .* hv, 'all');
            vnorm = sqrt(vv);
            ghg = ghg / gg;
            vhg = vhg / (vnorm * gnorm);
            vhv = vhv / vv;
            if abs(vhg) <= 1.0e-2 * max(abs(ghg), abs(vhv))
                vmu = ghg - vhv;
                wcos = 1.0;
                wsin = 0.0;
            else
                temp = 0.5 * (ghg - vhv);
                if temp < 0
                    vmu = temp - sqrt(temp ^ 2 + vhg ^ 2);
                else
                    vmu = temp + sqrt(temp ^ 2 + vhg ^ 2);
                end
                temp = sqrt(vmu ^ 2 + vhg ^ 2);
                wcos = vmu / temp;
                wsin = vhg / temp;
            end
            tempa = wcos / gnorm;
            tempb = wsin / vnorm;
            tempc = wcos / vnorm;
            tempd = wsin / gnorm;
            d = tempa * g + tempb * v;
            v = tempc * v - tempd * g;

            % The final D is a multiple of the current D, V, D + V or D - V. We make the choice from these
            % possibilities that is optimal.
            dlin = wcos * gnorm / delbar;
            vlin = -wsin * gnorm / delbar;
            tempa = abs(dlin) + 0.5 * abs(vmu + vhv);
            tempb = abs(vlin) + 0.5 * abs(ghg - vmu);
            tempc = sqrt(0.5) * (abs(dlin) + abs(vlin)) + 0.25 * abs(ghg + vhv);
            if tempa >= tempb && tempa >= tempc
                if dlin * (vmu + vhv) < 0
                    tempd = -delbar;
                else
                    tempd = delbar;
                end
                tempv = 0.0;
            elseif tempb >= tempc
                tempd = 0.0;
                if vlin * (ghg - vmu) < 0
                    tempv = -delbar;
                else
                    tempv = delbar;
                end
            else
                if dlin * (ghg + vhv) < 0
                    tempd = -sqrt(0.5) * delbar;
                else
                    tempd = sqrt(0.5) * delbar;
                end
                if vlin * (ghg + vhv) < 0
                    tempv = -sqrt(0.5) * delbar;
                else
                    tempv = sqrt(0.5) * delbar;
                end
            end
            d = tempd * d + tempv * v;

            % Replace D with DCAUCHY if needed. Powell's code does not have this part. Indeed, only the KNEW-th
            % entries of VLAG and VLAGC are needed.
            vlag(:) = powalg_obj.calvlag_qint(pl, d, xopt, kopt);
            vlagc(:) = powalg_obj.calvlag_qint(pl, dcauchy, xopt, kopt);
            if abs(vlagc(knew)) > 2.0 * abs(vlag(knew)) || isnan(vlag(knew))
                d = dcauchy;
            end

            %====================%
            %  Calculation ends  %
            %====================%


        end

    end
end