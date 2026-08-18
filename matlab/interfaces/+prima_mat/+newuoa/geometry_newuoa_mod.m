classdef geometry_newuoa_mod
    %--------------------------------------------------------------------------------------------------%
    % This module contains subroutines concerning the geometry-improving of the interpolation set XPT.
    %
    % Coded by Zaikun ZHANG (www.zhangzk.net) based on Powell's code and the NEWUOA paper.
    %
    % Dedicated to the late Professor M. J. D. Powell FRS (1936--2015).
    %
    % Started: July 2020
    %
    % Last Modified: Sunday, April 21, 2024 PM03:20:57
    %--------------------------------------------------------------------------------------------------%

    methods
        function knew = setdrop_tr(~, idz, kopt, ximproved, bmat, d, delta, rho, xpt, zmat)
            %--------------------------------------------------------------------------------------------------%
            % This subroutine sets KNEW to the index of the interpolation point to be deleted AFTER A TRUST
            % REGION STEP. KNEW will be set in a way ensuring that the geometry of XPT is "optimal" after
            % XPT(:, KNEW) is replaced with XNEW = XOPT + D, where D is the trust-region step.
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
                %%MATLAB: distsq = sum((xpt - (xpt(:, kopt) + d)).^2)  % d should be a column! Implicit expansion

            else
                distsq(:) = sum((xpt - xpt(:, kopt)) .^ 2, 1);
                %%MATLAB: distsq = sum((xpt - xpt(:, kopt)).^2)  % Implicit expansion
            end

            weight = max(1.0, distsq ./ max(0.1 * delta, rho) ^ 2) .^ 3; % Powell's code.
            % Other possible definitions of WEIGHT.
            % %weight = max(ONE, distsq / max(TENTH * delta, rho)**2)**3.5  ! This sometimes works better
            % %weight = max(ONE, distsq / rho**2)**3  ! This works almost the same as Powell's code
            % %weight = max(ONE, distsq / delta**2)**3  ! BOBYQA code. It does not work as well as the above.
            % %weight = max(ONE, distsq / max(TENTH * delta, rho)**2)**4  ! It does not work as well as the above.
            % %weight = max(ONE, distsq / max(TENTH * delta, rho)**2)**2  ! This works poorly.

            den = powalg_obj.calden(kopt, bmat, d, xpt, zmat, 'idz', idz);
            score = weight .* abs(den);

            % If the new F is not better than FVAL(KOPT), we set SCORE(KOPT) = -1 to avoid KNEW = KOPT.
            if ~ximproved
                score(kopt) = -1.0;
            end

            % SCORE(K) is NaN implies ABS(DEN(K)) is NaN, but we want ABS(DEN) to be big. So we exclude such K.
            score(isnan(score)) = -1.0;

            knew = 0;
            % The following IF works a bit better than `IF (ANY(SCORE > 0))` from Powell's BOBYQA/LINCOA code.
            if any(score > 1, 'all') || (ximproved && any(score > 0, 'all'))
                % Powell's UOBYQA and NEWUOA code
                % See (7.5) of the NEWUOA paper for the definition of KNEW in this case.
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
        function d = geostep(obj, idz, knew, kopt, bmat, delbar, xpt, zmat)
            %--------------------------------------------------------------------------------------------------%
            % This subroutine finds a step D that intends to improve the geometry of the interpolation set when
            % XPT(:, KNEW) is changed to XOPT + D, where XOPT = XPT(:, KOPT).
            %
            % XPT contains the current interpolation points.
            % BMAT provides the last N ROWs of H.
            % ZMAT and IDZ give a factorization of the first NPT by NPT sub-matrix of H.
            % KNEW is the index of the interpolation point to be dropped.
            % DELBAR is the trust region bound for the geometry step
            % D will be set to the step from X to the new point.
            %--------------------------------------------------------------------------------------------------%



            powalg_obj = prima_mat.common.powalg_mod();


            % BMAT(N, NPT + N)

            % XPT(N, NPT)
            % ZMAT(NPT, NPT - N - 1)


            d = NaN(size(xpt, 1), 1); % D(N)



            pqlag = NaN(size(xpt, 2), 1);


            %====================%
            % Calculation starts %
            %====================%

            d = obj.biglag(idz, knew, bmat, delbar, xpt(:, kopt), xpt, zmat);

            % PQLAG contains the leading NPT elements of the KNEW-th column of H, and it provides the second
            % derivative parameters of LFUNC.
            pqlag(:) = powalg_obj.omega_col(idz, zmat, knew);
            alpha = pqlag(knew); % ALPHA is the KNEW-th diagonal entry of H, i.e., that of Omega.

            % Calculate VLAG and BETA for D. Indeed, only VLAG(KNEW) is needed.
            vlag = powalg_obj.calvlag_lfqint(kopt, bmat, d, xpt, zmat, 'idz', idz);
            beta = powalg_obj.calbeta(kopt, bmat, d, xpt, zmat, 'idz', idz);
            denom = alpha * beta + vlag(knew) ^ 2;

            % If the cancellation in DENOM is unacceptable, then BIGDEN calculates an alternative model step D.
            % As in (6.17) of the NEWUOA paper, DENRAT = |ALPHA*BETA + TAU^2| / TAU^2 with TAU = VLAG(KNEW).
            % Powell's code does not check whether VLAG(KNEW)**2 > 0, which holds in theory. VLAG(KNEW) can
            % become 0 or NaN numerically, which did happen in tests, indicating a failure of BIGLAG, because
            % BIGLAG should maximize |VLAG(KNEW)|. Upon this failure, it is reasonable to call BIGDEN. For the
            % same reason, we check whether BETA is NaN. Why not check ALPHA? Because BIGDEN cannot improve ALPHA.
            % Powell's code takes DDEN once it is calculated. We take it only if it renders a bigger denominator.
            denrat = -1.0;
            if vlag(knew) ^ 2 > 0 && ~isnan(beta)
                denrat = abs(1.0 + alpha * beta / vlag(knew) ^ 2);
            end
            % If DENRAT is NaN at this point, then ALPHA is NaN, and there is no need to call BIGDEN.
            if denrat <= 0.8
                dden = obj.bigden(idz, knew, kopt, bmat, d, xpt, zmat);
                vlag = powalg_obj.calvlag_lfqint(kopt, bmat, dden, xpt, zmat, 'idz', idz);
                beta = powalg_obj.calbeta(kopt, bmat, dden, xpt, zmat, 'idz', idz);
                if abs(alpha * beta + vlag(knew) ^ 2) >= abs(denom) || isnan(denom)
                    d = dden;
                end
            end

            % In case D is zero or contains Inf/NaN, replace it with a displacement from XPT(:, KNEW) to
            % XOPT. Powell's code does not have this.
            if sum(abs(d), 'all') <= 0 || ~isfinite(sum(abs(d), 'all'))
                d(:) = xpt(:, knew) - xpt(:, kopt);
                scaling = delbar / norm(d);
                d = max(0.6 * scaling, min(0.5, scaling)) * d; % 0.6: ensure |D| > DELBAR/2

            end

            %====================%
            %  Calculation ends  %
            %====================%



        end
        function d = biglag(obj, idz, knew, bmat, delbar, x, xpt, zmat)
            %--------------------------------------------------------------------------------------------------%
            % This subroutine calculates a D by approximately solving
            %
            % max |LFUNC(X + D)|, subject to ||D|| <= DELBAR,
            %
            % where LFUNC is the KNEW-th Lagrange function. See Section 6 of the NEWUOA paper.
            %--------------------------------------------------------------------------------------------------%



            linalg_obj = prima_mat.common.linalg_mod();
            powalg_obj = prima_mat.common.powalg_mod();
            univar_obj = prima_mat.common.univar_mod();


            % BMAT(N, NPT + N)

            % X(N)
            % XPT(N, NPT)
            % ZMAT(NPT, NPT - N - 1)


            d = NaN(size(xpt, 1), 1); % D(N)



            angle = NaN;
            cf = NaN(5, 1);
            cth = NaN;


            dold = NaN(numel(x), 1);
            gc = NaN(numel(x), 1);


            pqlag = NaN(size(xpt, 2), 1);
            s = NaN(numel(x), 1);


            sth = NaN;

            % LFUNC(X)

            w = NaN(numel(x), 1);


            n = size(xpt, 1);


            %====================%
            % Calculation starts %
            %====================%

            % PQLAG contains the leading NPT elements of the KNEW-th column of H, and it provides the second
            % derivative parameters of LFUNC.
            pqlag(:) = powalg_obj.omega_col(idz, zmat, knew);

            % Set the unscaled initial D. Form the gradient of LFUNC at X, and multiply D by the Hessian of LFUNC.
            d(:) = xpt(:, knew) - x;
            dd = sum(d .* d, 'all');
            gd = powalg_obj.hess_mul(d, xpt, pqlag); % GD = MATPROD(XPT, PQLAG * MATPROD(D, XPT))

            gc(:) = bmat(:, knew) + powalg_obj.hess_mul(x, xpt, pqlag); % GC = BMAT(:,KNEW) + MATPROD(XPT,PQLAG*MATPROD(X,XPT))

            % Scale D and GD, with a sign change if needed. Set S to another vector in the initial 2-D subspace.
            gg = sum(gc .* gc, 'all');
            sp = sum(d .* gc, 'all');
            dhd = sum(d .* gd, 'all');
            scaling = delbar / sqrt(dd);
            if sp * dhd < 0
                scaling = -scaling;
            end
            t = 0.0;
            if sp ^ 2 > 0.99 * dd * gg
                t = 1.0;
            end
            tau = scaling * (abs(sp) + 0.5 * scaling * abs(dhd));
            if gg * delbar ^ 2 < 1.0e-2 * tau ^ 2
                t = 1.0;
            end
            if isfinite(sum(abs(scaling * d), 'all'))
                d = scaling * d;
                gd = scaling * gd;
                s = gc + t * gd;
                maxiter = n;
            else
                maxiter = 0; % Return immediately to avoid producing a D containing NaN/Inf.
            end

            tol = min(0.1, max(eps(1.0) ^ 0.25, 1.0e-4));
            for iter = 1:maxiter
                % Begin the iteration by overwriting S with a vector that has the required length and direction,
                % except that termination occurs if the given D and S are nearly parallel.
                % TOL is the tolerance for telling whether S and D are nearly parallel. In Powell's code, the
                % tolerance is 1.0D-4. We adapt it to the following value in case single precision is in use.

                % Powell's code calculates S as follows. In precise arithmetic, INPROD(S, D) = 0, ||S|| = ||D||.
                % However, when DD*SS - DS**2 is tiny, the error in S can be large and hence damage these
                % equalities significantly. This did happen in tests, especially when using the single precision.
                % %ds = inprod(d, s)
                % %ss = inprod(s, s)
                % %if (dd * ss - ds**2 <= 1.0E-8_RP * dd * ss) then
                % %    exit
                % %end if
                % %denom = sqrt(dd * ss - ds**2)
                % %s = (dd * s - ds * d) / denom

                % We calculate S as follows. It did improve the performance of NEWUOA in our test.
                ss = sum(s .* s, 'all');
                s = s - linalg_obj.project1(s, d); % PROJECT(X, V) is the projection of X to SPAN(V): X'*(V/||V||)*(V/||V||)
                % N.B.:
                % 1. The condition ||S||<=TOL*SQRT(SS) below is equivalent to DS^2>=(1-TOL^2)*DD*SS in theory.
                % As shown above, Powell's code triggers an exit if DS^2>=(1-1.0E-8)*DD*SS. So our condition is
                % the same except that we take EPS into account in case single precision is in use.
                % 2. The condition below should be non-strict so that ||S|| = 0 can trigger the exit.
                if norm(s) <= tol * sqrt(ss)
                    break
                end
                s = (norm(d) / norm(s)) * s;

                % In precise arithmetic, INPROD(S, D) = 0 and ||S|| = ||D|| = DELBAR.
                if abs(sum(d .* s, 'all')) >= 0.1 * norm(d) * norm(s) || norm(s) >= 2.0 * delbar
                    break
                end

                w = powalg_obj.hess_mul(s, xpt, pqlag); % W = MATPROD(XPT, PQLAG * MATPROD(S, XPT))

                % Seek the value of the angle that maximizes ||TAU||.
                % First, calculate the coefficients of the objective function on the circle.
                cf(1) = 0.5 * sum(s .* w, 'all');
                cf(2) = sum(d .* gc, 'all');
                cf(3) = sum(s .* gc, 'all');
                cf(4) = 0.5 * sum(d .* gd, 'all') - cf(1);
                cf(5) = sum(s .* gd, 'all');
                % The 50 in the line below was chosen by Powell. It works the best in tests, MAGICALLY. Larger
                % (e.g., 60, 100) or smaller (e.g., 20, 40) values will worsen the performance of NEWUOA. Why??
                angle = univar_obj.circle_maxabs(@(varargin) obj.circle_fun_biglag(varargin{:}), cf, 50);

                % Calculate the new D and GD.
                cth = cos(angle);
                sth = sin(angle);
                dold(:) = d;
                d = cth * d + sth * s;

                % Exit in case of Inf/NaN in D.
                if ~isfinite(sum(abs(d), 'all'))
                    d(:) = dold;
                    break
                end

                % Test for convergence.
                if abs(obj.circle_fun_biglag(angle, cf)) <= 1.1 * abs(obj.circle_fun_biglag(0.0, cf))
                    break
                end

                % Calculate GD and S.
                gd = cth * gd + sth * w;
                s = gc + gd;
            end

            %====================%
            %  Calculation ends  %
            %====================%



        end
        function d = bigden(obj, idz, knew, kopt, bmat, d0, xpt, zmat)
            %--------------------------------------------------------------------------------------------------%
            % BIGDEN calculates a D by approximately solving
            %
            % max |SIGMA(XOPT + D)|, subject to ||D|| <= DELBAR,
            %
            % where SIGMA is the denominator sigma in the updating formula (4.11)--(4.12) for H, which is the
            % inverse of the coefficient matrix for the interpolation system (see (3.12)). Indeed, each column
            % of H corresponds to a Lagrange basis function of the interpolation problem.  See Section 6 of the
            % NEWUOA paper.
            % N.B.:
            % In Powell's code, BIGDEN calculates also the VLAG and BETA for the selected D. Here, to reduce the
            % coupling of code, we return only D but compute VLAG and BETA outside by calling VLAGBETA. It makes
            % no difference mathematically, but the computed VLAG/BETA will change slightly due to rounding.
            %--------------------------------------------------------------------------------------------------%



            linalg_obj = prima_mat.common.linalg_mod();
            powalg_obj = prima_mat.common.powalg_mod();
            univar_obj = prima_mat.common.univar_mod();


            % BMAT(N, NPT+N)

            % XPT(N, NPT)
            % ZMAT(NPT, NPT - N - 1)


            d = NaN(size(xpt, 1), 1); % D(N)

            % Local variable


            j = NaN;
            k = NaN;


            nw = NaN;

            angle = NaN;


            den = NaN(9, 1);
            denex = NaN(9, 1);
            denmax = NaN;

            dold = NaN(size(xpt, 1), 1);

            dstemp = NaN(size(xpt, 2), 1);

            par = NaN(5, 1);
            pqlag = NaN(size(xpt, 2), 1);
            prod_custom = NaN(size(xpt, 1) + size(xpt, 2), 5);
            s = NaN(size(xpt, 1), 1);

            sstemp = NaN(size(xpt, 2), 1);
            tau = NaN;
            tempa = NaN;
            tempb = NaN;
            tempc = NaN;

            v = NaN(size(xpt, 2), 1);
            vlag = NaN(size(xpt, 1) + size(xpt, 2), 1);
            w = NaN(size(xpt, 1) + size(xpt, 2), 5);
            x = NaN(size(xpt, 1), 1);
            xd = NaN;
            xptemp = NaN(size(xpt, 1), size(xpt, 2));
            xs = NaN;

            y = NaN(size(xpt, 1), 1);
            yd = NaN;
            ysq = NaN;


            n = size(xpt, 1);
            npt = size(xpt, 2);


            %====================%
            % Calculation starts %
            %====================%

            x(:) = xpt(:, kopt); % For simplicity, we use X to denote XOPT.

            delbar = norm(d0); % In theory, ||D0|| = DELBAR.

            % PQLAG contains the leading NPT elements of the KNEW-th column of H, and it provides the second
            % derivative parameters of LFUNC.
            pqlag(:) = powalg_obj.omega_col(idz, zmat, knew);
            alpha = pqlag(knew); % ALPHA is the KNEW-th diagonal entry of H, i.e., that of Omega.

            % The initial search direction D is taken from the last call of BIGLAG, and the initial S is set
            % below, usually to the direction from X to X_KNEW, but a different direction to an interpolation
            % point may be chosen, in order to prevent S from being nearly parallel to D.
            d(:) = d0;
            dd = sum(d .* d, 'all');
            s(:) = xpt(:, knew) - x;
            ds = sum(d .* s, 'all');
            ss = sum(s .* s, 'all');
            xsq = sum(x .* x, 'all');

            if ~(ds ^ 2 <= 0.99 * dd * ss)
                % `.NOT. (A <= B)` differs from `A > B`.  The former holds iff A > B or {A, B} contains NaN.
                dtest = ds ^ 2 / ss;
                xptemp(:, :) = xpt - x;
                %%MATLAB: xptemp = xpt - x  % x should be a column! Implicit expansion
                %----------------------------------------------------------------%
                %---------!dstemp = matprod(d, xpt) - inprod(x, d) !-------------%
                dstemp(:) = xptemp.' * d;
                %----------------------------------------------------------------%
                sstemp(:) = sum((xptemp) .^ 2, 1);

                dstemp(kopt) = 2.0 * ds + 1.0;
                sstemp(kopt) = ss;
                k = fortran.minloc(dstemp .^ 2 ./ sstemp, 'dim', 1);
                % K can be 0 due to NaN. In that case, set K = KNEW. Otherwise, memory errors will occur.
                if k == 0
                    k = knew;
                end
                if (~(dstemp(k) ^ 2 / sstemp(k) >= dtest)) && k ~= kopt
                    % `.NOT. (A >= B)` differs from `A < B`.  The former holds iff A < B or {A, B} contains NaN.
                    % Although unlikely, if NaN occurs, it may happen that K = KOPT.
                    s(:) = xpt(:, k) - x;
                end
            end

            densav = 0.0;

            tol = min(0.1, max(eps(1.0) ^ 0.25, 1.0e-4));
            for iter = 1:n
                % Begin the iteration by overwriting S with a vector that has the required length and direction.
                % TOL is the tolerance for telling whether S and D are nearly parallel. In Powell's code, the
                % tolerance is 1.0D-4. We adapt it to the following value in case single precision is in use.

                % Powell's code calculates S as follows. In precise arithmetic, INPROD(S, D) = 0, ||S|| = ||D||.
                % However, when DD*SS - DS**2 is tiny, the error in S can be large and hence damage these
                % equalities significantly. This did happen in tests, especially when using the single precision.
                % %ds = inprod(d, s)
                % %ss = inprod(s, s)
                % %ssden = dd * ss - ds**2
                % %if (ssden < 1.0E-8_RP * dd * ss) then
                % %    exit
                % %end if
                % %s = (ONE / sqrt(ssden)) * (dd * s - ds * d)

                % We calculate S as below. It did improve the performance of NEWUOA in our test.
                ss = sum(s .* s, 'all');
                s = s - linalg_obj.project1(s, d); % PROJECT(X, V) is the projection of X to SPAN(V): X'*(V/||V||)*(V/||V||)
                % N.B.:
                % 1. The condition ||S||<=TOL*SQRT(SS) below is equivalent to DS^2>=(1-TOL^2)*DD*SS in theory.
                % As shown above, Powell's code triggers an exit if DS^2>=(1-1.0E-8)*DD*SS. So our condition is
                % the same except that we take EPS into account in case single precision is in use.
                % 2. The condition below should be non-strict so that ||S|| = 0 can trigger the exit.
                if norm(s) <= tol * sqrt(ss)
                    break
                end
                s = (s ./ norm(s)) * norm(d);
                % In precise arithmetic, INPROD(S, D) = 0 and ||S|| = ||D|| = DELBAR = ||D0||.
                if abs(sum(d .* s, 'all')) >= 0.1 * norm(d) * norm(s) || norm(s) >= 2.0 * delbar
                    break
                end

                % Set the coefficients of the first two terms of BETA.
                xd = sum(x .* d, 'all');
                xs = sum(x .* s, 'all');
                dd = sum(d .* d, 'all');
                tempa = 0.5 * xd * xd;
                tempb = 0.5 * xs * xs;
                den(1) = dd * (xsq + 0.5 * dd) + tempa + tempb;
                den(2) = 2.0 * xd * dd;
                den(3) = 2.0 * xs * dd;
                den(4) = tempa - tempb;
                den(5) = xd * xs;
                den(6:9) = 0.0;

                % Put the coefficients of WCHECK in W.
                for k = 1:npt
                    tempa = sum(xpt(:, k) .* d, 'all');
                    tempb = sum(xpt(:, k) .* s, 'all');
                    tempc = sum(xpt(:, k) .* x, 'all');
                    w(k, 1) = 0.25 * (tempa ^ 2 + tempb ^ 2);
                    w(k, 2) = tempa * tempc;
                    w(k, 3) = tempb * tempc;
                    w(k, 4) = 0.25 * (tempa ^ 2 - tempb ^ 2);
                    w(k, 5) = 0.5 * tempa * tempb;
                end
                w(npt + 1:npt + n, 1:5) = 0.0;
                w(npt + 1:npt + n, 2) = d;
                w(npt + 1:npt + n, 3) = s;

                % Put the coefficients of THETA*WCHECK in PROD.
                for j = 1:5
                    prod_custom(1:npt, j) = powalg_obj.omega_mul(idz, zmat, w(1:npt, j));
                    nw = npt;
                    if j == 2 || j == 3
                        prod_custom(1:npt, j) = prod_custom(1:npt, j) + bmat(:, 1:npt).' * w(npt + 1:npt + n, j);
                        nw = npt + n;
                    end
                    prod_custom(npt + 1:npt + n, j) = bmat(:, 1:nw) * w(1:nw, j);
                end

                % Include in DEN the part of BETA that depends on THETA.
                for k = 1:npt + n
                    par(1:5) = 0.5 * prod_custom(k, 1:5) .* w(k, 1:5);
                    den(1) = den(1) - par(1) - sum(par(1:5), 'all');
                    tempa = prod_custom(k, 1) * w(k, 2) + prod_custom(k, 2) * w(k, 1);
                    tempb = prod_custom(k, 2) * w(k, 4) + prod_custom(k, 4) * w(k, 2);
                    tempc = prod_custom(k, 3) * w(k, 5) + prod_custom(k, 5) * w(k, 3);
                    den(2) = den(2) - tempa - 0.5 * (tempb + tempc);
                    den(6) = den(6) - 0.5 * (tempb - tempc);
                    tempa = prod_custom(k, 1) * w(k, 3) + prod_custom(k, 3) * w(k, 1);
                    tempb = prod_custom(k, 2) * w(k, 5) + prod_custom(k, 5) * w(k, 2);
                    tempc = prod_custom(k, 3) * w(k, 4) + prod_custom(k, 4) * w(k, 3);
                    den(3) = den(3) - tempa - 0.5 * (tempb - tempc);
                    den(7) = den(7) - 0.5 * (tempb + tempc);
                    tempa = prod_custom(k, 1) * w(k, 4) + prod_custom(k, 4) * w(k, 1);
                    den(4) = den(4) - tempa - par(2) + par(3);
                    tempa = prod_custom(k, 1) * w(k, 5) + prod_custom(k, 5) * w(k, 1);
                    tempb = prod_custom(k, 2) * w(k, 3) + prod_custom(k, 3) * w(k, 2);
                    den(5) = den(5) - tempa - 0.5 * tempb;
                    den(8) = den(8) - par(4) + par(5);
                    tempa = prod_custom(k, 4) * w(k, 5) + prod_custom(k, 5) * w(k, 4);
                    den(9) = den(9) - 0.5 * tempa;
                end

                par(1:5) = 0.5 * prod_custom(knew, 1:5) .^ 2;
                denex(1) = alpha * den(1) + par(1) + sum(par(1:5), 'all');
                tempa = 2.0 * prod_custom(knew, 1) * prod_custom(knew, 2);
                tempb = prod_custom(knew, 2) * prod_custom(knew, 4);
                tempc = prod_custom(knew, 3) * prod_custom(knew, 5);
                denex(2) = alpha * den(2) + tempa + tempb + tempc;
                denex(6) = alpha * den(6) + tempb - tempc;
                tempa = 2.0 * prod_custom(knew, 1) * prod_custom(knew, 3);
                tempb = prod_custom(knew, 2) * prod_custom(knew, 5);
                tempc = prod_custom(knew, 3) * prod_custom(knew, 4);
                denex(3) = alpha * den(3) + tempa + tempb - tempc;
                denex(7) = alpha * den(7) + tempb + tempc;
                tempa = 2.0 * prod_custom(knew, 1) * prod_custom(knew, 4);
                denex(4) = alpha * den(4) + tempa + par(2) - par(3);
                tempa = 2.0 * prod_custom(knew, 1) * prod_custom(knew, 5);
                denex(5) = alpha * den(5) + tempa + prod_custom(knew, 2) * prod_custom(knew, 3);
                denex(8) = alpha * den(8) + par(4) - par(5);
                denex(9) = alpha * den(9) + prod_custom(knew, 4) * prod_custom(knew, 5);

                % Seek the value of the angle that maximizes the |DENOM|.
                angle = univar_obj.circle_maxabs(@(varargin) obj.circle_fun_bigden(varargin{:}), denex, 50);

                % Calculate the new D.
                dold = d;
                d = cos(angle) * d + sin(angle) * s;

                % Exit in case of Inf/NaN in D.
                if ~isfinite(sum(abs(d), 'all'))
                    d = dold;
                    break
                end

                % Test for convergence.
                if iter > 1
                    densav = max(densav, obj.circle_fun_bigden(0.0, denex));
                end
                denmax = obj.circle_fun_bigden(angle, denex);
                if abs(denmax) <= 1.1 * abs(densav)
                    break
                end
                densav = denmax;

                % Set S to HALF the gradient of the denominator with respect to D. First, calculate the new VLAG.
                par(:) = [1.0, cos(angle), sin(angle), cos(2.0 * angle), sin(2.0 * angle)];
                vlag(:) = prod_custom * par;
                tau = vlag(knew);
                y = x + d;
                yd = sum(y .* d, 'all');
                ysq = sum(y .* y, 'all');
                v = (tau * pqlag - alpha * vlag(1:npt)) .* (xpt.' * y);
                s(:) = tau * bmat(:, knew) + alpha * (yd * x + ysq * d - vlag(npt + 1:npt + n));
                s = s + xpt * v;
            end

            %====================%
            %  Calculation ends  %
            %====================%



        end
        function f = circle_fun_biglag(~, theta, args)
            %--------------------------------------------------------------------------------------------------%
            % This function defines the objective function of the 2-dimensional search on a circle in BIGLAG.
            %--------------------------------------------------------------------------------------------------%



            f = NaN;


            %====================%
            % Calculation starts %
            %====================%

            cth = cos(theta);
            sth = sin(theta);
            f = args(1) + (args(2) + args(4) * cth) * cth + (args(3) + args(5) * cth) * sth;

            %====================%
            %  Calculation ends  %
            %====================%
        end
        function f = circle_fun_bigden(~, theta, args)
            %--------------------------------------------------------------------------------------------------%
            % This function defines the objective function of the 2-dimensional search on a circle in BIGDEN.
            %--------------------------------------------------------------------------------------------------%



            f = NaN;


            par = NaN(numel(args), 1);


            %====================%
            % Calculation starts %
            %====================%

            par(1) = 1.0;
            par(2:2:8) = cos(theta * [1.0, 2.0, 3.0, 4.0]);
            par(3:2:9) = sin(theta * [1.0, 2.0, 3.0, 4.0]);
            f = sum(args .* par, 'all');

            %====================%
            %  Calculation ends  %
            %====================%
        end

    end
end