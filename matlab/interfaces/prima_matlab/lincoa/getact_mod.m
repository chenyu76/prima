classdef getact_mod
    %--------------------------------------------------------------------------------------------------%
    % This module provides the GETACT subroutine of LINCOA. It is used only in the trust-region
    % subproblem solver. We do not put it in trustregion.f90 as it is too long.
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
    % Last Modified: Saturday, March 09, 2024 PM12:12:21
    %--------------------------------------------------------------------------------------------------%

    methods
        function [iact, nact, qfac, resact, resnew, rfac, psd] = getact(obj, amat, delta, g, iact, nact, qfac, resact, resnew, rfac, psd)
            %--------------------------------------------------------------------------------------------------%
            %-------------------------------------------------------------%
            % THE FOLLOWING DESCRIPTION NEEDS VERIFICATION!               %
            % Note that the set JJ gets updated within this subroutine,   %
            % which seems inconsistent with the description below.        %
            % See the lines below "Pick the next integer L or terminate". %
            %-------------------------------------------------------------%
            %
            % This subroutine solves a linearly constrained projected problem (LCPP)
            %
            % min ||D + G|| subject to AMAT(:, j)^T * D <= 0 for j in JJ.
            %
            % The solution is PSD, which is a projected steepest descent direction PSD for a linearly
            % constrained trust-region subproblem (LCTRS)
            %
            % min Q(X_k + D)  subject to ||D|| <= Delta and AMAT^T*(X_k + D) <= B,
            %
            % where X_k is in R^N, B is in R^M, and AMAT is in R^{NxM}.
            %
            % In (LCPP), JJ is the index set defined in (3.3) of Powell (2015) as
            %
            % JJ = {j : B_j - A_j^T*Y <= 0.2*Delta*||A_j||, 1 <= j <= M} with A_j = AMAT(:, j),
            %
            % i.e., the index set of the nearly active constraints of (LCTRS) (Powell wrote that j is in JJ if
            % and only if the distance from Y to the boundary of the j-th constraint is at most 0.2*Delta).
            % Here, Y is the point where G is taken, namely G = nabla Q(Y). Y is not necessarily X_k, but an
            % iterate of the algorithm (e.g., truncated conjugate gradient) that solves (LCTRS). In LINCOA,
            % ||A_j|| is 1 as the gradients of the linear constraints are normalized before LINCOA starts.
            %
            % The subroutine solves (LCPP) by the active set method of Goldfarb-Idnani (1983). It does not only
            % calculate PSD, but also identify the active set of (LCPP) at the solution PSD, namely
            %
            % II = {j in JJ : AMAT(:, j)^T*PSD = 0} (see (3.5) of Powell (2015)),
            %
            % and maintains a QR factorization of A corresponding to the active set. More specifically,
            % IACT(1:NACT) is a set of indices such that the columns of AMAT(:, IACT(1:NACT)) constitute a basis
            % of the "active constraint" gradients (i.e., those corresponding to the set II mentioned above, but
            % not JJ!), and QFAC*RFAC(:, 1:NACT) is the QR factorization of! AMAT(:, IACT(1:NACT)) such that
            %
            % SIZE(QFAC) = [N, N], SIZE(RFAC, 1) = N, diag(RFAC(:, 1:NACT)) > 0.
            %
            % NACT, IACT, QFAC and RFAC are maintained up to date across invocations of GETACT for warm starts.
            %
            % DELTA, RESNEW, RESACT, and G are the same as the terms with these names in SUBROUTINE TRSTEP.
            % The elements of RESNEW and RESACT are also kept up to date. See Section 3 of Powell (2015).
            % Note that the updates only permute RESACT but do not change the values inside.
            %
            % VLAM is the vector of Lagrange multipliers of the calculation.
            %
            % See Section 3 of Powell (2015) for more information.
            %--------------------------------------------------------------------------------------------------%

            % Common modules
            consts_obj = consts_mod();
            debug_obj = debug_mod();
            infnan_obj = infnan_mod();
            linalg_obj = linalg_mod();


            % Inputs
            % AMAT(N, M)

            % G(N)

            % In-outputs
            % IACT(M)

            % QFAC(N, N)
            % RESACT(M)
            % RESNEW(M)
            % RFAC(N, N)

            % Outputs
            % PSD(N)

            % Local variables
            srname = "GETACT";
            icon = NaN;
            iter = NaN;
            l = NaN;
            m = NaN;
            maxiter = NaN;
            n = NaN;
            mask = false(size(amat, 2), 1);
            apsd = NaN(size(amat, 2), 1);
            dd = NaN;
            ddsav = NaN;
            dnorm = NaN;
            gg = NaN;
            frac = NaN(numel(g), 1);
            psdsav = NaN(numel(psd), 1);
            tdel = NaN;
            tol = NaN;
            v = NaN(numel(g), 1);
            violmx = NaN;
            vlam = NaN(numel(g), 1);
            vmu = NaN(numel(g), 1);
            vmult = NaN;

            % Sizes.
            m = size(amat, 2);
            n = fix(numel(g));

            % Preconditions
            if consts_obj.DEBUGGING
                debug_obj.assert(m >= 0, "M >= 0", srname);
                debug_obj.assert(n >= 1, "N >= 1", srname);
                debug_obj.assert(size(amat, 1) == n && size(amat, 2) == m, "SIZE(AMAT) == [N, M]", srname);
                debug_obj.assert(all(infnan_obj.is_finite(g), 'all'), "G is finite", srname);
                debug_obj.assert(nact >= 0 && nact <= min(m, n), "0 <= NACT <= MIN(M, N)", srname);
                debug_obj.assert(numel(iact) == m, "SIZE(IACT) == M", srname);
                debug_obj.assert(all(iact(1:nact) >= 1 & iact(1:nact) <= m, 'all'), "1 <= IACT <= M", srname);
                debug_obj.assert(numel(resact) == m, "SIZE(RESACT) == M", srname);
                debug_obj.assert(numel(resnew) == m, "SIZE(RESNEW) == M", srname);
                debug_obj.assert(size(qfac, 1) == n && size(qfac, 2) == n, "SIZE(QFAC) == [N, N]", srname);
                tol = max(consts_obj.TEN ^ max(-10, -consts_obj.MAXPOW10), min(0.1, consts_obj.TEN ^ min(8, consts_obj.MAXPOW10) * consts_obj.EPS * double(n)));
                debug_obj.assert(linalg_obj.isorth(qfac, 'tol', tol), "QFAC is orthogonal", srname);
                debug_obj.assert(size(rfac, 1) == n && size(rfac, 2) == n, "SIZE(RFAC) == [N, N]", srname);
                debug_obj.assert(linalg_obj.istriu(rfac), "RFAC is upper triangular", srname);
                debug_obj.assert(numel(psd) == n, "SIZE(PSD) == N", srname);
            end

            %====================%
            % Calculation starts %
            %====================%

            % Quick return when M = 0.
            if m <= 0
                nact = 0;
                qfac(:, :) = linalg_obj.eye1(n);
                psd(:) = -g;
                return
            end

            % Set some constants.
            gg = linalg_obj.inprod(g, g);
            tdel = 0.2 * delta; % Changing TDEL to 0.1_RP*DELTA does not improve the performance of LINCOA.

            % Set the initial QFAC to the identity matrix in the case NACT = 0.
            if nact == 0
                qfac(:, :) = linalg_obj.eye1(n);
            end

            % Remove any constraints from the initial active set whose residuals exceed TDEL.
            % Compilers may complain if VLAM is not set. The value does not matter, as it will be overwritten.
            vlam = repmat(consts_obj.ZERO, size(vlam));
            for icon = nact:-1:1
                if resact(icon) > tdel
                    % Delete constraint IACT(ICON) from the active set, and set NACT = NACT - 1.
                    [iact, nact, qfac, resact, resnew, rfac, vlam] = obj.delact(icon, iact, nact, qfac, resact, resnew, rfac, vlam);
                end
            end

            % Remove any constraints from the initial active set whose Lagrange multipliers are nonnegative,
            % and set the surviving multipliers.
            % The following loop will run for at most NACT times, since each call of DELACT reduces NACT by 1.
            while nact > 0
                vlam(1:nact) = linalg_obj.lsqr_Rfull(g, qfac(:, 1:nact), rfac(1:nact, 1:nact));
                if ~any(vlam(1:nact) >= 0, 'all')
                    break
                end
                icon = max(linalg_obj.trueloc(vlam(1:nact) >= 0), [], 'all');
                %%MATLAB: icon = max(find(vlam(1:nact) >= 0)); % OR: icon = find(vlam(1:nact) >= 0, 1, 'last')
                [iact, nact, qfac, resact, resnew, rfac, vlam] = obj.delact(icon, iact, nact, qfac, resact, resnew, rfac, vlam);
            end
            % Zaikun 20220330: What if NACT = 0 at this point?

            % Set the new search direction D. Terminate if the 2-norm of D is ZERO or does not decrease, or if
            % NACT=N holds. The situation NACT=N occurs for sufficiently large DELTA if the origin is in the
            % convex hull of the constraint gradients.
            % Start with initialization of PSDSAV and DDSAV.
            psdsav = repmat(consts_obj.ZERO, size(psdsav)); % Must be set, in case the loop exits due to abnormality at iteration 1.
            ddsav = consts_obj.TWO * gg; % By Powell. This value is used at iteration 1 to test whether DD >= DDSAV. Why?

            % What is the theoretical maximal number of iterations in the following procedure? Powell's code for
            % this part is essentially a `DO WHILE (NACT < N) ... END DO` loop. We enforce the following maximal
            % number of iterations, which is never reached in our tests (indeed, even 2*N cannot be reached).
            % N.B.: 1. The formulation of MAXITER below contains a precaution against overflow. In
            % MATLAB/Python/Julia/R, we can write maxiter = min(10000, 2*(m + n))
            % 2. The iteration counter ITER never appears in the code of the iterations, as its purpose is
            % merely to impose an upper bound on the number of iterations.
            maxiter = fix(min(10 ^ min(4, floor(log10(double(intmax('int64'))))), 2 * fix(m + n)));
            for iter = 1:maxiter
                % When NACT == N, exit with PSD = 0. Indeed, with a correctly implemented matrix product, the
                % lines below this IF should render DD = 0 and trigger an exit. We make it explicit for clarity.
                if nact >= n
                    % Indeed, NACT > N should never happen.
                    psd = repmat(consts_obj.ZERO, size(psd));
                    break
                end

                % Set PSD to the projection of -G to range(QFAC(:,NACT+1:N))
                psd(:) = -linalg_obj.matprod21(qfac(:, nact + 1:n), linalg_obj.matprod12(g, qfac(:, nact + 1:n)));
                %%MATLAB: psd = -qfac(:, nact + 1:n) * (g' * qfac(:, nact + 1:n))';
                %----------------------------------------------------------------------------------------------%
                % Zaikun: The schemes below work evidently worse than the one above in a test on 20220417. Why?
                %-------------------------------------------------------------------------%
                % VERSION 1:
                % %psd = matprod(qfac(:, 1:nact), matprod(g, qfac(:, 1:nact))) - g
                %-------------------------------------------------------------------------%
                % VERSION 2:
                % %if (2 * nact < n) then
                % %    psd = matprod(qfac(:, 1:nact), matprod(g, qfac(:, 1:nact))) - g
                % %else
                % %    psd = -matprod(qfac(:, nact + 1:n), matprod(g, qfac(:, nact + 1:n)))
                % %end if
                %-------------------------------------------------------------------------%
                %----------------------------------------------------------------------------------------------%

                dd = linalg_obj.inprod(psd, psd);
                dnorm = sqrt(dd);

                if dnorm <= consts_obj.EPS || infnan_obj.is_nan_sp(dnorm)
                    break
                end

                if dd >= ddsav
                    psd = repmat(consts_obj.ZERO, size(psd)); % Zaikun 20220329: Powell wrote this. Why?
                    %psd = psdsav  ! This does not seem to improve the performance.
                    break
                end

                %---------------------------------------------------------------------------------------%
                % Powell's code does not handle the following pathological cases.
                if linalg_obj.inprod(psd, g) > 0 || ~infnan_obj.is_finite(sum(abs(psd), 'all'))
                    psd(:) = psdsav;
                    break
                end
                % In our tests, tolerating the following cases seems to render better numerical results.
                % %if (dd > gg) then
                % %    psd = (sqrt(gg) / dnorm) * psd
                % %    exit
                % %end if
                % %if (inprod(psd, g) < -gg) then
                % %    exit
                % %end if
                %---------------------------------------------------------------------------------------%

                psdsav(:) = psd;
                ddsav = dd;

                % Pick the next integer L or terminate; a positive L is the index of the most violated constraint.
                apsd(:) = linalg_obj.matprod12(psd, amat);
                mask(:) = (resnew > 0 & resnew <= tdel & apsd > (dnorm / delta) * resnew);
                %----------------------------------------------------------------------------------------------%
                % N.B.: the definition of L and VIOLMX can be simplified as follows, but we prefer explicitness.
                %L = INT(MAXLOC(APSD, MASK=MASK, DIM=1), IK) ! MAXLOC(...) = 0 if MASK is all FALSE.
                %VIOLMX = MAXVAL(APSD, MASK=MASK)  ! MAXVAL(...) = -HUGE(APSD) if MASK is all FALSE.
                if any(mask, 'all')
                    l = fix(fortran.maxloc(apsd, 'mask', mask, 'dim', 1));
                    violmx = apsd(l);
                else
                    l = 0;
                    violmx = -consts_obj.REALMAX;
                end
                %%MATLAB: apsd(mask) = -Inf; [violmx, l] = max(apsd);
                % N.B.: the value of L will differ from the Fortran version if MASK is all FALSE, but this is OK
                % because VIOLMX will be -Inf, which will trigger the `exit` below. This is tricky. Be cautious!
                %----------------------------------------------------------------------------------------------%

                % Terminate if VIOLMX <= 0 (when MASK contains only FALSE) or a positive value of VIOLMX may be
                % due to computer rounding errors.
                % N.B.: 1. Theoretically (but not numerically), APSD(IACT(1:NACT)) = 0 or empty.
                % 2. CAUTION: the Inf-norm of APSD(IACT(1:NACT)) is NOT always MAXVAL(ABS(APSD(IACT(1:NACT)))),
                % as the latter returns -HUGE(APSD) instead of 0 when NACT = 0! In MATLAB, max([]) = []; in
                % Python, R, and Julia, the maximum of an empty array raises errors/warnings (as of 20220318).
                % Powell's condition for the IF is as follows. Very often, the threshold is almost zero.
                % %if (all(.not. mask) .or. violmx <= min(0.01_RP * dnorm, TEN * norm(apsd(iact(1:nact)), 'inf'))) then
                % The following condition works essentially the same as Powell's. However, it ensures that
                % VIOLMX > EPS * DNORM when the EXIT is not triggered, which implies that AMAT(:, L) is not in
                % the range of QFAC(:, 1:NACT).
                if all(~mask, 'all') || violmx <= max(consts_obj.EPS * dnorm, consts_obj.TEN * linalg_obj.named_norm_vec(apsd(iact(1:nact)), "inf"))
                    break
                end

                % Add constraint L to the active set. ADDACT sets NACT = NACT + 1 and VLAM(NACT) = 0.
                [iact, nact, qfac, resact, resnew, rfac, vlam] = obj.addact(l, amat(:, l), iact, nact, qfac, resact, resnew, rfac, vlam);

                % Set the components of the vector VMU if VIOLMX is positive.
                % N.B.: 1. In theory, NACT > 0 is not needed in the condition below, because VIOLMX must be 0
                % when NACT is 0. We keep NACT > 0 for security: when NACT <= 0, RFAC(NACT, NACT) is invalid.
                % 2. The loop will run for at most NACT <= N times: if VIOLMX > 0, then ICON > 0, and hence
                % VLAM(ICON) = 0, which implies that DELACT will be called to reduce NACT by 1.
                while violmx > 0 && nact > 0
                    v(1:nact - 1) = consts_obj.ZERO;
                    v(nact) = consts_obj.ONE / rfac(nact, nact); % This is why we must ensure NACT > 0.
                    % Solve the linear system RFAC(1:NACT, 1:NACT) * VMU(1:NACT) = V(1:NACT) .
                    vmu(1:nact) = linalg_obj.solve(rfac(1:nact, 1:nact), v(1:nact)); % VMU(NACT) = V(NACT)/RFAC(NACT,NACT)>0
                    %%MATLAB: vmu(1:nact) = rfac(1:nact, 1:nact) \ v(1:nact);

                    % Calculate the multiple of VMU to subtract from VLAM, and update VLAM.
                    % N.B.: 1. VLAM(1:NACT-1) < 0 and VLAM(NACT) <= 0 by the updates of VLAM. 2. VMU(NACT) > 0.
                    % 3. Only the places where VMU(1:NACT) < 0 is relevant below, if any.
                    frac = repmat(consts_obj.REALMAX, size(frac));
                    frac(vmu(1:nact) < 0 & vlam(1:nact) < 0) = vlam(vmu(1:nact) < 0 & vlam(1:nact) < 0) ./ vmu(vmu(1:nact) < 0 & vlam(1:nact) < 0);
                    %%MATLAB: frac = vlam / vmu; frac(vmu >= 0 | vlam >= 0) = Inf;
                    vmult = min([violmx; reshape(frac(1:nact), [], 1)], [], 'all');
                    icon = max([0; reshape(linalg_obj.trueloc(frac(1:nact) <= vmult), [], 1)], [], 'all');
                    %%MATLAB: icon = max([0; find(frac(1:nact) <= vmult)]); % find(frac(1:nact)<=vmult) can be empty

                    % N.B.: 0. The definition of ICON given above is mathematically equivalent to the following.
                    % %ICON = MAXVAL(TRUELOC([VIOLMX, FRACMULT(1:NACT)] <= VMULT)) - 1_IK, OR
                    % %ICON = INT(MINLOC([VIOLMX, FRACMULT(1:NACT)], DIM=1, BACK=.TRUE.), IK) - 1_IK
                    % However, such implementations are problematic in the unlikely case of VMULT = NaN: ICON
                    % will be -Inf in the first and unspecified in the second. The MATLAB counterpart of the
                    % first implementation will render ICON = [] as `find` (the MATLAB version of TRUELOC)
                    % returns [].
                    % 1. The BACK argument in MINLOC is available in F2008. Not supported by Absoft as of 2022.
                    % 2. A motivation for backward MINLOC is to save computation in DELACT below (what else?).

                    violmx = max(violmx - vmult, consts_obj.ZERO);
                    vlam(1:nact) = vlam(1:nact) - vmult * vmu(1:nact);
                    if icon > 0 && icon <= nact
                        % Powell: IF (ICON>0). We check ICON<=NACT for safety.
                        vlam(icon) = consts_obj.ZERO;
                    end

                    % Reduce the active set if necessary, so that all components of the new VLAM are negative,
                    % with resetting of the residuals of the constraints that become inactive.
                    for icon = nact:-1:1
                        if vlam(icon) >= 0
                            % Powell's version: IF (.NOT. VLAM(ICON) < 0) THEN
                            % Delete the constraint with index IACT(ICON) from the active set; set NACT = NACT-1.
                            [iact, nact, qfac, resact, resnew, rfac, vlam] = obj.delact(icon, iact, nact, qfac, resact, resnew, rfac, vlam);
                        end
                    end
                end % End of DO WHILE (VIOLMX > 0 .AND. NACT > 0)

                %----------------------------------------------------------------------------------------------%
                % NACT can become 0 at this point iff VLAM(1:NACT) >= 0 before calling DELACT, which is true
                % if NACT happens to be 1 when the WHILE loop starts. However, we have never observed a failure
                % of the assertion below as of 20220329. Why?
                %-----------------------------------------%
                debug_obj.assert(nact > 0, "NACT > 0", srname); %
                %-----------------------------------------%
                if nact == 0
                    break
                end
                %----------------------------------------------------------------------------------------------%
            end % End of DO WHILE (NACT < N)

            % It is possible to have NACT == 0 here. The following lines improve the performance of LINCOA.
            % Powell's code does not take care of this case explicitly.
            if nact == 0
                qfac(:, :) = linalg_obj.eye1(n);
                psd(:) = -g;
            end

            %====================%
            %  Calculation ends  %
            %====================%

            % Postconditions
            if consts_obj.DEBUGGING
                % During the development, we want to get alerted if ITER reaches MAXITER.
                debug_obj.assert(iter < maxiter, "ITER < MAXITER", srname);
                debug_obj.assert(nact >= 0 && nact <= min(m, n), "0 <= NACT <= MIN(M, N)", srname); % Can NACT be 0?
                debug_obj.assert(numel(iact) == m, "SIZE(IACT) == M", srname);
                debug_obj.assert(all(iact(1:nact) >= 1 & iact(1:nact) <= m, 'all'), "1 <= IACT <= M", srname);
                debug_obj.assert(size(qfac, 1) == n && size(qfac, 2) == n, "SIZE(QFAC) == [N, N]", srname);
                debug_obj.assert(linalg_obj.isorth(qfac, 'tol', tol), "QFAC is orthogonal", srname);
                debug_obj.assert(size(rfac, 1) == n && size(rfac, 2) == n, "SIZE(RFAC) == [N, N]", srname);
                debug_obj.assert(linalg_obj.istriu(rfac), "RFAC is upper triangular", srname);
                debug_obj.assert(numel(psd) == n, "SIZE(PSD) == N", srname);
                % PSD = -G when NACT == 0; G may contain Inf/NaN.
                debug_obj.assert(all(infnan_obj.is_finite(psd), 'all') || nact == 0, "PSD is finite unless NACT == 0", srname);
                % In theory, ||PSD||^2 <= GG and -GG <= PSD^T*G <= 0.
                % N.B. 1. Do not use DD, which may not be up to date. 2. PSD^T*G can be NaN if G is huge.
                debug_obj.assert(linalg_obj.inprod(psd, psd) <= consts_obj.TWO * gg, "||PSD||^2 <= 2*GG", srname);
                debug_obj.assert(~(linalg_obj.inprod(psd, g) > 100.0 * consts_obj.EPS * gg || linalg_obj.inprod(psd, g) < -consts_obj.TWO * gg), "-2*GG <= PSD^T*G <= 0", srname);
            end

        end
        function [iact, nact, qfac, resact, resnew, rfac, vlam] = addact(~, l, c, iact, nact, qfac, resact, resnew, rfac, vlam)
            %--------------------------------------------------------------------------------------------------%
            % This subroutine adds the constraint with index L to the active set as the (NACT+ )-th active
            % constraint, updates IACT, QFAC, etc accordingly, and increments NACT to NACT+1. Here, C is the
            % gradient of the new active constraint.
            %--------------------------------------------------------------------------------------------------%

            consts_obj = consts_mod();
            debug_obj = debug_mod();
            linalg_obj = linalg_mod();
            powalg_obj = powalg_mod();

            % Inputs

            % C(N)

            % In-outputs
            % IACT(M)

            % QFAC(N, N)
            % RESACT(M)
            % RESNEW(M)
            % RFAC(N, N)
            % VLAM(N)

            % Local variables (debugging only)
            srname = "ADD_ACT";
            m = NaN;
            n = NaN;
            nsave = NaN;
            tol = NaN;

            % Sizes
            m = fix(numel(iact));
            n = fix(numel(vlam));

            % Preconditions
            if consts_obj.DEBUGGING
                debug_obj.assert(m >= 1, "M >= 1", srname); % Should not be called when M == 0.
                debug_obj.assert(n >= 1, "N >= 1", srname);
                debug_obj.assert(nact >= 0 && nact <= min(m, n) - 1, "0 <= NACT <= MIN(M, N)-1", srname);
                debug_obj.assert(l >= 1 && l <= m, "1 <= L <= M", srname);
                debug_obj.assert(all(iact(1:nact) >= 1 & iact(1:nact) <= m, 'all'), "1 <= IACT <= M", srname);
                debug_obj.assert(~any(iact(1:nact) == l, 'all'), "L is not in IACT(1:NACT)", srname);
                debug_obj.assert(size(qfac, 1) == n && size(qfac, 2) == n, "SIZE(QFAC) == [N, N]", srname);
                tol = max(consts_obj.TEN ^ max(-10, -consts_obj.MAXPOW10), min(0.1, consts_obj.TEN ^ min(8, consts_obj.MAXPOW10) * consts_obj.EPS * double(n)));
                debug_obj.assert(linalg_obj.isorth(qfac, 'tol', tol), "QFAC is orthogonal", srname);
                debug_obj.assert(size(rfac, 1) == n && size(rfac, 2) == n, "SIZE(RFAC) == [N, N]", srname);
                debug_obj.assert(linalg_obj.istriu(rfac), "RFAC is upper triangular", srname);
                debug_obj.assert(numel(resact) == m, "SIZE(RESACT) == M", srname);
                debug_obj.assert(numel(resnew) == m, "SIZE(RESNEW) == M", srname);
                nsave = nact; % For debugging only

            end

            %====================%
            % Calculation starts %
            %====================%

            % QRADD applies Givens rotations to the last (N-NACT) columns of QFAC so that the first (NACT+1)
            % columns of QFAC are the ones required for the addition of the L-th constraint, and add the
            % appropriate column to RFAC.
            % N.B.: QRADD always augment NACT by 1, which differs from the corresponding subroutine in COBYLA.
            % It is ensured that C cannot be represented by the gradients of the existing active constraints.
            [qfac, rfac, nact] = powalg_obj.qradd_Rfull(c, qfac, rfac, nact); % NACT is increased by 1!
            % Indeed, it suffices to pass RFAC(:, 1:NACT+1) to QRADD as follows.
            % %call qradd(c, qfac, rfac(:, 1:nact + 1), nact)  ! NACT is increased by 1!

            % Update IACT, RESACT, RESNEW, and VLAM. N.B.: NACT has been increased by 1 in QRADD.
            iact(nact) = l;
            resact(nact) = resnew(l); % RESACT(NACT) = RESNEW(IACT(NACT))
            resnew(l) = consts_obj.ZERO; % RESNEW(IACT(NACT)) = ZERO  ! Why not TINYCV? See DECACT.
            vlam(nact) = consts_obj.ZERO;

            %====================%
            %  Calculation ends  %
            %====================%

            % Postconditions
            if consts_obj.DEBUGGING
                debug_obj.assert(nact == nsave + 1, "NACT = NSAVE + 1", srname);
                debug_obj.assert(nact >= 1 && nact <= min(m, n), "1 <= NACT <= MIN(M, N)", srname);
                debug_obj.assert(all(iact(1:nact) >= 1 & iact(1:nact) <= m, 'all'), "1 <= IACT <= M", srname);
                debug_obj.assert(size(qfac, 1) == n && size(qfac, 2) == n, "SIZE(QFAC) == [N, N]", srname);
                debug_obj.assert(linalg_obj.isorth(qfac, 'tol', tol), "QFAC is orthogonal", srname);
                debug_obj.assert(size(rfac, 1) == n && size(rfac, 2) == n, "SIZE(RFAC) == [N, N]", srname);
                debug_obj.assert(linalg_obj.istriu(rfac), "RFAC is upper triangular", srname);
                debug_obj.assert(numel(resact) == m, "SIZE(RESACT) == M", srname);
                debug_obj.assert(numel(resnew) == m, "SIZE(RESNEW) == M", srname);
            end

        end
        function [iact, nact, qfac, resact, resnew, rfac, vlam] = delact(~, icon, iact, nact, qfac, resact, resnew, rfac, vlam)
            %--------------------------------------------------------------------------------------------------%
            % This subroutine deletes the constraint with index IACT(ICON) from the active set, updates IACT,
            % QFAC, etc accordingly, and reduces NACT to NACT-1.
            %--------------------------------------------------------------------------------------------------%

            consts_obj = consts_mod();
            debug_obj = debug_mod();
            linalg_obj = linalg_mod();
            powalg_obj = powalg_mod();

            % Inputs


            % In-outputs
            % IACT(M)

            % QFAC(N, N)
            % RESACT(M)
            % RESNEW(M)
            % RFAC(N, N)
            % VLAM(N)

            % Local variables (debugging only)
            srname = "DELACT";
            l = NaN;
            m = NaN;
            n = NaN;
            nsave = NaN;
            tol = NaN;

            % Sizes
            m = fix(numel(iact));
            n = fix(numel(vlam));

            % Preconditions
            % Preconditions
            if consts_obj.DEBUGGING
                debug_obj.assert(m >= 1, "M >= 1", srname); % Should not be called when M == 0.
                debug_obj.assert(n >= 1, "N >= 1", srname);
                debug_obj.assert(nact >= 1 && nact <= min(m, n), "1 <= NACT <= MIN(M, N)", srname);
                debug_obj.assert(icon >= 1 && icon <= nact, "1 <= ICON <= NACT", srname);
                debug_obj.assert(all(iact(1:nact) >= 1 & iact(1:nact) <= m, 'all'), "1 <= IACT <= M", srname);
                debug_obj.assert(size(qfac, 1) == n && size(qfac, 2) == n, "SIZE(QFAC) == [N, N]", srname);
                tol = max(consts_obj.TEN ^ max(-10, -consts_obj.MAXPOW10), min(0.1, consts_obj.TEN ^ min(8, consts_obj.MAXPOW10) * consts_obj.EPS * double(n)));
                debug_obj.assert(linalg_obj.isorth(qfac, 'tol', tol), "QFAC is orthogonal", srname);
                debug_obj.assert(size(rfac, 1) == n && size(rfac, 2) == n, "SIZE(RFAC) == [N, N]", srname);
                debug_obj.assert(linalg_obj.istriu(rfac), "RFAC is upper triangular", srname);
                debug_obj.assert(numel(resact) == m, "SIZE(RESACT) == M", srname);
                debug_obj.assert(numel(resnew) == m, "SIZE(RESNEW) == M", srname);
                nsave = nact; % For debugging only
                l = iact(icon); % For debugging only

            end

            %====================%
            % Calculation starts %
            %====================%

            % The following instructions rearrange the active constraints so that the new value of IACT(NACT) is
            % the old value of IACT(ICON). QREXC implements the updates of QFAC and RFAC by a sequence of Givens
            % rotations. Then NACT is reduced by one.

            [qfac, R_slice] = powalg_obj.qrexc_Rfull(qfac, rfac(:, 1:nact), icon); rfac(:, 1:nact) = R_slice; % QREXC does nothing if ICON == NACT.
            % Indeed, it suffices to pass QFAC(:, 1:NACT) and RFAC(1:NACT, 1:NACT) to QREXC as follows. However,
            % compilers may create a temporary copy of RFAC(1:NACT, 1:NACT), which is not contiguous in memory.
            % %call qrexc(qfac(:, 1:nact), rfac(1:nact, 1:nact), icon)

            iact(icon:nact) = [reshape(iact(icon + 1:nact), [], 1); iact(icon)];
            resact(icon:nact) = [reshape(resact(icon + 1:nact), [], 1); resact(icon)];
            resnew(iact(nact)) = max(resact(nact), consts_obj.TINYCV);
            vlam(icon:nact) = [reshape(vlam(icon + 1:nact), [], 1); vlam(icon)];
            nact = nact - 1;

            %====================%
            %  Calculation ends  %
            %====================%

            % Postconditions
            if consts_obj.DEBUGGING
                debug_obj.assert(nact == nsave - 1, "NACT = NSAVE - 1", srname);
                debug_obj.assert(nact >= 0 && nact <= min(m, n) - 1, "1 <= NACT <= MIN(M, N)-1", srname);
                debug_obj.assert(all(iact(1:nact) >= 1 & iact(1:nact) <= m, 'all'), "1 <= IACT <= M", srname);
                debug_obj.assert(~any(iact(1:nact) == l, 'all'), "L is not in IACT(1:NACT)", srname);
                debug_obj.assert(size(qfac, 1) == n && size(qfac, 2) == n, "SIZE(QFAC) == [N, N]", srname);
                debug_obj.assert(linalg_obj.isorth(qfac, 'tol', tol), "QFAC is orthogonal", srname);
                debug_obj.assert(size(rfac, 1) == n && size(rfac, 2) == n, "SIZE(RFAC) == [N, N]", srname);
                debug_obj.assert(linalg_obj.istriu(rfac), "RFAC is upper triangular", srname);
                debug_obj.assert(numel(resact) == m, "SIZE(RESACT) == M", srname);
                debug_obj.assert(numel(resnew) == m, "SIZE(RESNEW) == M", srname);
            end

        end

    end
end