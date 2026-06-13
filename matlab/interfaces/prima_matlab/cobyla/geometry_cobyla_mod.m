classdef geometry_cobyla_mod
    %--------------------------------------------------------------------------------------------------%
    % This module contains subroutines concerning the geometry-improving of the interpolation set.
    %
    % Coded by Zaikun ZHANG (www.zhangzk.net) based on Powell's code and the COBYLA paper.
    %
    % Dedicated to the late Professor M. J. D. Powell FRS (1936--2015).
    %
    % Started: July 2021
    %
    % Last Modified: Sunday, April 21, 2024 PM03:25:55
    %--------------------------------------------------------------------------------------------------%

    methods
        function jdrop = setdrop_tr(~, ximproved, d, delta, rho, sim, simi)
            %--------------------------------------------------------------------------------------------------%
            % This subroutine finds (the index) of a current interpolation point to be replaced with the
            % trust-region trial point. See (19)--(22) of the COBYLA paper.
            % N.B.:
            % 1. If XIMPROVED == TRUE, then JDROP > 0 so that D is included into XPT. Otherwise, it is a bug.
            % 2. COBYLA never sets JDROP = N+1.
            % TODO: Check whether it improves the performance if JDROP = N+1 is allowed when XIMPROVED is TRUE.
            % Note that UPDATEXFC should be revised accordingly.
            %--------------------------------------------------------------------------------------------------%

            % Common modules
            consts_obj = consts_mod();
            linalg_obj = linalg_mod();
            infnan_obj = infnan_mod();
            debug_obj = debug_mod();


            % Inputs

            % D(N)


            % SIM(N, N+1)
            % SIMI(N, N)

            % Outputs
            jdrop = NaN;

            % Local variables
            srname = "SETDROP_TR";
            n = NaN;
            distsq = NaN(size(sim, 2), 1);
            weight = NaN(size(sim, 2), 1);
            score = NaN(size(sim, 2), 1);
            simid = NaN(size(simi, 1), 1);
            %real(RP) :: sigbar(size(sim, 1))
            %real(RP) :: veta(size(sim, 1))
            %real(RP) :: vsig(size(sim, 1))
            itol = consts_obj.TENTH;

            % Sizes
            n = size(sim, 1);

            % Preconditions
            if consts_obj.DEBUGGING
                debug_obj.assert(n >= 1, "N >= 1", srname);
                debug_obj.assert(numel(d) == n && all(infnan_obj.is_finite(d), 'all'), "SIZE(D) == N, D is finite", srname);
                debug_obj.assert(delta >= rho && rho > 0, "DELTA >= RHO > 0", srname);
                debug_obj.assert(size(sim, 1) == n && size(sim, 2) == n + 1, "SIZE(SIM) == [N, N+1]", srname);
                debug_obj.assert(all(infnan_obj.is_finite(sim), 'all'), "SIM is finite", srname);
                debug_obj.assert(all(sum(abs(sim(:, 1:n)), 1) > 0, 'all'), "SIM(:, 1:N) has no zero column", srname);
                debug_obj.assert(size(simi, 1) == n && size(simi, 2) == n, "SIZE(SIMI) == [N, N]", srname);
                debug_obj.assert(all(infnan_obj.is_finite(simi), 'all'), "SIMI is finite", srname);
                debug_obj.assert(linalg_obj.isinv(sim(:, 1:n), simi, 'tol', itol), "SIMI = SIM(:, 1:N)^{-1}", srname);
            end

            %====================%
            % Calculation starts %
            %====================%

            %--------------------------------------------------------------------------------------------------%
            % The following code is Powell's scheme for defining JDROP.
            %--------------------------------------------------------------------------------------------------%
            %% JDROP = 0 by default. It cannot be removed, as JDROP may not be set below in some cases (e.g.,
            %% when XIMPROVED == FALSE, MAXVAL(ABS(SIMID)) <= 1, and MAXVAL(VETA) <= EDGMAX).
            %jdrop = 0
            %
            %% SIMID(J) is the value of the J-th Lagrange function at D. It is the counterpart of VLAG in UOBYQA
            %% and DEN in NEWUOA/BOBYQA/LINCOA, but it excludes the value of the (N+1)-th Lagrange function.
            %simid = matprod(simi, d)
            %if (any(abs(simid) > 1) .or. (ximproved .and. any(.not. is_nan(simid)))) then
            %    jdrop = int(maxloc(abs(simid), mask=(.not. is_nan(simid)), dim=1), kind(jdrop))
            %    %%MATLAB: [~, jdrop] = max(simid, [], 'omitnan');
            %end if
            %
            %% VETA(J) is the distance from the J-th vertex of the simplex to the best vertex, taking the trial
            %% point SIM(:, N+1) + D into account.
            %if (ximproved) then
            %    veta = sqrt(sum((sim(:, 1:n) - spread(d, dim=2, ncopies=n))**2, dim=1))
            %    %%MATLAB: veta = sqrt(sum((sim(:, 1:n) - d).^2));  % d should be a column! Implicit expansion
            %else
            %    veta = sqrt(sum(sim(:, 1:n)**2, dim=1))
            %end if
            %
            %% VSIG(J) (J=1, .., N) is the Euclidean distance from vertex J to the opposite face of the simplex.
            %vsig = ONE / sqrt(sum(simi**2, dim=2))
            %sigbar = abs(simid) * vsig
            %
            %% The following JDROP will overwrite the previous one if its premise holds. FACTOR_DELTA = 1.1
            %% and FACTOR_ALPHA = 0.25.
            %mask = (veta > factor_delta * delta .and. (sigbar >= factor_alpha * delta .or. sigbar >= vsig))
            %if (any(mask)) then
            %    jdrop = int(maxloc(veta, mask=mask, dim=1), kind(jdrop))
            %    %%MATLAB: etamax = max(veta(mask)); jdrop = find(mask & ~(veta < etamax), 1, 'first');
            %end if
            %
            %% Powell's code does not include the following instructions. With Powell's code, if SIMID consists
            %% of only NaN, then JDROP can be 0 even when XIMPROVED == TRUE (i.e., D reduces the merit function).
            %% With the following code, JDROP cannot be 0 when XIMPROVED == TRUE, unless VETA is all NaN, which
            %% should not happen if X0 does not contain NaN, the trust-region/geometry steps never contain NaN,
            %% and we exit once encountering an iterate containing Inf (due to overflow).
            %if (ximproved .and. jdrop <= 0) then  ! Write JDROP <= 0 instead of JDROP == 0 for robustness.
            %    jdrop = int(maxloc(veta, mask=(.not. is_nan(veta)), dim=1), kind(jdrop))
            %    %%MATLAB: [~, jdrop] = max(veta, [], 'omitnan');
            %end if
            %--------------------------------------------------------------------------------------------------%
            % Powell's scheme ends here.
            %--------------------------------------------------------------------------------------------------%


            % The following definition of JDROP is inspired by SETDROP_TR in UOBYQA/NEWUOA/BOBYQA/LINCOA.
            % It is simpler and works better than Powell's scheme. Note that we allow JDROP to be N+1 if
            % IMPROVEX is TRUE, whereas Powell's code does not.
            % See also (4.1) of Scheinberg-Toint-2010: Self-Correcting Geometry in Model-Based Algorithms for
            % Derivative-Free Unconstrained Optimization, which refers to the strategy here as the "combined
            % distance/poisedness criteria".

            % DISTQ(J) is the square of the distance from the J-th vertex of the simplex to the "best" point so
            % far, taking the trial point SIM(:, N+1) + D into account.
            if ximproved
                distsq(1:n) = sum(fortran.power((sim(:, 1:n) - fortran.spread(d, 'dim', 2, 'ncopies', n)), 2), 1);
                %%MATLAB: distsq = sum((sim(:, 1:n) - d).^2);  % d should be a column! Implicit expansion
                distsq(n + 1) = sum(fortran.power(d, 2), 'all');
            else
                distsq(1:n) = sum(fortran.power(sim(:, 1:n), 2), 1);
                distsq(n + 1) = consts_obj.ZERO;
            end

            weight(:) = max(consts_obj.ONE, distsq ./ max(rho, consts_obj.TENTH * delta) ^ 2); % Similar to Powell's NEWUOA code
            % Other possible definitions of WEIGHT.
            % %weight = distsq  ! Similar to Powell's LINCOA code, but WRONG. See comments in LINCOA/geometry.f90.
            % %weight = max(ONE, 25.0_RP * distsq / delta**2)  ! Similar to Powell's BOBYQA code, works well
            % %weight = max(ONE, TEN * distsq / delta**2)
            % %weight = max(ONE, 1.0E2_RP * distsq / delta**2)
            % %weight = max(ONE, distsq / rho**2)  ! Similar to Powell's UOBYQA

            % If 1 <= J <= N, SIMID(J) is the value of the J-th Lagrange function at D; the value of the
            % (N+1)-th Lagrange function is 1 - SUM(SIMID). [SIMID, 1 - SUM(SIMID)] is the counterpart of
            % VLAG in UOBYQA and DEN in NEWUOA/BOBYQA/LINCOA.
            simid(:) = linalg_obj.matprod21(simi, d);
            score(:) = weight .* abs([reshape(simid, [], 1); consts_obj.ONE - sum(simid, 'all')]);

            % If XIMPROVED = FALSE (D does not render a better X), set SCORE(N+1) = -1 to avoid JDROP = N+1.
            if ~ximproved
                score(n + 1) = -consts_obj.ONE;
            end

            % SCORE(J) is NaN implies SIMID(J) is NaN, but we want ABS(SIMID) to be big. So we exclude such J.
            score(linalg_obj.trueloc(infnan_obj.is_nan(score))) = -consts_obj.ONE;

            jdrop = 0;
            % The following IF works a bit better than `IF (ANY(SCORE > 1) .OR. ANY(SCORE > 0) .AND. XIMPROVED)`
            % from Powell's UOBYQA and NEWUOA code.
            if any(score > 0, 'all')
                % Powell's BOBYQA and LINCOA code
                jdrop = fix(fortran.maxloc(score, 'dim', 1));
                %%MATLAB: [~, jdrop] = max(score);

            end

            if (ximproved && jdrop == 0) || jdrop < 0
                % JDROP < 0 is impossible in theory.
                jdrop = fix(fortran.maxloc(distsq, 'dim', 1));
            end

            %====================%
            %  Calculation ends  %
            %====================%

            % Postconditions
            if consts_obj.DEBUGGING
                debug_obj.assert(jdrop >= 0 && jdrop <= n + 1, "0 <= JDROP <= N+1", srname);
                debug_obj.assert(jdrop <= n || ximproved, "JDROP <= n unless IMPROVEX = TRUE", srname);
                debug_obj.assert(jdrop >= 1 || ~ximproved, "JDROP >= 1 unless IMPROVEX = FALSE", srname);
                % JDROP >= 1 when XIMPROVED = TRUE unless NaN occurs in DISTSQ, which should not happen if the
                % starting point does not contain NaN and the trust-region/geometry steps never contain NaN.

            end

        end
        function d = geostep(~, jdrop, amat, bvec, conmat, cpen, cval, delbar, fval, simi)
            %--------------------------------------------------------------------------------------------------%
            % This function calculates a geometry step so that the geometry of the interpolation set is improved
            % when SIM(:, JDRO_GEO) is replaced with SIM(:, N+1) + D. See (15)--(17) of the COBYLA paper.
            %--------------------------------------------------------------------------------------------------%

            % Common modules
            consts_obj = consts_mod();
            debug_obj = debug_mod();
            infnan_obj = infnan_mod();
            linalg_obj = linalg_mod();


            % Inputs



            % CONMAT(M, N+1)

            % CVAL(N+1)

            % FVAL(N+1)
            % SIMI(N, N)

            % Outputs
            d = NaN(size(simi, 1), 1); % D(N)

            % Local variables
            srname = "GEOSTEP";
            m = NaN;
            m_lcon = NaN;
            n = NaN;
            A = NaN(size(simi, 1), size(conmat, 1));
            cvnd = NaN;
            cvpd = NaN;
            g = NaN(size(simi, 1), 1);

            % Sizes
            m_lcon = fix(numel(bvec));
            m = size(conmat, 1);
            n = size(simi, 1);

            % Preconditions
            if consts_obj.DEBUGGING
                debug_obj.assert(m >= m_lcon && m >= 0, "M >= 0", srname);
                debug_obj.assert(n >= 1, "N >= 1", srname);
                debug_obj.assert(delbar > 0, "DELBAR > 0", srname);
                debug_obj.assert(cpen > 0, "CPEN > 0", srname);
                debug_obj.assert(size(simi, 1) == n && size(simi, 2) == n, "SIZE(SIMI) == [N, N]", srname);
                debug_obj.assert(all(infnan_obj.is_finite(simi), 'all'), "SIMI is finite", srname);
                debug_obj.assert(numel(fval) == n + 1 && ~any(infnan_obj.is_nan(fval) | infnan_obj.is_posinf(fval), 'all'), "SIZE(FVAL) == NPT and FVAL is not NaN/+Inf", srname);
                debug_obj.assert(size(conmat, 1) == m && size(conmat, 2) == n + 1, "SIZE(CONMAT) == [M, N+1]", srname);
                debug_obj.assert(~any(infnan_obj.is_nan(conmat) | infnan_obj.is_posinf(conmat), 'all'), "CONMAT does not contain NaN/+Inf", srname);
                debug_obj.assert(numel(cval) == n + 1 && ~any(cval < 0 | infnan_obj.is_nan(cval) | infnan_obj.is_posinf(cval), 'all'), "SIZE(CVAL) == NPT and CVAL does not contain negative NaN/+Inf", srname);
                debug_obj.assert(jdrop >= 1 && jdrop <= n, "1 <= JDROP <= N", srname);
            end

            %====================%
            % Calculation starts %
            %====================%

            % SIMI(JDROP, :) is a vector perpendicular to the face of the simplex to the opposite of vertex
            % JDROP. Set D to the vector in this direction and with length DELBAR.
            d(:) = simi(jdrop, :);
            d(:) = delbar * (d ./ linalg_obj.p_norm(d));

            % The code below chooses the direction of D according to an approximation of the merit function.
            % See (17) of the COBYLA paper and  line 225 of Powell's cobylb.f.

            % Calculate the coefficients of the linear approximations to the objective and constraint functions.
            % N.B.: CONMAT and SIMI have been updated after the last trust-region step, but G and A have not.
            % So we cannot pass G and A from outside.
            g(:) = linalg_obj.matprod12(fval(1:n) - fval(n + 1), simi);
            A(:, 1:m_lcon) = amat;
            A(:, m_lcon + 1:m) = linalg_obj.matprod22(conmat(m_lcon + 1:m, 1:n) - fortran.spread(conmat(m_lcon + 1:m, n + 1), 'dim', 2, 'ncopies', n), simi)';
            %%MATLAB: A(:, m_lcon+1:m) = simi'*(conmat(m_lcon+1:m, 1:n) - conmat(m_lcon+1:m, n+1))' % Implicit expansion for subtraction
            % CVPD and CVND are the predicted constraint violation of D and -D by the linear models.
            cvpd = linalg_obj.maximum1([consts_obj.ZERO; reshape(conmat(:, n + 1) + linalg_obj.matprod12(d, A), [], 1)]);
            cvnd = linalg_obj.maximum1([consts_obj.ZERO; reshape(conmat(:, n + 1) - linalg_obj.matprod12(d, A), [], 1)]);
            % Take -D if the linear models predict that its merit function value is lower.
            if -linalg_obj.inprod(d, g) + cpen * cvnd < linalg_obj.inprod(d, g) + cpen * cvpd
                d(:) = -d;
            end

            %====================%
            %  Calculation ends  %
            %====================%

            % Postconditions
            if consts_obj.DEBUGGING
                debug_obj.assert(numel(d) == n && all(infnan_obj.is_finite(d), 'all'), "SIZE(D) == N, D is finite", srname);
                % In theory, ||S|| == DELBAR, which may be false due to rounding, but not too far.
                % It is crucial to ensure that the geometry step is nonzero, which holds in theory.
                debug_obj.assert(linalg_obj.p_norm(d) > 0.9 * delbar && linalg_obj.p_norm(d) <= 1.1 * delbar, "||D|| == DELBAR", srname);
            end
        end

    end
end