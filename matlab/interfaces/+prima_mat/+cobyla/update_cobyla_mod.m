classdef update_cobyla_mod
    %--------------------------------------------------------------------------------------------------%
    % This module contains subroutines concerning the update of the interpolation set.
    %
    % Coded by Zaikun ZHANG (www.zhangzk.net) based on Powell's code and the COBYLA paper.
    %
    % Dedicated to the late Professor M. J. D. Powell FRS (1936--2015).
    %
    % Started: July 2021
    %
    % Last Modified: Thu 14 Aug 2025 07:34:04 AM CST
    %--------------------------------------------------------------------------------------------------%

    methods
        function [conmat, cval, fval, sim, simi, info] = updatexfc(obj, jdrop, constr, cpen, cstrv, d, f, conmat, cval, fval, sim, simi)
            %--------------------------------------------------------------------------------------------------%
            % This subroutine revises the simplex by updating the elements of SIM, SIMI, FVAL, CONMAT, and CVAL.
            %--------------------------------------------------------------------------------------------------%

            % Common modules
            consts_obj = prima_mat.common.consts_mod();
            infnan_obj = prima_mat.common.infnan_mod();
            infos_obj = prima_mat.common.infos_mod();
            linalg_obj = prima_mat.common.linalg_mod();
            debug_obj = prima_mat.common.debug_mod();


            % Inputs

            % CONSTR(M)


            % D(N)


            % In-outputs
            % CONMAT(M, N+1)
            % CVAL(N+1)
            % FVAL(N+1)
            % SIM(N, N+1)
            % SIMI(N, N)

            % Outputs


            % Local variables
            srname = "UPDATEXFC";


            sim_old = NaN(size(sim, 1), size(sim, 2));
            simi_jdrop = NaN(size(simi, 2), 1);
            simi_old = NaN(size(simi, 1), size(simi, 2));
            simi_test = NaN(size(simi, 1), size(simi, 2));
            simid = NaN(size(simi, 1), 1);
            sum_simi = NaN(size(simi, 2), 1);
            itol = consts_obj.ONE;

            % Sizes
            m = numel(constr);
            n = size(sim, 1);

            % Preconditions
            if consts_obj.DEBUGGING
                debug_obj.assert(m >= 0, "M >= 0", srname);
                debug_obj.assert(n >= 1, "N >= 1", srname);
                debug_obj.assert(jdrop >= 0 && jdrop <= n + 1, "1 <= JDROP <= N+1", srname);
                debug_obj.assert(~any(infnan_obj.is_nan_sp(constr) | infnan_obj.is_posinf(constr), 'all'), "CONSTR does not contain NaN/+Inf", srname);
                debug_obj.assert(~(infnan_obj.is_nan_sp(cstrv) || infnan_obj.is_posinf(cstrv)), "CSTRV is not NaN/+Inf", srname);
                debug_obj.assert(numel(d) == n && all(infnan_obj.is_finite(d), 'all'), "SIZE(D) == N, D is finite", srname);
                debug_obj.assert(~(infnan_obj.is_nan_sp(f) || infnan_obj.is_posinf(f)), "F is not NaN/+Inf", srname);
                debug_obj.assert(size(conmat, 1) == m && size(conmat, 2) == n + 1, "SIZE(CONMAT) = [M, N+1]", srname);
                debug_obj.assert(~any(infnan_obj.is_nan_sp(conmat) | infnan_obj.is_posinf(conmat), 'all'), "CONMAT does not contain NaN/+Inf", srname);
                debug_obj.assert(numel(cval) == n + 1 && ~any(cval < 0 | infnan_obj.is_nan_sp(cval) | infnan_obj.is_posinf(cval), 'all'), "SIZE(CVAL) == N+1 and CVAL does not contain negative values or NaN/+Inf", srname);
                debug_obj.assert(numel(fval) == n + 1 && ~any(infnan_obj.is_nan_sp(fval) | infnan_obj.is_posinf(fval), 'all'), "SIZE(FVAL) == N+1 and FVAL is not NaN/+Inf", srname);
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

            % Do nothing when JDROP is 0. This can only happen after a trust-region step.
            if jdrop <= 0
                % JDROP < 0 is impossible if the input is correct.
                info = infos_obj.INFO_DFT; % INFO must be set, as it is an output!
                return
            end

            sim_old(:, :) = sim;
            simi_old(:, :) = simi;
            % N.B.: The use of OUTPROD is expensive memory-wise, but it is not our concern in this implementation.
            if jdrop <= n
                sim(:, jdrop) = d;
                simi_jdrop(:) = simi(jdrop, :) ./ linalg_obj.inprod(simi(jdrop, :), d);
                simi(:, :) = simi - linalg_obj.outprod(linalg_obj.matprod21(simi, d), simi_jdrop);
                simi(jdrop, :) = simi_jdrop;
            else                % JDROP = N+1
                sim(:, n + 1) = sim(:, n + 1) + d;
                sim(:, 1:n) = sim(:, 1:n) - reshape(d, [], 1);
                simid(:) = linalg_obj.matprod21(simi, d);
                sum_simi(:) = sum(simi, 1);
                simi(:, :) = simi + linalg_obj.outprod(simid, sum_simi ./ (consts_obj.ONE - sum(simid, 'all')));
            end

            % Check whether SIMI is a poor approximation to the inverse of SIM(:, 1:N).
            % Calculate SIMI from scratch if the current one is damaged by rounding errors.
            erri = linalg_obj.maximum2(abs(linalg_obj.matprod22(simi, sim(:, 1:n)) - linalg_obj.eye1(n))); % MAXIMUM(X) returns NaN if X contains NaN
            if erri > consts_obj.TENTH * itol || infnan_obj.is_nan_sp(erri)
                simi_test(:, :) = linalg_obj.inv(sim(:, 1:n));
                erri_test = linalg_obj.maximum2(abs(linalg_obj.matprod22(simi_test, sim(:, 1:n)) - linalg_obj.eye1(n)));
                if erri_test < erri || (infnan_obj.is_nan_sp(erri) && ~infnan_obj.is_nan_sp(erri_test))
                    simi(:, :) = simi_test;
                    erri = erri_test;
                end
            end

            % If SIMI is satisfactory, then update FVAL, CONMAT, CVAL, and the pole position. Otherwise, restore
            % SIM and SIMI, and return with INFO = DAMAGING_ROUNDING.
            if erri <= itol
                fval(jdrop) = f;
                conmat(:, jdrop) = constr;
                cval(jdrop) = cstrv;
                % Switch the best vertex to the pole position SIM(:, N+1) if it is not there already.
                [conmat, cval, fval, sim, simi, info] = obj.updatepole(cpen, conmat, cval, fval, sim, simi);
            else                % ERRI > ITOL or ERRI is NaN
                info = infos_obj.DAMAGING_ROUNDING;
                sim(:, :) = sim_old;
                simi(:, :) = simi_old;
            end

            %====================%
            %  Calculation ends  %
            %====================%

            % Postconditions
            if consts_obj.DEBUGGING
                debug_obj.assert(size(conmat, 1) == m && size(conmat, 2) == n + 1, "SIZE(CONMAT) = [M, N+1]", srname);
                debug_obj.assert(~any(infnan_obj.is_nan_sp(conmat) | infnan_obj.is_posinf(conmat), 'all'), "CONMAT does not contain NaN/+Inf", srname);
                debug_obj.assert(numel(cval) == n + 1 && ~any(cval < 0 | infnan_obj.is_nan_sp(cval) | infnan_obj.is_posinf(cval), 'all'), "SIZE(CVAL) == N+1 and CVAL does not contain negative values or NaN/+Inf", srname);
                debug_obj.assert(numel(fval) == n + 1 && ~any(infnan_obj.is_nan_sp(fval) | infnan_obj.is_posinf(fval), 'all'), "SIZE(FVAL) == N+1 and FVAL is not NaN/+Inf", srname);
                debug_obj.assert(size(sim, 1) == n && size(sim, 2) == n + 1, "SIZE(SIM) == [N, N+1]", srname);
                debug_obj.assert(all(infnan_obj.is_finite(sim), 'all'), "SIM is finite", srname);
                debug_obj.assert(all(sum(abs(sim(:, 1:n)), 1) > 0, 'all'), "SIM(:, 1:N) has no zero column", srname);
                debug_obj.assert(size(simi, 1) == n && size(simi, 2) == n, "SIZE(SIMI) == [N, N]", srname);
                debug_obj.assert(all(infnan_obj.is_finite(simi), 'all'), "SIMI is finite", srname);
                debug_obj.assert(linalg_obj.isinv(sim(:, 1:n), simi, 'tol', itol) || info == infos_obj.DAMAGING_ROUNDING, "SIMI = SIM(:, 1:N)^{-1} unless the rounding is damaging", srname);
            end
        end
        function [conmat, cval, fval, sim, simi, info] = updatepole(obj, cpen, conmat, cval, fval, sim, simi)
            %--------------------------------------------------------------------------------------------------%
            % This subroutine identifies the best vertex of the current simplex with respect to the merit
            % function PHI = F + CPEN * CSTRV, and then switch this vertex to SIM(:, N + 1), which Powell called
            % the "pole position" in his comments. CONMAT, CVAL, FVAL, and SIMI are updated accordingly.
            %
            % N.B. 1: In precise arithmetic, the following two procedures produce the same results:
            % 1) apply UPDATEPOLE to SIM twice, first with CPEN = CPEN1 and then with CPEN = CPEN2;
            % 2) apply UPDATEPOLE to SIM with CPEN = CPEN2.
            % In finite-precision arithmetic, however, they may produce different results unless CPEN1 = CPEN2.
            %
            % N.B. 2: When JOPT == N+1, the best vertex is already at the pole position, so there is nothing to
            % switch. However, as in Powell's code, the code below will check whether SIMI is good enough to
            % work as the inverse of SIM(:, 1:N) or not. If not, Powell's code would invoke an error return of
            % COBYLB; our implementation, however, will try calculating SIMI from scratch; if the recalculated
            % SIMI is still of poor quality, then UPDATEPOLE will return with INFO = DAMAGING_ROUNDING,
            % informing COBYLB that SIMI is poor due to damaging rounding errors.
            %
            % N.B. 3: UPDATEPOLE should be called when and only when FINDPOLE can potentially returns a value
            % other than N+1. The value of FINDPOLE is determined by CPEN, CVAL, and FVAL, the latter two being
            % decided by SIM. Thus UPDATEPOLE should be called after CPEN or SIM changes. COBYLA updates CPEN at
            % only two places: the beginning of each trust-region iteration, and when REDRHO is called;
            % SIM is updated only by UPDATEXFC, which itself calls UPDATEPOLE internally. Therefore, we only
            % need to call UPDATEPOLE after updating CPEN at the beginning of each trust-region iteration and
            % after each invocation of REDRHO.
            %--------------------------------------------------------------------------------------------------%

            % Common modules
            consts_obj = prima_mat.common.consts_mod();
            infos_obj = prima_mat.common.infos_mod();
            debug_obj = prima_mat.common.debug_mod();
            infnan_obj = prima_mat.common.infnan_mod();
            linalg_obj = prima_mat.common.linalg_mod();


            % Inputs


            % In-outputs
            % CONMAT(M, N+1)
            % CVAL(N+1)
            % FVAL(N+1)
            % SIM(N, N+1)
            % SIMI(N, N)

            % Outputs


            % Local variables
            srname = "UPDATEPOLE";


            sim_jopt = NaN(size(sim, 1), 1);
            sim_old = NaN(size(sim, 1), size(sim, 2));
            simi_old = NaN(size(simi, 1), size(simi, 2));
            simi_test = NaN(size(simi, 1), size(simi, 2));
            itol = consts_obj.ONE;

            % Sizes
            m = size(conmat, 1);
            n = size(sim, 1);

            % Preconditions
            if consts_obj.DEBUGGING
                debug_obj.assert(m >= 0, "M >= 0", srname);
                debug_obj.assert(n >= 1, "N >= 1", srname);
                debug_obj.assert(cpen > 0, "CPEN > 0", srname);
                debug_obj.assert(size(conmat, 1) == m && size(conmat, 2) == n + 1, "SIZE(CONMAT) = [M, N+1]", srname);
                debug_obj.assert(~any(infnan_obj.is_nan_sp(conmat) | infnan_obj.is_posinf(conmat), 'all'), "CONMAT does not contain NaN/+Inf", srname);
                debug_obj.assert(numel(cval) == n + 1 && ~any(cval < 0 | infnan_obj.is_nan_sp(cval) | infnan_obj.is_posinf(cval), 'all'), "SIZE(CVAL) == N+1 and CVAL does not contain negative values or NaN/+Inf", srname);
                debug_obj.assert(numel(fval) == n + 1 && ~any(infnan_obj.is_nan_sp(fval) | infnan_obj.is_posinf(fval), 'all'), "SIZE(FVAL) == N+1 and FVAL is not NaN/+Inf", srname);
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

            % INFO must be set, as it is an output.
            info = infos_obj.INFO_DFT;

            % Identify the optimal vertex of the current simplex.
            jopt = obj.findpole(cpen, cval, fval);

            % Switch the best vertex to the pole position SIM(:, N+1) if it is not there already, and update
            % SIMI. Before the update, save a copy of SIM and SIMI. If the update is unsuccessful due to
            % damaging rounding errors, we restore them and return with INFO = DAMAGING_ROUNDING.
            sim_old(:, :) = sim;
            simi_old(:, :) = simi;
            if jopt >= 1 && jopt <= n
                % Unless there is a bug in FINDPOLE, it is guaranteed that JOPT >= 1.
                % When JOPT == N + 1, there is nothing to switch; in addition, SIMI(JOPT, :) will be illegal.
                sim(:, n + 1) = sim(:, n + 1) + sim(:, jopt);
                sim_jopt(:) = sim(:, jopt);
                sim(:, jopt) = consts_obj.ZERO;
                sim(:, 1:n) = sim(:, 1:n) - reshape(sim_jopt, [], 1);
                %%MATLAB: sim(:, 1:n) = sim(:, 1:n) - sim_jopt; % sim_jopt should be a column! Implicit expansion
                % The above update is equivalent to multiply SIM(:, 1:N) from the right side by a matrix whose
                % JOPT-th row is [-1, -1, ..., -1], while all the other rows are the same as those of the
                % identity matrix. It is easy to check that the inverse of this matrix is itself. Therefore,
                % SIMI should be updated by a multiplication with this matrix (i.e., its inverse) from the left
                % side, as is done in the following line. The JOPT-th row of the updated SIMI is minus the sum
                % of all rows of the original SIMI, whereas all the other rows remain unchanged.
                simi(jopt, :) = -sum(simi, 1); % Must ensure that 1 <= JOPT <= N!

            end

            % Check whether SIMI is a poor approximation to the inverse of SIM(:, 1:N).
            % Calculate SIMI from scratch if the current one is damaged by rounding errors.
            erri = linalg_obj.maximum2(abs(linalg_obj.matprod22(simi, sim(:, 1:n)) - linalg_obj.eye1(n))); % MAXIMUM(X) returns NaN if X contains NaN
            if erri > consts_obj.TENTH * itol || infnan_obj.is_nan_sp(erri)
                simi_test(:, :) = linalg_obj.inv(sim(:, 1:n));
                erri_test = linalg_obj.maximum2(abs(linalg_obj.matprod22(simi_test, sim(:, 1:n)) - linalg_obj.eye1(n)));
                if erri_test < erri || (infnan_obj.is_nan_sp(erri) && ~infnan_obj.is_nan_sp(erri_test))
                    simi(:, :) = simi_test;
                    erri = erri_test;
                end
            end

            % If SIMI is satisfactory, then update FVAL, CONMAT, and CVAL. Otherwise, restore SIM and SIMI, and
            % return with INFO = DAMAGING_ROUNDING.
            if erri <= itol
                if jopt >= 1 && jopt <= n
                    fval([jopt, n + 1]) = fval([n + 1, jopt]);
                    conmat(:, [jopt, n + 1]) = conmat(:, [n + 1, jopt]);
                    cval([jopt, n + 1]) = cval([n + 1, jopt]);
                end
            else                % ERRI > ITOL or ERRI is NaN
                info = infos_obj.DAMAGING_ROUNDING;
                sim(:, :) = sim_old;
                simi(:, :) = simi_old;
            end

            %====================%
            %  Calculation ends  %
            %====================%

            % Postconditions
            if consts_obj.DEBUGGING
                debug_obj.assert(obj.findpole(cpen, cval, fval) == n + 1 || info == infos_obj.DAMAGING_ROUNDING, "The best point is SIM(:, N+1) unless the rounding is damaging", srname);
                debug_obj.assert(size(conmat, 1) == m && size(conmat, 2) == n + 1, "SIZE(CONMAT) = [M, N+1]", srname);
                debug_obj.assert(~any(infnan_obj.is_nan_sp(conmat) | infnan_obj.is_posinf(conmat), 'all'), "CONMAT does not contain NaN/+Inf", srname);
                debug_obj.assert(numel(cval) == n + 1 && ~any(cval < 0 | infnan_obj.is_nan_sp(cval) | infnan_obj.is_posinf(cval), 'all'), "SIZE(CVAL) == N+1 and CVAL does not contain negative values or NaN/+Inf", srname);
                debug_obj.assert(numel(fval) == n + 1 && ~any(infnan_obj.is_nan_sp(fval) | infnan_obj.is_posinf(fval), 'all'), "SIZE(FVAL) == N+1 and FVAL is not NaN/+Inf", srname);
                debug_obj.assert(size(sim, 1) == n && size(sim, 2) == n + 1, "SIZE(SIM) == [N, N+1]", srname);
                debug_obj.assert(all(infnan_obj.is_finite(sim), 'all'), "SIM is finite", srname);
                debug_obj.assert(all(sum(abs(sim(:, 1:n)), 1) > 0, 'all'), "SIM(:, 1:N) has no zero column", srname);
                debug_obj.assert(size(simi, 1) == n && size(simi, 2) == n, "SIZE(SIMI) == [N, N]", srname);
                debug_obj.assert(all(infnan_obj.is_finite(simi), 'all'), "SIMI is finite", srname);
                % Do not check SIMI = SIM(:, 1:N)^{-1}, as it may not be true due to damaging rounding.
                debug_obj.assert(linalg_obj.isinv(sim(:, 1:n), simi, 'tol', itol) || info == infos_obj.DAMAGING_ROUNDING, "SIMI = SIM(:, 1:N)^{-1} unless the rounding is damaging", srname);
            end

        end
        function jopt = findpole(~, cpen, cval, fval)
            %--------------------------------------------------------------------------------------------------%
            % This subroutine identifies the best vertex of the current simplex with respect to the merit
            % function PHI = F + CPEN * CSTRV.
            %--------------------------------------------------------------------------------------------------%

            % Common modules
            consts_obj = prima_mat.common.consts_mod();
            debug_obj = prima_mat.common.debug_mod();
            infnan_obj = prima_mat.common.infnan_mod();


            % Inputs

            % CVAL(N+1)
            % FVAL(N+1)

            % Outputs
            jopt = NaN;

            % Local variables
            srname = "FINDPOLE";

            phi = NaN(numel(cval), 1);


            % Size
            n = numel(fval) - 1;

            % Preconditions
            if consts_obj.DEBUGGING
                debug_obj.assert(cpen > 0, "CPEN > 0", srname);
                debug_obj.assert(numel(cval) == n + 1 && ~any(cval < 0 | infnan_obj.is_nan_sp(cval) | infnan_obj.is_posinf(cval), 'all'), "SIZE(CVAL) == N+1 and CVAL does not contain negative values or NaN/+Inf", srname);
                debug_obj.assert(numel(fval) == n + 1 && ~any(infnan_obj.is_nan_sp(fval) | infnan_obj.is_posinf(fval), 'all'), "SIZE(FVAL) == N+1 and FVAL is not NaN/+Inf", srname);
            end

            %====================%
            % Calculation starts %
            %====================%

            % Identify the optimal vertex of the current simplex.
            jopt = numel(fval); % We use N + 1 as the default value of JOPT.
            phi(:) = fval + cpen * cval;
            phimin = min(phi, [], 'all');
            % Essentially, JOPT = MINLOC(PHI). However, we keep JOPT = N + 1 unless there is a strictly better
            % choice. When there are multiple choices, we choose the JOPT with the smallest value of CVAL.
            if phimin < phi(jopt) || any(cval < cval(jopt) & phi <= phi(jopt), 'all')
                jopt = fortran.minloc(cval, 'mask', (phi <= phimin), 'dim', 1);
                %%MATLAB: cmin = min(cval(phi <= phimin)); jopt = find(phi <= phimin & cval <= cmin, 1, 'first');

            end

            %====================%
            %  Calculation ends  %
            %====================%

            % Postconditions
            if consts_obj.DEBUGGING
                debug_obj.assert(jopt >= 1 && jopt <= n + 1, "1 <= JOPT <= N+1", srname);
                debug_obj.assert(jopt == n + 1 || phi(jopt) < phi(n + 1) || (phi(jopt) <= phi(n + 1) && cval(jopt) < cval(n + 1)), "JOPT = N+1 unless PHI(JOPT) < PHI(N+1) or PHI(JOPT) <= PHI(N+1) and CVAL(JOPT) < CVAL(N+1)", srname);
            end
        end

    end
end