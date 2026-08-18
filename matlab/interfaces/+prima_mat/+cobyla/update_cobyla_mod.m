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


            % CONSTR(M)


            % D(N)


            % CONMAT(M, N+1)
            % CVAL(N+1)
            % FVAL(N+1)
            % SIM(N, N+1)
            % SIMI(N, N)


            sim_old = NaN(size(sim, 1), size(sim, 2));
            simi_jdrop = NaN(size(simi, 2), 1);
            simi_old = NaN(size(simi, 1), size(simi, 2));
            simi_test = NaN(size(simi, 1), size(simi, 2));
            simid = NaN(size(simi, 1), 1);
            sum_simi = NaN(size(simi, 2), 1);
            itol = 1.0;

            n = size(sim, 1);

            %====================%
            % Calculation starts %
            %====================%

            % Do nothing when JDROP is 0. This can only happen after a trust-region step.
            if jdrop <= 0
                % JDROP < 0 is impossible if the input is correct.
                info = 0; % INFO must be set, as it is an output!
                return
            end

            sim_old(:, :) = sim;
            simi_old(:, :) = simi;
            % N.B.: The use of OUTPROD is expensive memory-wise, but it is not our concern in this implementation.
            if jdrop <= n
                sim(:, jdrop) = d;
                simi_jdrop(:) = simi(jdrop, :) ./ sum(simi(jdrop, :).' .* d, 'all');
                simi(:, :) = simi - simi * d * simi_jdrop.';
                simi(jdrop, :) = simi_jdrop;
            else                % JDROP = N+1
                sim(:, n + 1) = sim(:, n + 1) + d;
                sim(:, 1:n) = sim(:, 1:n) - d;
                simid(:) = simi * d;
                sum_simi(:) = sum(simi, 1);
                simi(:, :) = simi + simid * (sum_simi ./ (1.0 - sum(simid, 'all'))).';
            end

            % Check whether SIMI is a poor approximation to the inverse of SIM(:, 1:N).
            % Calculate SIMI from scratch if the current one is damaged by rounding errors.
            erri = max(abs(simi * sim(:, 1:n) - eye(n)), [], 'all'); % MAXIMUM(X) returns NaN if X contains NaN
            if erri > 0.1 * itol || isnan(erri)
                simi_test(:, :) = inv(sim(:, 1:n));
                erri_test = max(abs(simi_test * sim(:, 1:n) - eye(n)), [], 'all');
                if erri_test < erri || (isnan(erri) && ~isnan(erri_test))
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
                info = 7;
                sim(:, :) = sim_old;
                simi(:, :) = simi_old;
            end

            %====================%
            %  Calculation ends  %
            %====================%


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


            % CONMAT(M, N+1)
            % CVAL(N+1)
            % FVAL(N+1)
            % SIM(N, N+1)
            % SIMI(N, N)


            sim_jopt = NaN(size(sim, 1), 1);
            sim_old = NaN(size(sim, 1), size(sim, 2));
            simi_old = NaN(size(simi, 1), size(simi, 2));
            simi_test = NaN(size(simi, 1), size(simi, 2));
            itol = 1.0;

            n = size(sim, 1);

            %====================%
            % Calculation starts %
            %====================%

            % INFO must be set, as it is an output.
            info = 0;

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
                sim(:, jopt) = 0.0;
                sim(:, 1:n) = sim(:, 1:n) - sim_jopt;
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
            erri = max(abs(simi * sim(:, 1:n) - eye(n)), [], 'all'); % MAXIMUM(X) returns NaN if X contains NaN
            if erri > 0.1 * itol || isnan(erri)
                simi_test(:, :) = inv(sim(:, 1:n));
                erri_test = max(abs(simi_test * sim(:, 1:n) - eye(n)), [], 'all');
                if erri_test < erri || (isnan(erri) && ~isnan(erri_test))
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
                info = 7;
                sim(:, :) = sim_old;
                simi(:, :) = simi_old;
            end

            %====================%
            %  Calculation ends  %
            %====================%


        end
        function jopt = findpole(~, cpen, cval, fval)
            %--------------------------------------------------------------------------------------------------%
            % This subroutine identifies the best vertex of the current simplex with respect to the merit
            % function PHI = F + CPEN * CSTRV.
            %--------------------------------------------------------------------------------------------------%


            % CVAL(N+1)
            % FVAL(N+1)


            jopt = NaN;

            phi = NaN(numel(cval), 1);

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


        end

    end
end