classdef update_uobyqa_mod
    %--------------------------------------------------------------------------------------------------%
    % This module provides subroutines concerning the updates when XPT(:, KNEW) becomes XNEW = XOPT + D.
    %
    % Coded by Zaikun ZHANG (www.zhangzk.net) based on Powell's code and the UOBYQA paper.
    %
    % Dedicated to the late Professor M. J. D. Powell FRS (1936--2015).
    %
    % Started: July 2020
    %
    % Last Modified: Thu 14 Aug 2025 07:37:19 AM CST
    %--------------------------------------------------------------------------------------------------%

    methods
        function [kopt, fval, pl, pq, xpt] = update(~, knew, d, f, moderr, kopt, fval, pl, pq, xpt)
            %--------------------------------------------------------------------------------------------------%
            % This subroutine updates PL, PQ, XPT, KOPT, and FVAL when XPT(:, KNEW) becomes XNEW.
            % See Section 4 of the UOBYQA paper.
            %--------------------------------------------------------------------------------------------------%

            % Common modules
            consts_obj = prima_mat.common.consts_mod();
            debug_obj = prima_mat.common.debug_mod();
            infnan_obj = prima_mat.common.infnan_mod();
            linalg_obj = prima_mat.common.linalg_mod();
            powalg_obj = prima_mat.common.powalg_mod();

            % Inputs

            % D(N)



            % In-outputs

            % FVAL(NPT)
            % PL(NPT-1, NPT)
            % PQ(NPT-1)
            % XPT(N, NPT)

            % Local variables
            srname = "UPDATE";
            n = NaN;
            npt = NaN;
            plnew = NaN(size(pl, 1), 1);
            vlag = NaN(size(xpt, 2), 1);

            % Sizes
            n = size(xpt, 1);
            npt = size(xpt, 2);

            % Preconditions
            if consts_obj.DEBUGGING
                debug_obj.assert(npt == (n + 1) * (n + 2) / 2, "NPT = (N+1)(N+2)/2", srname);
                debug_obj.assert(knew >= 0 && knew <= npt, "0 <= KNEW <= NPT", srname);
                debug_obj.assert(kopt >= 1 && kopt <= npt, "1 <= KOPT <= NPT", srname);
                debug_obj.assert(knew >= 1 || f >= fval(kopt), "KNEW >= 1 unless X is not improved", srname);
                debug_obj.assert(knew ~= kopt || f < fval(kopt), "KNEW /= KOPT unless X is improved", srname);
                debug_obj.assert(numel(d) == n && all(infnan_obj.is_finite(d), 'all'), "SIZE(D) == N, D is finite", srname);
                debug_obj.assert(~(infnan_obj.is_nan_sp(f) || infnan_obj.is_posinf(f)), "F is not NaN or +Inf", srname);
                debug_obj.assert(~any(fval < fval(kopt), 'all'), "FVAL(KOPT) = MINVAL(FVAL)", srname);
                debug_obj.assert(all(infnan_obj.is_finite(xpt), 'all'), "XPT is finite", srname);
                debug_obj.assert(size(pl, 1) == npt - 1 && size(pl, 2) == npt, "SIZE(PL) == [NPT-1, NPT]", srname);
                debug_obj.assert(numel(pq) == npt - 1, "SIZE(PQ) == NPT-1", srname);
            end

            %====================%
            % Calculation starts %
            %====================%

            % Do essentially nothing when KNEW is 0. This can only happen after a trust-region step.
            if knew <= 0
                % KNEW < 0 is impossible if the input is correct.
                return
            end

            % Update the Lagrange functions.
            vlag(:) = powalg_obj.calvlag_qint(pl, d, xpt(:, kopt), kopt);
            pl(:, knew) = pl(:, knew) ./ vlag(knew);
            plnew(:) = pl(:, knew);
            pl(:, :) = pl - linalg_obj.outprod(plnew, vlag);
            % N.B.: The use of OUTPROD is expensive memory-wise, but it is not our concern in this implementation.
            pl(:, knew) = plnew;

            % Update the quadratic model.
            pq(:) = pq + moderr * plnew;

            % Replace the interpolation point that has index KNEW by the point XNEW.
            xpt(:, knew) = xpt(:, kopt) + d;
            fval(knew) = f;

            % KOPT is NOT identical to MINLOC(FVAL). Indeed, if FVAL(KNEW) = FVAL(KOPT) and KNEW < KOPT, then
            % MINLOC(FVAL) = KNEW /= KOPT. Do not change KOPT in this case.
            if f < fval(kopt)
                kopt = knew;
            end

            %====================%
            %  Calculation ends  %
            %====================%

            % Postconditions
            if consts_obj.DEBUGGING
                debug_obj.assert(size(xpt, 1) == n && size(xpt, 2) == npt && all(infnan_obj.is_finite(xpt), 'all'), "SIZE(XPT) == [N, NPT], XPT is finite", srname);
                debug_obj.assert(kopt >= 1 && kopt <= npt, "1 <= KOPT <= NPT", srname);
                debug_obj.assert(~any(fval < fval(kopt), 'all'), "FVAL(KOPT) = MINVAL(FVAL)", srname);
                debug_obj.assert(size(pl, 1) == npt - 1 && size(pl, 2) == npt, "SIZE(PL) == [NPT-1, NPT]", srname);
                debug_obj.assert(numel(pq) == npt - 1, "SIZE(PQ) == NPT-1", srname);
            end

        end

    end
end