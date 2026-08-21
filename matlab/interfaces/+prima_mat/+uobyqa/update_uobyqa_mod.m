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


            powalg_obj = prima_mat.common.powalg_mod();

            vlag = NaN(size(xpt, 2), 1);

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
            plnew = pl(:, knew);
            pl = pl - plnew * vlag.';
            % N.B.: The use of OUTPROD is expensive memory-wise, but it is not our concern in this implementation.
            pl(:, knew) = plnew;

            % Update the quadratic model.
            pq = pq + moderr * plnew;

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


        end

    end
end