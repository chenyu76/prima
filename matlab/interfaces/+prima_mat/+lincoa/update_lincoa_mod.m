classdef update_lincoa_mod
    %--------------------------------------------------------------------------------------------------%
    % This module provides subroutines concerning the updates when XPT(:, KNEW) becomes XNEW = XOPT + D.
    %
    % Coded by Zaikun ZHANG (www.zhangzk.net) based on Powell's LINCOA code.
    %
    % Dedicated to the late Professor M. J. D. Powell FRS (1936--2015).
    %
    % Started: February 2022
    %
    % Last Modified: Friday, March 15, 2024 PM03:35:14
    %--------------------------------------------------------------------------------------------------%

    methods
        function [kopt, fval, xpt] = updatexf(~, knew, ximproved, f, xnew, kopt, fval, xpt)
            %--------------------------------------------------------------------------------------------------%
            % This subroutine updates [XPT, FVAL, KOPT] so that XPT(:, KNEW) is updated to XNEW.
            %--------------------------------------------------------------------------------------------------%


            %====================%
            % Calculation starts %
            %====================%

            % Do essentially nothing when KNEW is 0. This can only happen after a trust-region step.
            if knew <= 0
                % KNEW < 0 is impossible if the input is correct.
                return
            end

            xpt(:, knew) = xnew;
            fval(knew) = f;

            if ximproved
                kopt = knew;
            end

            %====================%
            %  Calculation ends  %
            %====================%


        end
        function [gopt, hq, pq] = updateq(~, idz, knew, ximproved, bmat, d, moderr, xdrop, xosav, xpt, zmat, gopt, hq, pq)
            %--------------------------------------------------------------------------------------------------%
            % This subroutine updates GOPT, HQ, and PQ when XPT(:, KNEW) changes from XDROP to XNEW = XOSAV + D,
            % where XOSAV is the unupdated XOPT, namely the XOPT before UPDATEXF is called.
            % See Section 4 of the NEWUOA paper and that of the BOBYQA paper (there is no LINCOA paper).
            % N.B.:
            % XNEW is encoded in [BMAT, ZMAT, IDZ] after UPDATEH being called, and it also equals XPT(:, KNEW)
            % after UPDATEXF being called. Indeed, we only need BMAT(:, KNEW) instead of the entire matrix.
            %--------------------------------------------------------------------------------------------------%


            linalg_obj = prima_mat.common.linalg_mod();
            powalg_obj = prima_mat.common.powalg_mod();

            pqinc = NaN(size(pq));

            %====================%
            % Calculation starts %
            %====================%

            % Do nothing when KNEW is 0. This can only happen after a trust-region step.
            if knew <= 0
                % KNEW < 0 is impossible if the input is correct.
                return
            end

            % The unupdated model corresponding to [GOPT, HQ, PQ] interpolates F at all points in XPT except for
            % XNEW. The error is MODERR = [F(XNEW)-F(XOPT)] - [Q(XNEW)-Q(XOPT)].

            % Absorb PQ(KNEW)*XDROP*XDROP^T into the explicit part of the Hessian.
            % Implement R1UPDATE properly so that it ensures that HQ is symmetric.
            hq = linalg_obj.r1_sym(hq, pq(knew), xdrop);
            pq(knew) = 0.0;

            % Update the implicit part of the Hessian.
            pqinc(:) = moderr * powalg_obj.omega_col(idz, zmat, knew);
            pq = pq + pqinc;

            % Update the gradient, which needs the updated XPT.
            gopt = gopt + moderr * bmat(:, knew) + powalg_obj.hess_mul(xosav, xpt, pqinc);

            % Further update GOPT if XIMPROVED is TRUE, as XOPT changes from XOSAV to XNEW = XOSAV + D.
            if ximproved
                gopt = gopt + powalg_obj.hess_mul(d, xpt, pq, 'hq', hq);
            end

            %====================%
            %  Calculation ends  %
            %====================%


        end
        function [qalt_better, gopt, pq, hq, galt, pqalt] = tryqalt(~, idz, bmat, fval, xopt, xpt, zmat, qalt_better, gopt, pq, hq)
            %--------------------------------------------------------------------------------------------------%
            % This subroutine tests whether to replace Q by the alternative model, namely the model that
            % minimizes the F-norm of the Hessian subject to the interpolation conditions. It first calculates
            % the alternative model represented by [GALT, PQALT], and sets [GOPT, PQ, HQ] = [GALT, PQALT, 0]
            % if the recent few (three) alternative models are more accurate in predicting the function value of
            % XOPT + D, i.e., if ALL(QALT_BETTER) = TRUE.
            %--------------------------------------------------------------------------------------------------%


            powalg_obj = prima_mat.common.powalg_mod();

            % In-outptuts


            npt = size(xpt, 2);

            %====================%
            % Calculation starts %
            %====================%

            % Establish the alternative model, which is the least Frobenius norm interpolant.
            pqalt = powalg_obj.omega_mul(idz, zmat, fval);
            galt = bmat(:, 1:npt) * fval + powalg_obj.hess_mul(xopt, xpt, pqalt);

            % Replace the current model with the alternative model if ALL(QALT_BETTER) = TRUE, i.e., the
            % recent few alternative models are more accurate in predicting the function value of XOPT + D.
            if all(qalt_better, 'all')
                pq = pqalt;
                hq = zeros(size(hq));
                gopt = galt;
                qalt_better(:) = false;
            end

            %====================%
            %  Calculation ends  %
            %====================%


        end
        function rescon = updateres(~, ximproved, amat, b, delta, dnorm, xopt, rescon)
            %--------------------------------------------------------------------------------------------------%
            % This subroutine updates RESCON when XOPT has been updated by a step D.
            % RESCON holds information about the constraint residuals at the current trust region center XOPT.
            % 1. If if B(J) - AMAT(:, J)^T*XOPT <= DELTA, then RESCON(J) = B(J) - AMAT(:, J)^T*XOPT. Note that
            % RESCON >= 0 in this case, because the algorithm keeps XOPT to be feasible.
            % 2. Otherwise, RESCON(J) is a negative value that B(J) - AMAT(:, J)^T*XOPT >= |RESCON(J)| >= DELTA.
            % RESCON can be updated without calculating the constraints that are far from being active, so that
            % we only need to evaluate the constraints that are nearly active.
            %--------------------------------------------------------------------------------------------------%


            % Norm of D
            % XOPT(N); the updated value of XOPT


            ax = NaN(size(b));

            %====================%
            % Calculation starts %
            %====================%

            % Zaikun 20221115: Currently, UPDATERES does not update RESCON unless XIMPROVED is TRUE. Shouldn't
            % we do it whenever DELTA is updated? Have we MISUNDERSTOOD RESCON?
            if ~ximproved
                return
            end

            mask = (abs(rescon) < dnorm + delta);
            ax(find(mask)) = amat(:, find(mask)).' * xopt;
            mask00 = mask;
            rescon(mask00) = max(b(mask00) - ax(mask00), 0.0);
            mask01 = ~mask00;
            rescon(mask01) = min(-abs(rescon(mask01)) + dnorm, -delta);

            rescon(rescon >= delta) = -rescon(rescon >= delta);

            %%MATLAB:
            %%mask = (abs(rescon) < delta + dnorm);
            %%rescon(mask) = max(b(mask) - (xopt'*amat(:, mask))', 0);
            %%rescon(~mask) = max(rescon(~mask) - dnorm, delta);
            %%rescon(rescon >= delta) = -rescon(rescon >= delta);

            %====================%
            %  Calculation ends  %
            %====================%


        end

    end
end