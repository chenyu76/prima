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

            % Common modules
            consts_obj = consts_mod();
            debug_obj = debug_mod();
            infnan_obj = infnan_mod();


            % Inputs


            % XNEW(N)

            % In-outputs


            % FVAL(NPT)
            % XPT(N, NPT)

            % Local variables
            srname = "UPDATEXF";
            n = NaN;
            npt = NaN;

            % Sizes
            n = size(xpt, 1);
            npt = size(xpt, 2);

            % Preconditions
            if consts_obj.DEBUGGING
                debug_obj.assert(n >= 1 && npt >= n + 2, "N >= 1, NPT >= N + 2", srname);
                debug_obj.assert(knew >= 0 && knew <= npt, "0 <= KNEW <= NPT", srname);
                debug_obj.assert(kopt >= 1 && kopt <= npt, "1 <= KOPT <= NPT", srname);
                debug_obj.assert(knew >= 1 || ~ximproved, "KNEW >= 1 unless X is not improved", srname);
                debug_obj.assert(knew ~= kopt || ximproved, "KNEW /= KOPT unless X is improved", srname);
                debug_obj.assert(numel(xnew) == n && all(infnan_obj.is_finite(xnew), 'all'), "SIZE(XNEW) == N, XNEW is finite", srname);
                debug_obj.assert(~(infnan_obj.is_nan_sp(f) || infnan_obj.is_posinf(f)), "F is not NaN or +Inf", srname);
                debug_obj.assert(all(infnan_obj.is_finite(xpt), 'all'), "XPT is finite", srname);
                debug_obj.assert(numel(fval) == npt && ~any(infnan_obj.is_nan(fval) | infnan_obj.is_posinf(fval), 'all'), "SIZE(FVAL) == NPT and FVAL is not NaN or +Inf", srname);
            end

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

            % Postconditions
            if consts_obj.DEBUGGING
                debug_obj.assert(size(xpt, 1) == n && size(xpt, 2) == npt && all(infnan_obj.is_finite(xpt), 'all'), "SIZE(XPT) == [N, NPT], XPT is finite", srname);
                debug_obj.assert(kopt >= 1 && kopt <= npt, "1 <= KOPT <= NPT", srname);
            end

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

            % Common modules
            consts_obj = consts_mod();
            debug_obj = debug_mod();
            infnan_obj = infnan_mod();
            linalg_obj = linalg_mod();
            powalg_obj = powalg_mod();


            % Inputs



            % BMAT(N, NPT + N)
            % D(:)

            % XDROP(N)
            % XOSAV(N)
            % XPT(N, NPT)
            % ZMAT(NPT, NPT - N - 1)

            % In-outputs
            % GOPT(N)
            % HQ(N, N)
            % PQ(NPT)

            % Local variables
            srname = "UPDATEQ";
            n = NaN;
            npt = NaN;
            pqinc = NaN(numel(pq), 1);

            % Sizes
            n = fix(numel(gopt));
            npt = fix(numel(pq));

            % Preconditions
            if consts_obj.DEBUGGING
                debug_obj.assert(n >= 1 && npt >= n + 2, "N >= 1, NPT >= N + 2", srname);
                debug_obj.assert(idz >= 1 && idz <= size(zmat, 2) + 1, "1 <= IDZ <= SIZE(ZMAT, 2) + 1", srname);
                debug_obj.assert(knew >= 0 && knew <= npt, "0 <= KNEW <= NPT", srname);
                debug_obj.assert(knew >= 1 || ~ximproved, "KNEW >= 1 unless X is not improved", srname);
                debug_obj.assert(numel(xdrop) == n && all(infnan_obj.is_finite(xdrop), 'all'), "SIZE(XDROP) == N, XDROP is finite", srname);
                debug_obj.assert(numel(xosav) == n && all(infnan_obj.is_finite(xosav), 'all'), "SIZE(XOSAV) == N, XOSAV is finite", srname);
                debug_obj.assert(size(bmat, 1) == n && size(bmat, 2) == npt + n, "SIZE(BMAT)==[N, NPT+N]", srname);
                debug_obj.assert(linalg_obj.issymmetric(bmat(:, npt + 1:npt + n)), "BMAT(:, NPT+1:NPT+N) is symmetric", srname);
                debug_obj.assert(size(zmat, 1) == npt && size(zmat, 2) == npt - n - 1, "SIZE(ZMAT) == [NPT, NPT - N - 1]", srname);
                debug_obj.assert(all(infnan_obj.is_finite(xpt), 'all'), "XPT is finite", srname);
                debug_obj.assert(size(hq, 1) == n && linalg_obj.issymmetric(hq), "HQ is an NxN symmetric matrix", srname);
            end

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
            pq(knew) = consts_obj.ZERO;

            % Update the implicit part of the Hessian.
            pqinc(:) = moderr * powalg_obj.omega_col(idz, zmat, knew);
            pq(:) = pq + pqinc;

            % Update the gradient, which needs the updated XPT.
            gopt(:) = gopt + moderr * bmat(:, knew) + powalg_obj.hess_mul(xosav, xpt, pqinc);

            % Further update GOPT if XIMPROVED is TRUE, as XOPT changes from XOSAV to XNEW = XOSAV + D.
            if ximproved
                gopt(:) = gopt + powalg_obj.hess_mul(d, xpt, pq, 'hq', hq);
            end

            %====================%
            %  Calculation ends  %
            %====================%

            % Postconditions
            if consts_obj.DEBUGGING
                debug_obj.assert(numel(gopt) == n, "SIZE(GOPT) = N", srname);
                debug_obj.assert(size(hq, 1) == n && linalg_obj.issymmetric(hq), "HQ is an NxN symmetric matrix", srname);
                debug_obj.assert(numel(pq) == npt, "SIZE(PQ) = NPT", srname);
            end

        end
        function [qalt_better, gopt, pq, hq, galt, pqalt] = tryqalt(~, idz, bmat, fval, xopt, xpt, zmat, qalt_better, gopt, pq, hq, galt, pqalt)
            %--------------------------------------------------------------------------------------------------%
            % This subroutine tests whether to replace Q by the alternative model, namely the model that
            % minimizes the F-norm of the Hessian subject to the interpolation conditions. It first calculates
            % the alternative model represented by [GALT, PQALT], and sets [GOPT, PQ, HQ] = [GALT, PQALT, 0]
            % if the recent few (three) alternative models are more accurate in predicting the function value of
            % XOPT + D, i.e., if ALL(QALT_BETTER) = TRUE.
            %--------------------------------------------------------------------------------------------------%
            % Common modules
            consts_obj = consts_mod();
            debug_obj = debug_mod();
            infnan_obj = infnan_mod();
            linalg_obj = linalg_mod();
            powalg_obj = powalg_mod();


            % Inputs

            % BMAT(N, NPT + N)
            % FVAL(NPT)
            % XOPT(N)
            % XPT(N, NPT)
            % ZMAT(NPT, NPT - N - 1)

            % In-outptuts
            % QALT_BETTER(3)
            % GOPT(N)
            % PQ(NPT)
            % HQ(N, N)

            % Outputs
            % GALT(N)
            % PQALT(NPT)

            % Local variables
            srname = "TRYQALT";
            n = NaN;
            npt = NaN;

            % Sizes
            n = size(xpt, 1);
            npt = size(xpt, 2);

            % Preconditions
            if consts_obj.DEBUGGING
                debug_obj.assert(n >= 1 && npt >= n + 2, "N >= 1, NPT >= N + 2", srname);
                debug_obj.assert(idz >= 1 && idz <= size(zmat, 2) + 1, "1 <= IDZ <= SIZE(ZMAT, 2) + 1", srname);
                debug_obj.assert(numel(xopt) == n && all(infnan_obj.is_finite(xopt), 'all'), "SIZE(XOPT) == N, XOPT is finite", srname);
                debug_obj.assert(all(infnan_obj.is_finite(xpt), 'all'), "XPT is finite", srname);
                debug_obj.assert(numel(fval) == npt && ~any(infnan_obj.is_nan(fval) | infnan_obj.is_posinf(fval), 'all'), "SIZE(FVAL) == NPT and FVAL is not NaN or +Inf", srname);
                debug_obj.assert(size(bmat, 1) == n && size(bmat, 2) == npt + n, "SIZE(BMAT)==[N, NPT+N]", srname);
                debug_obj.assert(linalg_obj.issymmetric(bmat(:, npt + 1:npt + n)), "BMAT(:, NPT+1:NPT+N) is symmetric", srname);
                debug_obj.assert(size(zmat, 1) == npt && size(zmat, 2) == npt - n - 1, "SIZE(ZMAT) == [NPT, NPT - N - 1]", srname);
                debug_obj.assert(numel(gopt) == n, "SIZE(GOPT) = N", srname);
                debug_obj.assert(size(hq, 1) == n && linalg_obj.issymmetric(hq), "HQ is an NxN symmetric matrix", srname);
                debug_obj.assert(numel(pq) == npt, "SIZE(PQ) = NPT", srname);
                debug_obj.assert(numel(galt) == n, "SIZE(GALT) = N", srname);
                debug_obj.assert(numel(pqalt) == npt, "SIZE(PQALT) = NPT", srname);
            end

            %====================%
            % Calculation starts %
            %====================%

            % Establish the alternative model, which is the least Frobenius norm interpolant.
            pqalt(:) = powalg_obj.omega_mul(idz, zmat, fval);
            galt(:) = linalg_obj.matprod21(bmat(:, 1:npt), fval) + powalg_obj.hess_mul(xopt, xpt, pqalt);

            % Replace the current model with the alternative model if ALL(QALT_BETTER) = TRUE, i.e., the
            % recent few alternative models are more accurate in predicting the function value of XOPT + D.
            if all(qalt_better, 'all')
                pq(:) = pqalt;
                hq = repmat(consts_obj.ZERO, size(hq));
                gopt(:) = galt;
                qalt_better = repmat(false, size(qalt_better));
            end

            %====================%
            %  Calculation ends  %
            %====================%

            % Postconditions
            if consts_obj.DEBUGGING
                debug_obj.assert(numel(gopt) == n, "SIZE(GOPT) = N", srname);
                debug_obj.assert(size(hq, 1) == n && linalg_obj.issymmetric(hq), "HQ is an NxN symmetric matrix", srname);
                debug_obj.assert(numel(pq) == npt, "SIZE(PQ) = NPT", srname);
                debug_obj.assert(numel(galt) == n, "SIZE(GALT) = N", srname);
                debug_obj.assert(numel(pqalt) == npt, "SIZE(PQALT) = NPT", srname);
            end

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

            % Common modules
            consts_obj = consts_mod();
            debug_obj = debug_mod();
            infnan_obj = infnan_mod();
            linalg_obj = linalg_mod();


            % Inputs

            % AMAT(N, M)
            % B(M)

            % Norm of D
            % XOPT(N); the updated value of XOPT

            % In-outputs
            % RESCON(M)

            % Local variables
            srname = "UPDATERES";
            m = NaN;
            n = NaN;
            mask = false(numel(b), 1);
            ax = NaN(numel(b), 1);

            % Sizes
            m = fix(numel(b));
            n = fix(numel(xopt));

            % Preconditions
            if consts_obj.DEBUGGING
                debug_obj.assert(size(amat, 1) == n && size(amat, 2) == m, "SIZE(AMAT) == [N, M]", srname);
                debug_obj.assert(delta > 0, "DELTA > 0", srname);
                debug_obj.assert(dnorm > 0, "DNORM > 0", srname);
                debug_obj.assert(all(infnan_obj.is_finite(xopt), 'all'), "XOPT is finite", srname);
                debug_obj.assert(numel(rescon) == m, "SIZE(RESCON) == M", srname);
                % Zaikun 20221115: The following cannot pass?! Is it due to the update of DELTA? Did we
                % misunderstand Powell's definition of RESCON?
                %call assert(all((rescon >= 0 .and. rescon <= delta) .or. rescon <= -delta), &
                %    & '0 <= RESCON <= DELTA or RESCON <= -DELTA', srname)

            end

            %====================%
            % Calculation starts %
            %====================%

            % Zaikun 20221115: Currently, UPDATERES does not update RESCON unless XIMPROVED is TRUE. Shouldn't
            % we do it whenever DELTA is updated? Have we MISUNDERSTOOD RESCON?
            if ~ximproved
                return
            end

            mask(:) = (abs(rescon) < dnorm + delta);
            ax(linalg_obj.trueloc(mask)) = linalg_obj.matprod12(xopt, amat(:, linalg_obj.trueloc(mask)));
            mask00 = mask; %Unsupported statement inside WHERE block: StmtLineBreak 1
            rescon(mask00) = max(b(mask00) - ax(mask00), consts_obj.ZERO); %Unsupported statement inside WHERE block: StmtLineBreak 1
            mask01 = ~mask00; %Unsupported statement inside WHERE block: StmtLineBreak 1
            rescon(mask01) = min(-abs(rescon(mask01)) + dnorm, -delta); %Unsupported statement inside WHERE block: StmtLineBreak 1

            rescon(linalg_obj.trueloc(rescon >= delta)) = -rescon(linalg_obj.trueloc(rescon >= delta));

            %%MATLAB:
            %%mask = (abs(rescon) < delta + dnorm);
            %%rescon(mask) = max(b(mask) - (xopt'*amat(:, mask))', 0);
            %%rescon(~mask) = max(rescon(~mask) - dnorm, delta);
            %%rescon(rescon >= delta) = -rescon(rescon >= delta);

            %====================%
            %  Calculation ends  %
            %====================%

            % Postconditions
            if consts_obj.DEBUGGING
                debug_obj.assert(numel(rescon) == m, "SIZE(RESCON) == M", srname);
                %call assert(all((rescon >= 0 .and. rescon <= delta) .or. rescon <= -delta), &
                %    & '0 <= RESCON <= DELTA or RESCON <= -DELTA', srname)

            end
        end

    end
end