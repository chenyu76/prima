classdef update_newuoa_mod
    %--------------------------------------------------------------------------------------------------%
    % This module provides subroutines concerning the updates when XPT(:, KNEW) becomes XNEW = XOPT + D.
    %
    % Coded by Zaikun ZHANG (www.zhangzk.net) based on Powell's code and the NEWUOA paper.
    %
    % Dedicated to the late Professor M. J. D. Powell FRS (1936--2015).
    %
    % Started: July 2020
    %
    % Last Modified: Friday, March 15, 2024 PM03:38:22
    %--------------------------------------------------------------------------------------------------%

    methods
        function [kopt, fval, xpt] = updatexf(~, knew, ximproved, f, xnew, kopt, fval, xpt)
            %--------------------------------------------------------------------------------------------------%
            % This subroutine updates [XPT, FVAL, KOPT] so that XPT(:, KNEW) is updated to XNEW.
            %--------------------------------------------------------------------------------------------------%

            % Common modules
            consts_obj = prima_mat.common.consts_mod();
            debug_obj = prima_mat.common.debug_mod();
            infnan_obj = prima_mat.common.infnan_mod();


            % Inputs


            % XNEW(N)

            % In-outputs


            % FVAL(NPT)
            % XPT(N, NPT)

            % Local variables
            srname = "UPDATEXF";


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
                debug_obj.assert(numel(fval) == npt && ~any(infnan_obj.is_nan_sp(fval) | infnan_obj.is_posinf(fval), 'all'), "SIZE(FVAL) == NPT and FVAL is not NaN or +Inf", srname);
                debug_obj.assert(~any(fval < fval(kopt), 'all'), "FVAL(KOPT) = MINVAL(FVAL)", srname);
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

            % KOPT is NOT identical to MINLOC(FVAL). Indeed, if FVAL(KNEW) = FVAL(KOPT) and KNEW < KOPT, then
            % MINLOC(FVAL) = KNEW /= KOPT. Do not change KOPT in this case.
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
                debug_obj.assert(~any(fval < fval(kopt), 'all'), "FVAL(KOPT) = MINVAL(FVAL)", srname);
            end

        end
        function [gopt, hq, pq] = updateq(~, idz, knew, ximproved, bmat, d, moderr, xdrop, xosav, xpt, zmat, gopt, hq, pq)
            %--------------------------------------------------------------------------------------------------%
            % This subroutine updates GOPT, HQ, and PQ when XPT(:, KNEW) changes from XDROP to XNEW = XOSAV + D,
            % where XOSAV is the unupdated XOPT, namely the XOPT before UPDATEXF is called.
            % See Section 4 of the NEWUOA paper and that of the BOBYQA paper (there is no LINCOA paper).
            % N.B.:
            % 1. XNEW is encoded in [BMAT, ZMAT, IDZ] after UPDATEH being called, and it also equals XPT(:, KNEW)
            % after UPDATEXF being called.
            % 2. Indeed, we only need BMAT(:, KNEW) instead of the entire matrix.
            % 3. In Powell's implementation of NEWUOA, the quadratic model is represented by [GQ, PQ, HQ], where
            % GQ is the gradient of the quadratic model at XBASE. However, Powell implemented BOBYQA and LINCOA
            % without GQ but with GOPT, which is the gradient at XBASE + XOPT. In our implementation, we also
            % use GOPT instead of GQ.
            %--------------------------------------------------------------------------------------------------%

            % Common modules
            consts_obj = prima_mat.common.consts_mod();
            debug_obj = prima_mat.common.debug_mod();
            infnan_obj = prima_mat.common.infnan_mod();
            linalg_obj = prima_mat.common.linalg_mod();
            powalg_obj = prima_mat.common.powalg_mod();


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
        function [itest, gopt, hq, pq] = tryqalt(~, idz, bmat, fval, ratio, xopt, xpt, zmat, itest, gopt, hq, pq)
            %--------------------------------------------------------------------------------------------------%
            % This subroutine tests whether to replace Q by the alternative model, namely the model that
            % minimizes the F-norm of the Hessian subject to the interpolation conditions. It does the
            % replacement if certain criteria are met (i.e., when ITEST = 3). See the paragraph around (8.4) of
            % the NEWUOA paper.
            %--------------------------------------------------------------------------------------------------%

            % Common modules
            consts_obj = prima_mat.common.consts_mod();
            debug_obj = prima_mat.common.debug_mod();
            infnan_obj = prima_mat.common.infnan_mod();
            linalg_obj = prima_mat.common.linalg_mod();
            powalg_obj = prima_mat.common.powalg_mod();


            % Inputs

            % BMAT(N, NPT+N)
            % FVAL(NPT)

            % XOPT(N)
            % XOPT(N, NPT)
            % ZMAT(NPT, NPT-N-1)

            % In-output

            % GOPT(N)
            % HQ(N, N)
            % PQ(NPT)
            % N.B.:
            % GOPT, HQ, and PQ should be INTENT(INOUT) instead of INTENT(OUT). According to the Fortran 2018
            % standard, an INTENT(OUT) dummy argument becomes undefined on invocation of the procedure.
            % Therefore, if the procedure does not define such an argument, its value becomes undefined,
            % which is the case for HQ and PQ when ITEST < 3 at exit. In addition, the information in GOPT is
            % needed for defining ITEST, so it must be INTENT(INOUT).

            % Local variables
            srname = "TRYQALT";


            galt = NaN(numel(gopt), 1);
            pqalt = NaN(numel(pq), 1);

            % Sizes
            n = fix(numel(gopt));
            npt = fix(numel(pq));

            % Preconditions
            if consts_obj.DEBUGGING
                debug_obj.assert(n >= 1 && npt >= n + 2, "N >= 1, NPT >= N + 2", srname);
                % By the definition of RATIO in ratio.f90, RATIO cannot be NaN unless the actual reduction is
                % NaN, which should NOT happen due to the moderated extreme barrier.
                debug_obj.assert(~infnan_obj.is_nan_sp(ratio), "RATIO is not NaN", srname);
                debug_obj.assert(numel(fval) == npt && ~any(infnan_obj.is_nan_sp(fval) | infnan_obj.is_posinf(fval), 'all'), "SIZE(FVAL) == NPT and FVAL is not NaN or +Inf", srname);
                debug_obj.assert(size(bmat, 1) == n && size(bmat, 2) == npt + n, "SIZE(BMAT)==[N, NPT+N]", srname);
                debug_obj.assert(linalg_obj.issymmetric(bmat(:, npt + 1:npt + n)), "BMAT(:, NPT+1:NPT+N) is symmetric", srname);
                debug_obj.assert(size(zmat, 1) == npt && size(zmat, 2) == npt - n - 1, "SIZE(ZMAT) == [NPT, NPT - N - 1]", srname);
                debug_obj.assert(numel(gopt) == n, "SIZE(GOPT) = N", srname);
                debug_obj.assert(size(hq, 1) == n && linalg_obj.issymmetric(hq), "HQ is an NxN symmetric matrix", srname);
                debug_obj.assert(numel(pq) == npt, "SIZE(PQ) = NPT", srname);
            end

            %====================%
            % Calculation starts %
            %====================%

            % Calculate the parameters of the least Frobenius norm interpolant to the current data.
            pqalt(:) = powalg_obj.omega_mul(idz, zmat, fval);
            galt(:) = linalg_obj.matprod21(bmat(:, 1:npt), fval) + powalg_obj.hess_mul(xopt, xpt, pqalt);

            % Test whether to replace the new quadratic model by the least Frobenius norm interpolant, making
            % the replacement if the test is satisfied. In the sequel, TEN seems to work a bit better than 100.
            % In addition, Powell checked the magnitude of ABS(RATIO) instead of RATIO.
            % %if (abs(ratio) > 0.01 .or. inprod(gopt, gopt) < 1.0E2_RP * inprod(galt, galt)) then ! Powell's code
            if ratio > consts_obj.TENTH || linalg_obj.inprod(gopt, gopt) < consts_obj.TEN * linalg_obj.inprod(galt, galt)
                itest = 0;
            else
                itest = itest + 1;
            end
            if itest >= 3
                gopt(:) = galt;
                pq(:) = pqalt;
                hq = repmat(consts_obj.ZERO, size(hq));
                itest = 0;
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

    end
end