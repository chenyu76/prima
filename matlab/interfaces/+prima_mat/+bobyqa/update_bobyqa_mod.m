classdef update_bobyqa_mod
    %--------------------------------------------------------------------------------------------------%
    % This module provides subroutines concerning the updates when XPT(:, KNEW) becomes XNEW = XOPT + D.
    %
    % Coded by Zaikun ZHANG (www.zhangzk.net) based on Powell's code and the BOBYQA paper.
    %
    % Dedicated to the late Professor M. J. D. Powell FRS (1936--2015).
    %
    % Started: February 2022
    %
    % Last Modified: Thu 14 Aug 2025 07:34:12 AM CST
    %--------------------------------------------------------------------------------------------------%

    methods
        function [bmat, zmat, info] = updateh(~, knew, kopt, d, xpt, bmat, zmat, varargin)
            %--------------------------------------------------------------------------------------------------%
            % This subroutine updates arrays BMAT and ZMAT in order to replace the interpolation point
            % XPT(:, KNEW) by XNEW = XPT(:, KOPT) + D. See Section 4 of the BOBYQA paper. [BMAT, ZMAT] describes
            % the matrix H in the BOBYQA paper (eq. 2.7), which is the inverse of the coefficient matrix of the
            % KKT system for the least-Frobenius norm interpolation problem: ZMAT holds a factorization of the
            % leading NPT*NPT submatrix OMEGA of H, the factorization being OMEGA = ZMAT*ZMAT^T; BMAT holds the
            % last N ROWs of H except for the (NPT+1)th column. Note that the (NPT + 1)th row and (NPT + 1)th
            % column of H are not stored as they are unnecessary for the calculation.
            %--------------------------------------------------------------------------------------------------%



            linalg_obj = prima_mat.common.linalg_mod();
            powalg_obj = prima_mat.common.powalg_mod();


            % D(N)
            % XPT(N, NPT)


            % BMAT(N, NPT + N)
            % ZMAT(NPT, NPT-N-1)



            hcol = NaN(size(bmat, 2), 1);


            v1 = NaN(size(bmat, 1), 1);
            v2 = NaN(size(bmat, 1), 1);
            vlag = NaN(size(bmat, 2), 1);


            n = size(xpt, 1);
            npt = size(xpt, 2);


            %====================%
            % Calculation starts %
            %====================%

            ipObj = inputParser();
            addParameter(ipObj, 'info', NaN);
            parse(ipObj, varargin{:});
            info = ipObj.Results.info;
            if nargout >= 3
                info = 0;
            end

            % Do anything if KNEW is 0. This can only happen sometimes after a trust-region step.
            if knew <= 0
                % KNEW < 0 is impossible if the input is correct.
                return
            end

            % Put the KNEW-th column of the unupdated H (except for the (NPT+1)th entry) into HCOL. Powell's
            % code does this after ZMAT is rotated below, and then HCOL(1:NPT) = ZMAT(KNEW, 1) * ZMAT(:, 1),
            % which saves flops but also introduces rounding errors due to the rotation.
            hcol(1:npt) = zmat * zmat(knew, :).';
            hcol(npt + 1:npt + n) = bmat(:, knew);

            % Calculate VLAG and BETA and other parameters for (4.9) and (4.14) of the BOBYQA paper.
            beta = powalg_obj.calbeta(kopt, bmat, d, xpt, zmat);
            vlag(:) = powalg_obj.calvlag_lfqint(kopt, bmat, d, xpt, zmat);

            % In theory, DENOM can also be calculated after ZMAT is rotated below. However, this worsened the
            % performance of BOBYQA in a test on 20220413.
            alpha = hcol(knew);
            tau = vlag(knew);
            denom = alpha * beta + tau ^ 2;

            % After the following line, VLAG = H*w - e_KNEW in the NEWUOA paper (where t = KNEW).
            vlag(knew) = vlag(knew) - 1.0;

            % Quite rarely, due to rounding errors, VLAG or BETA may not be finite, or DENOM may not be
            % positive. In such cases, [BMAT, ZMAT] would be destroyed by the update, and hence we would rather
            % not update them at all. Or should we simply terminate the algorithm?
            if ~(isfinite(sum(abs(hcol), 'all') + sum(abs(vlag), 'all') + abs(beta)) && denom > 0)
                if nargout >= 3
                    info = 7;
                end
                return
            end

            % Update the matrix BMAT. It implements the last N rows of (4.9) in the BOBYQA paper.
            v1(:) = (alpha * vlag(npt + 1:npt + n) - tau * hcol(npt + 1:npt + n)) ./ denom;
            v2(:) = (-beta * hcol(npt + 1:npt + n) - tau * vlag(npt + 1:npt + n)) ./ denom;
            bmat(:, :) = bmat + v1 * vlag.' + v2 * hcol.'; %call r2update(bmat, ONE, v1, vlag, ONE, v2, hcol)
            % N.B.: The use of OUTPROD is expensive memory-wise, but it is not our concern in this implementation.
            % Numerically, the update above does not guarantee BMAT(:, NPT+1 : NPT+N) to be symmetric.
            A_slice = linalg_obj.symmetrize(bmat(:, npt + 1:npt + n)); bmat(:, npt + 1:npt + n) = A_slice;

            % Apply Givens rotations to put zeros in the KNEW-th row of ZMAT. After this, ZMAT(KNEW, :) contains
            % only one nonzero at ZMAT(KNEW, 1). Entries of ZMAT are treated as 0 if the moduli are quite small.
            for j = 2:npt - n - 1
                if abs(zmat(knew, j)) > 1.0e-20 * max(abs(zmat), [], 'all')
                    % This threshold is by Powell
                    grot = linalg_obj.planerot(zmat(knew, [1, j]).');
                    zmat(:, [1, j]) = zmat(:, [1, j]) * grot.';
                end
                zmat(knew, j) = 0.0;
            end

            % Complete the updating of ZMAT. See (4.14) of the BOBYQA paper.
            sqrtdn = sqrt(denom);
            zmat(:, 1) = (tau / sqrtdn) * zmat(:, 1) - (zmat(knew, 1) / sqrtdn) * vlag(1:npt);
            % Zaikun 20231012: Either of the following two lines worsens the performance of BOBYQA when the
            % objective function is evaluated with 5 or less correct significance digits. Strange.
            % %zmat(:, 1) = (tau * zmat(:, 1) - zmat(knew, 1) * vlag(1:npt)) / sqrtdn
            % %zmat(knew, 1) = zknew1 / sqrtdn  ! ZKNEW1 is the unupdated ZMAT(KNEW, 1)

            %====================%
            %  Calculation ends  %
            %====================%



        end
        function [kopt, fval, xpt] = updatexf(~, knew, ximproved, f, xnew, kopt, fval, xpt)
            %--------------------------------------------------------------------------------------------------%
            % This subroutine updates [XPT, FVAL, KOPT] so that XPT(:, KNEW) is updated to XNEW.
            %--------------------------------------------------------------------------------------------------%



            % XNEW(N)



            % FVAL(NPT)
            % XPT(N, NPT)



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



        end
        function [gopt, hq, pq] = updateq(~, knew, ximproved, bmat, d, moderr, xdrop, xosav, xpt, zmat, gopt, hq, pq)
            %--------------------------------------------------------------------------------------------------%
            % This subroutine updates GOPT, HQ, and PQ when XPT(:, KNEW) changes from XDROP to XNEW = XOSAV + D,
            % where XOSAV is the unupdated XOPT, namely the XOPT before UPDATEXF is called.
            % See Section 4 of the NEWUOA paper and that of the BOBYQA paper (there is no LINCOA paper).
            % N.B.:
            % XNEW is encoded in [BMAT, ZMAT] after UPDATEH being called, and it also equals XPT(:, KNEW)
            % after UPDATEXF being called. Indeed, we only need BMAT(:, KNEW) instead of the entire matrix.
            %--------------------------------------------------------------------------------------------------%



            linalg_obj = prima_mat.common.linalg_mod();
            powalg_obj = prima_mat.common.powalg_mod();


            % BMAT(N, NPT + N)
            % D(:)

            % XDROP(N)
            % XOSAV(N)
            % XPT(N, NPT)
            % ZMAT(NPT, NPT - N - 1)


            % GOPT(N)
            % HQ(N, N)
            % PQ(NPT)



            pqinc = NaN(numel(pq), 1);


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
            pqinc(:) = moderr * (zmat * zmat(knew, :).'); % pqinc = moderr * omega_col(1_IK, zmat, knew)
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



        end
        function [itest, gopt, hq, pq] = tryqalt(~, bmat, fval, ratio, sl, su, xopt, xpt, zmat, itest, gopt, hq, pq)
            %--------------------------------------------------------------------------------------------------%
            % This subroutine tests whether to replace Q by the alternative model, namely the model that
            % minimizes the F-norm of the Hessian subject to the interpolation conditions. It does the
            % replacement if certain criteria are met (i.e., when ITEST = 3). See the paragraph around (6.12) of
            % the BOBYQA paper.
            %--------------------------------------------------------------------------------------------------%



            powalg_obj = prima_mat.common.powalg_mod();


            % BMAT(N, NPT+N)
            % FVAL(NPT)

            % SL(N)
            % SU(N)
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



            galt = NaN(numel(gopt), 1);

            pgopt = NaN(numel(gopt), 1);
            pqalt = NaN(numel(pq), 1);

            % Debugging variables
            %real(RP) :: intp_tol



            npt = numel(pq);


            %====================%
            % Calculation starts %
            %====================%

            % Calculate the norm square of the projected gradient.
            pgopt(:) = gopt;
            pgopt(xopt >= su) = max(0.0, gopt(xopt >= su));
            pgopt(xopt <= sl) = min(0.0, gopt(xopt <= sl));

            % Calculate the parameters of the least Frobenius norm interpolant to the current data.
            pqalt(:) = zmat * (zmat.' * fval);
            galt(:) = bmat(:, 1:npt) * fval + powalg_obj.hess_mul(xopt, xpt, pqalt);

            % Calculate the norm square of the projected alternative gradient.
            pgalt = galt;
            pgalt(xopt >= su) = max(0.0, galt(xopt >= su));
            pgalt(xopt <= sl) = min(0.0, galt(xopt <= sl));

            % Test whether to replace the new quadratic model by the least Frobenius norm interpolant,
            % making the replacement if the test is satisfied.
            % N.B.: In the following IF, Powell's condition does not check RATIO. The condition here (with RATIO
            % > TENTH)is adopted and adapted from NEWUOA, and it seems to improve the performance.
            % %if (inprod(pgopt, pgopt) < TEN * inprod(pgalt, pgalt)) then  ! Powell's code
            if ratio > 0.1 || sum(pgopt .* pgopt, 'all') < 10.0 * sum(pgalt .* pgalt, 'all')
                itest = 0;
            else
                itest = itest + 1;
            end
            if itest >= 3
                gopt(:) = galt;
                pq(:) = pqalt;
                hq = zeros(size(hq));
                itest = 0;
            end

            %====================%
            %  Calculation ends  %
            %====================%



        end

    end
end