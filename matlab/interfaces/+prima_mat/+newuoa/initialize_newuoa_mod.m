classdef initialize_newuoa_mod
    %--------------------------------------------------------------------------------------------------%
    % This module performs the initialization of NEWUOA, described in Section 3 of the NEWUOA paper.
    %
    % Coded by Zaikun ZHANG (www.zhangzk.net) based on Powell's code and the NEWUOA paper.
    %
    % Started: July 2020
    %
    % Dedicated to the late Professor M. J. D. Powell FRS (1936--2015).
    %
    % Last Modified: Tue 10 Feb 2026 01:56:05 PM CET
    %--------------------------------------------------------------------------------------------------%

    methods
        function [ij, kopt, nf, fhist, fval, xbase, xhist, xpt, info] = initxf(~, calfun, iprint, maxfun, ftarget, rhobeg, x0, ij, fhist, fval, xbase, xhist, xpt)
            %--------------------------------------------------------------------------------------------------%
            % This subroutine does the initialization about the interpolation points & their function values.
            %
            % N.B.:
            % 1. Remark on IJ:
            % If NPT <= 2*N + 1, then IJ is empty. Assume that NPT >= 2*N + 2. Then SIZE(IJ) = [2, NPT-2*N-1].
            % IJ contains integers between 1 and 2*N. For each K > 2*N + 1, XPT(:, K) is
            % XPT(:, IJ(1, K) + 1) + XPT(:, IJ(2, K) + 1). The 1 in IJ + 1 comes from the fact that XPT(:, 1)
            % corresponds to the base point XBASE. Let I = IJ(1, K) if such a number is <= N; otherwise, let
            % I = IJ(1, K) - N; define J by IJ(2, K) in a similar fashion. Then all the entries of XPT(:, K)
            % are zero except that the I and J entries are RHOBEG or -RHOBEG. Indeed, XPT(I, K) is RHOBEG if
            % IJ(1, K) <= N  and -RHOBEG otherwise; XPT(J, K) is similar. Consequently, the Hessian of the
            % quadratic model will get a possibly nonzero (I, J) entry. In the code, IJ is defined according to
            % Powell's original code as well as Section 3 of the NEWUOA paper and (2.4) of the BOBYQA paper.
            % 2. At return,
            % INFO = INFO_DFT: initialization finishes normally
            % INFO = FTARGET_ACHIEVED: return because F <= FTARGET
            % INFO = NAN_INF_X: return because X contains NaN
            % INFO = NAN_INF_F: return because F is either NaN or +Inf
            %--------------------------------------------------------------------------------------------------%


            checkexit_obj = prima_mat.common.checkexit_mod();

            evaluate_obj = prima_mat.common.evaluate_mod();
            history_obj = prima_mat.common.history_mod();

            message_obj = prima_mat.common.message_mod();

            powalg_obj = prima_mat.common.powalg_mod();

            % X0(N)


            % IJ(2, MAX(0_IK, NPT-2*N-1_IK))


            % FHIST(MAXFHIST)
            % FVAL(NPT)
            % XBASE(N)
            % XHIST(N, MAXXHIST)
            % XPT(N, NPT)


            solver = "NEWUOA";

            evaluated = false(size(fval));

            n = size(xpt, 1);
            npt = size(xpt, 2);

            %====================%
            % Calculation starts %
            %====================%

            % Initialize INFO to the default value. At return, an INFO different from this value will indicate
            % an abnormal return.
            info = 0;

            % Initialize XBASE to X0.
            xbase = x0;

            % EVALUATED is a boolean array with EVALUATED(I) indicating whether the function value of the I-th
            % interpolation point has been evaluated. We need it for a portable counting of the number of
            % function evaluations, especially if the loop is conducted asynchronously. However, the loop here
            % is not fully parallelizable if NPT>2N+1, as the definition XPT(;, 2N+2:end) involves FVAL(1:2N+1).
            evaluated(:) = false;

            % Initialize XHIST, FHIST, and FVAL. Otherwise, compilers may complain that they are not
            % (completely) initialized if the initialization aborts due to abnormality (see CHECKEXIT).
            % N.B.: 1. Initializing them to NaN would be more reasonable (NaN is not available in Fortran).
            % 2. Do not initialize the models if the current initialization aborts due to abnormality. Otherwise,
            % errors or exceptions may occur, as FVAL and XPT etc are uninitialized.
            xhist = repmat(-realmax, size(xhist));
            fhist(:) = realmax;
            fval(:) = realmax;

            % Initialize XPT(:, 1: MIN(2*N + 1, NPT)).
            xpt(:, 1) = 0.0;
            xpt(:, 2:n + 1) = rhobeg * eye(n);
            % After the following line, XPT(:, 2*N+2 : NPT) = ZERO if it is nonempty. It will be revised later
            % according to FVAL(2 : 2*N + 1).
            xpt(:, n + 2:npt) = -rhobeg * eye(n, npt - n - 1);

            % Set FVAL(1 : min(2*N + 1, NPT)) by evaluating F. Totally parallelizable except for FMSG.
            for k = 1:min(npt, 2 * n + 1)
                x = xpt(:, k) + xbase;
                f = evaluate_obj.evaluatef(calfun, x);

                % Print a message about the function evaluation according to IPRINT.
                message_obj.fmsg(solver, "Initialization", iprint, k, rhobeg, f, x);
                % Save X and F into the history.
                [xhist, fhist] = history_obj.savehist(k, x, xhist, f, fhist);

                evaluated(k) = true;
                fval(k) = f;

                % Check whether to exit.
                subinfo = checkexit_obj.checkexit_unc(maxfun, k, f, ftarget, x);
                if subinfo ~= 0
                    info = subinfo;
                    break
                end
            end

            % Set IJ.
            % In general, when NPT = (N+1)*(N+2)/2, we can set IJ(:, 1 : NPT - (2*N+1)) to ANY permutation
            % of {{I, J} : 1 <= I /= J <= N}; when NPT < (N+1)*(N+2)/2, we can set it to the first NPT - (2*N+1)
            % elements of such a permutation. The following IJ is defined according to Powell's code. See also
            % Section 3 of the NEWUOA paper and (2.4) of the BOBYQA paper.
            ij = powalg_obj.setij(n, npt);

            % Further revise IJ according to FVAL(2 : 2*N + 1).
            % N.B.:
            % 1. For each K below, the following lines revises IJ(:, K) as follows:
            % change IJ(1, K) to IJ(1, K) + N if FVAL(IJ(1, K) + N + 1) < FVAL(IJ(1, K) + 1);
            % change IJ(2, K) to IJ(2, K) + N if FVAL(IJ(2, K) + N + 1) < FVAL(IJ(2, K) + 1).
            % The 1 in IJ + 1 comes from the fact that XPT(:, 1) corresponds to XBASE.
            % 2. The idea of this revision is as follows: Let [I, J] = IJ(:, K) with the IJ BEFORE the revision;
            % XPT(:, K) is the sum of either {XPT(:, I+1) or XPT(:, I+N+1)} + {XPT(:, J+1) or XPT(:, J+N+1)},
            % each choice being made in favor of the point that has a lower function value, with the hope that
            % such a choice will more likely render an XPT(:, K) with a lower function value.
            % 3. This revision is OPTIONAL. Due to this revision, the definition of XPT(:, 2*N + 2 : NPT) relies
            % on FVAL(2 : 2*N + 1), and it is the sole origin of the such dependency. If we remove the revision
            % IJ, then the evaluations of FVAL(1 : NPT) can be merged, and they are totally PARALLELIZABLE; this
            % can be beneficial if the function evaluations are expensive, which is likely the case.
            ij(1, fval(ij(1, :) + n + 1) < fval(ij(1, :) + 1)) = ij(1, fval(ij(1, :) + n + 1) < fval(ij(1, :) + 1)) + n;
            ij(2, fval(ij(2, :) + n + 1) < fval(ij(2, :) + 1)) = ij(2, fval(ij(2, :) + n + 1) < fval(ij(2, :) + 1)) + n;
            % MATLAB (but not Fortran) can index a vector by a 2D array of indices, thus the MATLAB code is
            %%MATLAB: ij(fval(ij + n + 1) < fval(ij + 1)) = ij(fval(ij + n  + 1) < fval(ij + 1)) + n;

            % Set XPT(:, 2*N + 2 : NPT). It depends on IJ and hence on FVAL(2 : 2*N + 1). Indeed, XPT(:, K) has
            % only two nonzeros for each K >= 2*N+2.
            xpt(:, 2 * n + 2:npt) = xpt(:, ij(1, :) + 1) + xpt(:, ij(2, :) + 1);

            % Set FVAL(2*N + 2 : NPT) by evaluating F. Totally parallelizable except for FMSG.
            if info == 0
                for k = 2 * n + 2:npt
                    x = xpt(:, k) + xbase;
                    f = evaluate_obj.evaluatef(calfun, x);

                    % Print a message about the function evaluation according to IPRINT.
                    message_obj.fmsg(solver, "Initialization", iprint, k, rhobeg, f, x);
                    % Save X and F into the history.
                    [xhist, fhist] = history_obj.savehist(k, x, xhist, f, fhist);

                    evaluated(k) = true;
                    fval(k) = f;

                    % Check whether to exit.
                    subinfo = checkexit_obj.checkexit_unc(maxfun, k, f, ftarget, x);
                    if subinfo ~= 0
                        info = subinfo;
                        break
                    end
                end
            end

            % Set NF, KOPT
            nf = nnz(evaluated); %%MATLAB: nf = sum(evaluated);
            kopt = fortran.minloc(fval, 'mask', evaluated, 'dim', 1);
            %%MATLAB: fopt = min(fval(evaluated)); kopt = find(evaluated & ~(fval > fopt), 1, 'first')

            %====================%
            %  Calculation ends  %
            %====================%


        end
        function [gopt, hq, pq, info] = initq(~, ij, fval, xpt, gopt, hq, pq, varargin)
            %--------------------------------------------------------------------------------------------------%
            % This subroutine initializes the quadratic model represented by [GOPT, HQ, PQ] so that its gradient
            % at XBASE + XPT(:,KOPT) is GOPT; its Hessian is HQ + sum_{K=1}^NPT PQ(K)*XPT(:, K)*XPT(:, K)'.
            %--------------------------------------------------------------------------------------------------%


            % IJ(2, MAX(0_IK, NPT - 2_IK * N - 1_IK))
            % FVAL(NPT)
            % XPT(N, NPT)


            % GOPT(N)
            % HQ(N, N)
            % PQ(NPT)


            n = size(xpt, 1);
            npt = size(xpt, 2);

            %====================%
            % Calculation starts %
            %====================%

            rhobeg = max(abs(xpt(:, 2)), [], 'all'); % Read RHOBEG from XPT.
            fbase = fval(1); % FBASE is the function value at XBASE.

            % Set GOPT by the forward difference.
            gopt(1:n) = (fval(2:n + 1) - fbase) ./ rhobeg;

            % The interpolation conditions decide GOPT(1:NDIAG) and the first NDIAG diagonal 2nd derivatives of
            % the initial quadratic model by a quadratic interpolation on three points, which is equivalent to
            % the central finite difference.
            ndiag = min(npt - n - 1, n);

            % Revise GOPT(1:NDIAG) to the value provided by the central finite difference.
            gopt(1:ndiag) = 0.5 * (gopt(1:ndiag) + (fbase - fval(n + 2:n + 1 + ndiag)) ./ rhobeg);

            % Set the diagonal of HQ by the 2nd-order central finite difference. If we do this before the
            % revision of GOPT(1:NDIAG), we can avoid the calculation of FVAL(K + 1) - FBASE) / RHOBEG. But we
            % prefer to decouple the initialization of GOPT and HQ. We are not concerned by this amount of flops.
            hq = zeros(size(hq));
            for k = 1:ndiag
                hq(k, k) = ((fval(k + 1) - fbase) / rhobeg - (fbase - fval(k + n + 1)) / rhobeg) / rhobeg;
            end
            %%MATLAB:
            %%hdiag = ((fval(2 : ndiag+1) - fbase) / rhobeg - (fbase - fval(n+2 : n+ndiag+1)) / rhobeg) / rhobeg
            %%hq(1:ndiag, 1:ndiag) = diag(hdiag)

            % When NPT > 2*N + 1, set the off-diagonal entries of HQ.
            for k = 1:npt - 2 * n - 1
                % With the I, J, XI, and XJ defined below, we have
                % FVAL(K+2*n+1) = F(XBASE + XI*e_I + XJ*e_J),
                % FVAL(IJ(1, K) + 1) = F(XBASE + XI*e_I),
                % FVAL(IJ(2, K) + 1) = F(XBASE + XJ*e_J).
                % The 1 in IJ + 1 comes from the fact that XPT(:, 1) corresponds to XBASE.
                % Thus the HQ(I,J) defined below approximates frac{partial^2}{partial X_I partial X_J} F(XBASE).
                % N.B.: Here, exchanging I and J will not lead to any change in precise arithmetic. Powell's
                % code exchanges I and J if needed to ensure that  I > J. This is because Powell's code saves HQ
                % as a 1D array that contains the lower triangular part of this symmetric matrix.
                i = mod(ij(1, k) - 1, n) + 1;
                j = mod(ij(2, k) - 1, n) + 1;
                xi = xpt(i, k + 2 * n + 1);
                xj = xpt(j, k + 2 * n + 1);
                hq(i, j) = (fbase - fval(ij(1, k) + 1) - fval(ij(2, k) + 1) + fval(k + 2 * n + 1)) / (xi * xj);
                hq(j, i) = hq(i, j);
            end

            [~, kopt] = min(fval);
            if kopt ~= 1
                gopt = gopt + hq * xpt(:, kopt);
            end

            pq(:) = 0.0;

            ipObj = inputParser();
            addParameter(ipObj, 'info', NaN);
            parse(ipObj, varargin{:});
            info = ipObj.Results.info;
            if nargout >= 4
                if any(isnan(gopt), 'all') || any(isnan(hq), 'all')
                    info = -3;
                else
                    info = 0;
                end
            end

            %====================%
            %  Calculation ends  %
            %====================%


        end
        function [idz, bmat, zmat, info] = inith(~, ij, xpt, bmat, zmat, varargin)
            %--------------------------------------------------------------------------------------------------%
            % This subroutine initializes [IDZ, BMAT, ZMAT] which represents the matrix H in (3.12) of the
            % NEWUOA paper (see also (2.7) of the BOBYQA paper).
            %--------------------------------------------------------------------------------------------------%


            %use, non_intrinsic :: powalg_mod, only : errh


            % IJ(2, MAX(0_IK, NPT - 2_IK * N - 1_IK))
            % XPT(N, NPT)
            % N.B.: XPT is essentially only used for debugging, to test the error in the initial H. The initial
            % ZMAT and BMAT are completely defined by RHOBEG and IJ.


            % BMAT(N, NPT + N)
            % ZMAT(NPT, NPT - N - 1)


            n = size(xpt, 1);
            npt = size(xpt, 2);

            %====================%
            % Calculation starts %
            %====================%

            rhobeg = max(abs(xpt(:, 2)), [], 'all'); % Read RHOBEG from XPT.
            rhosq = rhobeg ^ 2;

            % Set BMAT.
            recip = 1.0 / rhobeg;
            reciq = 0.5 / rhobeg;
            bmat = zeros(size(bmat));
            if npt <= 2 * n + 1
                % Set BMAT(1 : NPT-N-1, :)
                bmat(1:npt - n - 1, 2:npt - n) = reciq * eye(npt - n - 1);
                bmat(1:npt - n - 1, n + 2:npt) = -reciq * eye(npt - n - 1);
                % Set BMAT(NPT-N : N, :)
                bmat(npt - n:n, 1) = -recip;
                bmat(npt - n:n, npt - n + 1:n + 1) = recip * eye(2 * n - npt + 1);
                bmat(npt - n:n, 2 * npt - n:npt + n) = -(0.5 * rhosq) * eye(2 * n - npt + 1);
            else
                bmat(:, 2:n + 1) = reciq * eye(n);
                bmat(:, n + 2:2 * n + 1) = -reciq * eye(n);
            end

            % Set ZMAT.
            recip = 1.0 / rhosq;
            reciq = sqrt(0.5) / rhosq;
            zmat = zeros(size(zmat));
            if npt <= 2 * n + 1
                zmat(1, :) = -reciq - reciq;
                zmat(2:npt - n, :) = reciq * eye(npt - n - 1);
                zmat(n + 2:npt, :) = reciq * eye(npt - n - 1);
            else
                % Set ZMAT(:, 1:N).
                zmat(1, 1:n) = -reciq - reciq;
                zmat(2:n + 1, 1:n) = reciq * eye(n);
                zmat(n + 2:2 * n + 1, 1:n) = reciq * eye(n);
                % Set ZMAT(:, N+1 : NPT-N-1).
                zmat(1, n + 1:npt - n - 1) = recip;
                zmat(2 * n + 2:npt, n + 1:npt - n - 1) = recip * eye(npt - 2 * n - 1);
                for k = 1:npt - 2 * n - 1
                    zmat(ij(:, k) + 1, k + n) = -recip;
                end
            end

            % Set IDZ.
            idz = 1;

            ipObj = inputParser();
            addParameter(ipObj, 'info', NaN);
            parse(ipObj, varargin{:});
            info = ipObj.Results.info;
            if nargout >= 4
                if any(isnan(bmat), 'all') || any(isnan(zmat), 'all')
                    info = -3;
                else
                    info = 0;
                end
            end

            %====================%
            %  Calculation ends  %
            %====================%


        end

    end
end