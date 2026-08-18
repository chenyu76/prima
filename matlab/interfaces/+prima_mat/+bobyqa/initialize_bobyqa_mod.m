classdef initialize_bobyqa_mod
    %--------------------------------------------------------------------------------------------------%
    % This module performs the initialization of BOBYQA, described in Section 2 of the BOBYQA paper.
    %
    % Coded by Zaikun ZHANG (www.zhangzk.net) based on Powell's code and the BOBYQA paper.
    %
    % Dedicated to the late Professor M. J. D. Powell FRS (1936--2015).
    %
    % Started: February 2022
    %
    % Last Modified: Tue 10 Feb 2026 01:53:25 PM CET
    %--------------------------------------------------------------------------------------------------%

    methods
        function [x0, ij, kopt, nf, fhist, fval, sl, su, xbase, xhist, xpt, info] = initxf(~, calfun, iprint, maxfun, ftarget, rhobeg, xl, xu, x0, ij, fhist, fval, sl, su, xbase, xhist, xpt)
            %--------------------------------------------------------------------------------------------------%
            % This subroutine does the initialization about the interpolation points & their function values.
            %
            % N.B.:
            % 1. Remark on IJ:
            % If NPT <= 2*N + 1, then IJ is empty. Assume that NPT >= 2*N + 2. Then SIZE(IJ) = [2, NPT-2*N-1].
            % IJ contains integers between 1 and N. For each K > 2*N + 1, XPT(:, K) is
            % XPT(:, IJ(1, K) + 1) + XPT(:, IJ(2, K) + 1). The 1 in IJ + 1 comes from the fact that XPT(:, 1)
            % corresponds to the base point XBASE. Let I = IJ(1, K) and J = IJ(2, K). Then all the
            % entries of XPT(:, K) are zero except for the I and J entries. Consequently, the Hessian of the
            % quadratic model will get a possibly nonzero (I, J) entry.
            % 2. At return,
            % INFO = INFO_DFT: initialization finishes normally
            % INFO = FTARGET_ACHIEVED: return because F <= FTARGET
            % INFO = NAN_INF_X: return because X contains NaN
            % INFO = NAN_INF_F: return because F is either NaN or +Inf
            %--------------------------------------------------------------------------------------------------%

            % Common modules
            checkexit_obj = prima_mat.common.checkexit_mod();


            evaluate_obj = prima_mat.common.evaluate_mod();
            history_obj = prima_mat.common.history_mod();


            message_obj = prima_mat.common.message_mod();

            powalg_obj = prima_mat.common.powalg_mod();
            xinbd_obj = prima_mat.common.xinbd_mod();


            % Inputs
            % N.B.: INTENT cannot be specified if a dummy procedure is not a POINTER



            % XL(N)
            % XU(N)

            % In-outputs
            % X(N)

            % Outputs
            % IJ(2, MAX(0_IK, NPT-2*N-1_IK))



            % FVAL(NPT)
            % SL(N)
            % SU(N)
            % XBASE(N)
            % XHIST(N, MAXXHIST)
            % XPT(N, NPT)

            % Local variables
            solver = "BOBYQA";


            evaluated = false(size(xpt, 2), 1);

            x = NaN(size(xpt, 1), 1);

            % Sizes.
            n = size(xpt, 1);
            npt = size(xpt, 2);


            % Preconditions


            %====================%
            % Calculation starts %
            %====================%

            % Initialize INFO to the default value. At return, an INFO different from this value will indicate
            % an abnormal return.
            info = 0;

            % SL and SU are the lower and upper bounds on feasible moves from X0.
            sl(:) = xl - x0;
            su(:) = xu - x0;
            % After the preprocessing subroutine PREPROC, SL <= 0 and the nonzero entries of SL should be less
            % than -RHOBEG, while SU >= 0 and the nonzeros of SU should be larger than RHOBEG. However, this may
            % not be true due to rounding. The following lines revise SL and SU to ensure it. X0 is also revised
            % accordingly. In precise arithmetic, the "revisions" do not change SL, SU, or X0.
            mask00 = sl < 0; %Unsupported statement inside WHERE block: StatementLineBreak 1
            sl(mask00) = min(sl(mask00), -rhobeg); %Unsupported statement inside WHERE block: StatementLineBreak 1
            mask01 = ~mask00; %Unsupported statement inside WHERE block: StatementLineBreak 1
            x0(mask01) = xl(mask01); %Unsupported statement inside WHERE block: StatementLineBreak 1
            sl(mask01) = 0.0; %Unsupported statement inside WHERE block: StatementLineBreak 1
            su(mask01) = xu(mask01) - xl(mask01); %Unsupported statement inside WHERE block: StatementLineBreak 1

            mask00 = su > 0; %Unsupported statement inside WHERE block: StatementLineBreak 1
            su(mask00) = max(su(mask00), rhobeg); %Unsupported statement inside WHERE block: StatementLineBreak 1
            mask01 = ~mask00; %Unsupported statement inside WHERE block: StatementLineBreak 1
            x0(mask01) = xu(mask01); %Unsupported statement inside WHERE block: StatementLineBreak 1
            sl(mask01) = xl(mask01) - xu(mask01); %Unsupported statement inside WHERE block: StatementLineBreak 1
            su(mask01) = 0.0; %Unsupported statement inside WHERE block: StatementLineBreak 1

            %%MATLAB code for revising X, SL, and SU:
            %%sl(sl < 0) = min(sl(sl < 0), -rhobeg);
            %%x0(sl >= 0) = xl(sl >= 0);
            %%sl(sl >= 0) = 0;
            %%su(sl >= 0) = xu(sl >= 0) - xl(sl >= 0);
            %%su(su > 0) = max(su(su > 0), rhobeg);
            %%x0(su <= 0) = xu(su <= 0);
            %%sl(su <= 0) = xl(su <= 0) - xu(su <= 0);
            %%su(su <= 0) = 0;

            % Initialize XBASE to X0.
            xbase(:) = x0;

            % EVALUATED is a boolean array with EVALUATED(I) indicating whether the function value of the I-th
            % interpolation point has been evaluated. We need it for a portable counting of the number of
            % function evaluations, especially if the loop is conducted asynchronously.
            evaluated(:) = false;

            % Initialize XHIST, FHIST, and FVAL. Otherwise, compilers may complain that they are not
            % (completely) initialized if the initialization aborts due to abnormality (see CHECKEXIT).
            % N.B.: 1. Initializing them to NaN would be more reasonable (NaN is not available in Fortran).
            % 2. Do not initialize the models if the current initialization aborts due to abnormality. Otherwise,
            % errors or exceptions may occur, as FVAL and XPT etc are uninitialized.
            xhist = repmat(-realmax, size(xhist));
            fhist(:) = realmax;
            fval(:) = realmax;

            % Set XPT(:, 2 : N+1)
            xpt = zeros(size(xpt));
            for k = 1:n
                xpt(k, k + 1) = rhobeg;
                if su(k) <= 0
                    % SU(K) == 0
                    xpt(k, k + 1) = -rhobeg;
                end
            end
            % Set XPT(:, N+2 : MIN(2*N + 1, NPT)).
            for k = 1:min(npt - n - 1, n)
                xpt(k, k + n + 1) = -rhobeg;
                if sl(k) >= 0
                    % SL(K) == 0
                    xpt(k, k + n + 1) = min(2.0 * rhobeg, su(k));
                end
                if su(k) <= 0
                    % SU(K) == 0
                    xpt(k, k + n + 1) = max(-2.0 * rhobeg, sl(k));
                end
            end

            % Set FVAL(1 : MIN(2*N + 1, NPT)) by evaluating F. Totally parallelizable except for FMSG.
            for k = 1:min(npt, 2 * n + 1)
                x(:) = xinbd_obj.xinbd(xbase, xpt(:, k), xl, xu, sl, su); % In precise arithmetic, X = XBASE + XPT(:, K).
                f = evaluate_obj.evaluatef(calfun, x);

                % Print a message about the function evaluation according to IPRINT.
                message_obj.fmsg(solver, "Initialization", iprint, k, rhobeg, f, x);
                % Save X, F into the history.
                [xhist, fhist] = history_obj.savehist(k, x, xhist, f, fhist);

                evaluated(k) = true;
                fval(k) = f;

                % Check whether to exit
                subinfo = checkexit_obj.checkexit_unc(maxfun, k, f, ftarget, x);
                if subinfo ~= 0
                    info = subinfo;
                    break
                end
            end

            % For the K between 2 and N + 1, switch XPT(:, K) and XPT(:, K+1) if XPT(K-1, K) and XPT(K-1, K+N)
            % have different signs and FVAL(K) <= FVAL(K+N). This provides a bias towards potentially lower
            % values of F when defining XPT(:, 2*N + 2 : NPT). We may drop the requirement on the signs, but
            % Powell's code has such a requirement.
            % N.B.:
            % 1. The switching is OPTIONAL. If we remove it, then the evaluations of FVAL(1 : NPT) can be
            % merged, and they are totally PARALLELIZABLE; this can be beneficial if the function evaluations
            % are expensive, which is likely the case.
            % 2. The initialization of NEWUOA revises IJ (see below) instead of XPT and FVAL. Theoretically, it
            % is equivalent; practically, after XPT is revised, the initialization of the quadratic model and
            % the Lagrange polynomials also needs revision, as, e.g., XPT(:, 2:N) is not RHOBEG*EYE(N) anymore.
            for k = 2:min(npt - n, n + 1)
                if xpt(k - 1, k) * xpt(k - 1, k + n) < 0 && fval(k + n) < fval(k)
                    fval([k, k + n]) = fval([k + n, k]);
                    xpt(:, [k, k + n]) = xpt(:, [k + n, k]);
                    % Indeed, only XPT(K-1, [K, K+N]) needs switching, as the other entries are zero.

                end
            end

            % Set IJ.
            % In general, when NPT = (N+1)*(N+2)/2, we can set IJ(:, 1 : NPT - (2*N+1)) to ANY permutation
            % of {{I, J} : 1 <= I /= J <= N}; when NPT < (N+1)*(N+2)/2, we can set it to the first NPT - (2*N+1)
            % elements of such a permutation. The following IJ is defined according to Powell's code. See also
            % Section 3 of the NEWUOA paper and (2.4) of the BOBYQA paper.
            ij(:, :) = powalg_obj.setij(n, npt);

            % Set XPT(:, 2*N + 2 : NPT). It depends on XPT(:, 1 : 2*N + 1) and hence on FVAL(1: 2*N + 1).
            % Indeed, XPT(:, K) has only two nonzeros for each K >= 2*N+2.
            % N.B.: The 1 in IJ + 1 comes from the fact that XPT(:, 1) corresponds to XBASE.
            xpt(:, 2 * n + 2:npt) = xpt(:, ij(1, :) + 1) + xpt(:, ij(2, :) + 1);

            % Set FVAL(2*N + 2 : NPT) by evaluating F. Totally parallelizable except for FMSG.
            if info == 0
                for k = 2 * n + 2:npt
                    x(:) = xinbd_obj.xinbd(xbase, xpt(:, k), xl, xu, sl, su); % In precise arithmetic, X = XBASE + XPT(:, K).
                    f = evaluate_obj.evaluatef(calfun, x);

                    % Print a message about the function evaluation according to IPRINT.
                    message_obj.fmsg(solver, "Initialization", iprint, k, rhobeg, f, x);
                    % Save X, F into the history.
                    [xhist, fhist] = history_obj.savehist(k, x, xhist, f, fhist);

                    evaluated(k) = true;
                    fval(k) = f;

                    % Check whether to exit
                    subinfo = checkexit_obj.checkexit_unc(maxfun, k, f, ftarget, x);
                    if subinfo ~= 0
                        info = subinfo;
                        break
                    end
                end
            end

            % Set NF, KOPT
            nf = nnz(evaluated);
            kopt = fortran.minloc(fval, 'mask', evaluated, 'dim', 1);
            %%MATLAB: fopt = min(fval(evaluated)); kopt = find(evaluated & ~(fval > fopt), 1, 'first');

            %====================%
            %  Calculation ends  %
            %====================%

            % Postconditions


        end
        function [gopt, hq, pq, info] = initq(~, ij, fval, xpt, gopt, hq, pq, varargin)
            %--------------------------------------------------------------------------------------------------%
            % This subroutine initializes the quadratic model represented by [GOPT, HQ, PQ] so that its gradient
            % at XBASE + XPT(:, KOPT) is GOPT; its Hessian is HQ + sum_{K=1}^NPT PQ(K)*XPT(:, K)*XPT(:, K)'.
            %--------------------------------------------------------------------------------------------------%

            % Common modules



            % Inputs
            % IJ(2, MAX(0_IK, NPT-2*N-1_IK))
            % FVAL(NPT)
            % XPT(N, NPT)

            % Outputs

            % GOPT(N)
            % HQ(N, N)
            % PQ(NPT)

            % Local variables



            xa = NaN(min(size(xpt, 1), size(xpt, 2) - size(xpt, 1) - 1), 1);
            xb = NaN(numel(xa), 1);


            % Sizes
            n = size(xpt, 1);
            npt = size(xpt, 2);

            % Preconditions


            %====================%
            % Calculation starts %
            %====================%

            fbase = fval(1); % FBASE is the function value at XBASE.

            % Set GOPT by the forward difference.
            gopt(:) = (fval(2:n + 1) - fbase) ./ diag(xpt(:, 2:n + 1));

            % The interpolation conditions decide GOPT(1:NDIAG) and the first NDIAG diagonal 2nd derivatives of
            % the initial quadratic model by a quadratic interpolation on three points.
            ndiag = min(n, npt - n - 1);
            xa(:) = diag(xpt(:, 2:ndiag + 1));
            xb(:) = diag(xpt(:, n + 2:n + ndiag + 1));

            % Revise GOPT(1:NDIAG) to the value provided by the three-point interpolation.
            gopt(1:ndiag) = (gopt(1:ndiag) .* xb - ((fval(n + 2:n + ndiag + 1) - fbase) ./ xb) .* xa) ./ (xb - xa);

            % Set the diagonal of HQ by the three-point interpolation. If we do this before the revision of
            % GOPT(1:NDIAG), we can avoid the calculation of FVAL(K + 1) - FBASE) / RHOBEG. But we prefer to
            % decouple the initialization of GOPT and HQ. We are not concerned by this amount of flops.
            hq = zeros(size(hq));
            for k = 1:ndiag
                hq(k, k) = 2.0 * ((fval(k + 1) - fbase) / xa(k) - (fval(n + k + 1) - fbase) / xb(k)) / (xa(k) - xb(k));
            end
            %%MATLAB:
            %%hdiag = 2*((fval(2 : ndiag+1) - fbase) / xa - (fval(n+2 : n+ndiag+1) - fbase) / xb) / (xa-xb)
            %%hq(1:ndiag, 1:ndiag) = diag(hdiag)

            % When NPT > 2*N + 1, set the off-diagonal entries of HQ.
            for k = 1:npt - 2 * n - 1
                i = ij(1, k);
                j = ij(2, k);
                xi = xpt(i, k + 2 * n + 1);
                xj = xpt(j, k + 2 * n + 1);
                % N.B.: The 1 in I+1 and J+1 comes from the fact that XPT(:, 1) corresponds to XBASE.
                hq(i, j) = (fbase - fval(i + 1) - fval(j + 1) + fval(k + 2 * n + 1)) / (xi * xj);
                hq(j, i) = hq(i, j);
            end

            [~, kopt] = min(fval);
            if kopt ~= 1
                gopt(:) = gopt + hq * xpt(:, kopt);
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

            % Postconditions


        end
        function [bmat, zmat, info] = inith(~, ij, xpt, bmat, zmat, varargin)
            %--------------------------------------------------------------------------------------------------%
            % This subroutine initializes [BMAT, ZMAT] which represents the matrix H in (2.7) of the BOBYQA
            % paper (see also (3.12) of the NEWUOA paper).
            %--------------------------------------------------------------------------------------------------%

            % Common modules



            %use, non_intrinsic :: powalg_mod, only : errh


            % Inputs
            % IJ(2, MAX(0_IK, NPT-2*N-1_IK))
            % XPT(N, NPT)

            % Outputs

            % BMAT(N, NPT + N)
            % ZMAT(NPT, NPT - N - 1)

            % Local variables



            xa = NaN(min(size(xpt, 1), size(xpt, 2) - size(xpt, 1) - 1), 1);
            xb = NaN(numel(xa), 1);

            % Sizes
            n = size(xpt, 1);
            npt = size(xpt, 2);

            % Preconditions


            %====================%
            % Calculation starts %
            %====================%

            % Some values to be used for setting BMAT and ZMAT.
            rhobeg = max(abs(xpt(:, 2)), [], 'all'); % Read RHOBEG from XPT. Note that XPT(:, 1) = 0.
            rhosq = rhobeg ^ 2;

            % The interpolation set decides the first NDIAG diagonal 2nd derivatives of the Lagrange polynomials.
            ndiag = min(n, npt - n - 1);
            xa(:) = diag(xpt(:, 2:ndiag + 1));
            xb(:) = diag(xpt(:, n + 2:n + ndiag + 1));

            bmat = zeros(size(bmat));
            % Set BMAT(1 : NDIAG, :)
            bmat(1:ndiag, 1) = -(xa + xb) ./ (xa .* xb);
            for k = 1:ndiag
                bmat(k, k + n + 1) = -0.5 / xpt(k, k + 1);
                bmat(k, k + 1) = -bmat(k, 1) - bmat(k, k + n + 1);
            end
            % Set BMAT(NDIAG+1 : N, :)
            for k = ndiag + 1:n
                bmat(k, 1) = -1.0 / xpt(k, k + 1);
                bmat(k, k + 1) = -bmat(k, 1);
                bmat(k, npt + k) = -0.5 * rhosq;
            end

            zmat = zeros(size(zmat));
            % Set ZMAT(:, 1 : NDIAG)
            zmat(1, 1:ndiag) = sqrt(2.0) ./ (xa .* xb);
            for k = 1:ndiag
                zmat(k + 1, k) = -zmat(1, k) - sqrt(0.5) / rhosq;
                zmat(k + n + 1, k) = sqrt(0.5) / rhosq;
            end
            % Set ZMAT(:, NDIAG+1 : NPT-N-1)
            for k = ndiag + 1:npt - n - 1
                zmat(1, k) = 1.0 / rhosq;
                zmat(k + n + 1, k) = 1.0 / rhosq;
                zmat(ij(:, k - n) + 1, k) = -1.0 / rhosq;
            end

            ipObj = inputParser();
            addParameter(ipObj, 'info', NaN);
            parse(ipObj, varargin{:});
            info = ipObj.Results.info;
            if nargout >= 3
                if any(isnan(bmat), 'all') || any(isnan(zmat), 'all')
                    info = -3;
                else
                    info = 0;
                end
            end

            %====================%
            %  Calculation ends  %
            %====================%

            % Postconditions


        end

    end
end