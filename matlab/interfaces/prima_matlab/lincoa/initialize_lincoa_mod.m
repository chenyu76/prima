% FIXME: The definitions of CVAL, FEASIBLE, and KOPT are questionable.
classdef initialize_lincoa_mod
    %--------------------------------------------------------------------------------------------------%
    % This module performs the initialization of LINCOA.
    %
    % Coded by Zaikun ZHANG (www.zhangzk.net) based on Powell's code.
    %
    % Dedicated to the late Professor M. J. D. Powell FRS (1936--2015).
    %
    % Started: February 2022
    %
    % Last Modified: Tue 10 Feb 2026 02:41:49 PM CET
    %--------------------------------------------------------------------------------------------------%

    methods
        function [b, ij, kopt, nf, chist, cval, fhist, fval, xbase, xhist, xpt, evaluated, info] = initxf(~, calfun, iprint, maxfun, Aeq, Aineq, amat, beq, bineq, ctol, ftarget, rhobeg, xl, xu, x0, b, ij, chist, cval, fhist, fval, xbase, xhist, xpt, evaluated)
            %--------------------------------------------------------------------------------------------------%
            % This subroutine does the initialization about the interpolation points & their function values.
            %
            % N.B.:
            % 1. Remark on IJ:
            % If NPT <= 2*N + 1, then IJ is empty. Assume that NPT >= 2*N + 2. Then SIZE(IJ) = [2, NPT-2*N-1].
            % IJ contains integers between 1 and N. For each K > 2*N + 1, XPT(:, K) is
            % XPT(:, IJ(1, K) + 1) + XPT(:, IJ(2, K) + 1). The 1 in IJ + 1 comes from the fact that XPT(:, 1)
            % corresponds to the base point XBASE. Let I = IJ(1, K) and J = IJ(2, K). Then all the entries of
            % XPT(:, K) are zero except for the I and J entries. Consequently, the Hessian of the quadratic
            % model will get a possibly nonzero (I, J) entry.
            % 2. At return,
            % INFO = INFO_DFT: initialization finishes normally
            % INFO = FTARGET_ACHIEVED: return because F <= FTARGET
            % INFO = NAN_INF_X: return because X contains NaN
            % INFO = NAN_INF_F: return because F is either NaN or +Inf
            %--------------------------------------------------------------------------------------------------%

            % Common modules
            checkexit_obj = checkexit_mod();
            consts_obj = consts_mod();
            debug_obj = debug_mod();
            evaluate_obj = evaluate_mod();
            history_obj = history_mod();
            infnan_obj = infnan_mod();
            infos_obj = infos_mod();
            linalg_obj = linalg_mod();
            memory_obj = memory_mod();
            message_obj = message_mod();
            pintrf_obj = pintrf_mod();
            powalg_obj = powalg_mod();


            % Inputs
            % N.B.: INTENT cannot be specified if a dummy procedure is not a POINTER


            % AMAT(Meq, N)
            % AMAT(Mineq, N)
            % AMAT(N, M)
            % Beq(M)
            % Bineq(M)



            % XL(N)
            % XU(N)
            % X0(N)

            % In-outputs
            % B(M)

            % Outputs

            % IJ(2, MAX(0_IK, NPT-2*N-1))


            % EVALUATED(NPT)
            % CHIST(MAXCHIST)
            % CVAL(NPT)
            % FHIST(MAXFHIST)
            % FVAL(NPT)
            % XBASE(N)
            % XHIST(N, MAXXHIST)
            % XPT(N, NPT)

            % Local variables
            solver = "LINCOA";
            srname = "INITXF";
            k = NaN;
            m = NaN;
            maxchist = NaN;
            maxfhist = NaN;
            maxhist = NaN;
            maxxhist = NaN;
            n = NaN;
            npt = NaN;
            subinfo = NaN;
            ixl = NaN(1);
            ixu = NaN(1);
            feasible = false(size(xpt, 2), 1);
            constr = NaN(nnz(xl > -consts_obj.BOUNDMAX) + nnz(xu < consts_obj.BOUNDMAX) + 2 * numel(beq) + numel(bineq), 1);
            constr_leq = NaN(numel(beq), 1);
            cstrv = NaN;
            f = NaN;
            x = NaN(numel(x0), 1);

            % Sizes.
            m = fix(numel(b));
            n = size(xpt, 1);
            npt = size(xpt, 2);
            maxxhist = size(xhist, 2);
            maxfhist = fix(numel(fhist));
            maxchist = fix(numel(chist));
            maxhist = fix(max(maxxhist, max(maxfhist, maxchist)));

            % Preconditions
            if consts_obj.DEBUGGING
                debug_obj.assert(abs(iprint) <= 3, "IPRINT is 0, 1, -1, 2, -2, 3, or -3", srname);
                debug_obj.assert(m >= 0, "M >= 0", srname);
                debug_obj.assert(n >= 1, "N >= 1", srname);
                debug_obj.assert(npt >= n + 2, "NPT >= N+2", srname);
                debug_obj.assert(size(Aeq, 1) == numel(beq) && size(Aeq, 2) == n, "SIZE(Aeq) == [SIZE(Beq), M]", srname);
                debug_obj.assert(size(Aineq, 1) == numel(bineq) && size(Aineq, 2) == n, "SIZE(Aineq) == [SIZE(Bineq), M]", srname);
                debug_obj.assert(size(amat, 1) == n && size(amat, 2) == m, "SIZE(AMAT) == [N, M]", srname);
                debug_obj.assert(rhobeg > 0, "RHOBEG > 0", srname);
                debug_obj.assert(numel(xbase) == n, "SIZE(XBASE) == N", srname);
                debug_obj.assert(numel(xl) == n && numel(xu) == n, "SIZE(XL) == N == SIZE(XU)", srname);
                debug_obj.assert(numel(x0) == n && all(infnan_obj.is_finite(x0), 'all'), "SIZE(X0) == N, X0 is finite", srname);
                debug_obj.assert(numel(fval) == npt, "SIZE(FVAL) == NPT", srname);
                debug_obj.assert(numel(cval) == npt, "SIZE(CVAL) == NPT", srname);
                debug_obj.assert(numel(evaluated) == npt, "SIZE(EVALUATED) == NPT", srname);
                debug_obj.assert(size(xhist, 1) == n && maxxhist * (maxxhist - maxhist) == 0, "SIZE(XHIST, 1) == N, SIZE(XHIST, 2) == 0 or MAXHIST", srname);
                debug_obj.assert(maxfhist * (maxfhist - maxhist) == 0, "SIZE(FHIST) == 0 or MAXHIST", srname);
                debug_obj.assert(maxchist * (maxchist - maxhist) == 0, "SIZE(CHIST) == 0 or MAXHIST", srname);
            end

            %====================%
            % Calculation starts %
            %====================%

            % Initialize INFO to the default value. At return, an INFO different from this value will indicate
            % an abnormal return.
            info = infos_obj.INFO_DFT;

            % Initialize XBASE to X0.
            xbase(:) = x0;

            % EVALUATED is a boolean array with EVALUATED(I) indicating whether the function value of the I-th
            % interpolation point has been evaluated. We need it for a portable counting of the number of
            % function evaluations, especially if the loop is conducted asynchronously. However, the loop here
            % is not fully parallelizable if NPT>2N+1, as the definition XPT(;, 2N+2:end) involves FVAL(1:2N+1).
            evaluated = repmat(false, size(evaluated));

            % Initialize XHIST, FHIST, CHIST, FVAL, and CVAL. Otherwise, compilers may complain that they are
            % not (completely) initialized if the initialization aborts due to abnormality (see CHECKEXIT).
            % N.B.: 1. Initializing them to NaN would be more reasonable (NaN is not available in Fortran).
            % 2. Do not initialize the models if the current initialization aborts due to abnormality. Otherwise,
            % errors or exceptions may occur, as FVAL and XPT etc are uninitialized.
            xhist = repmat(-consts_obj.REALMAX, size(xhist));
            fhist = repmat(consts_obj.REALMAX, size(fhist));
            chist = repmat(consts_obj.REALMAX, size(chist));
            fval = repmat(consts_obj.REALMAX, size(fval));
            cval = repmat(consts_obj.REALMAX, size(cval));

            % Set the nonzero coordinates of XPT(K,.), K=1,2,...,min[2*N+1,NPT], but they may be altered
            % later to make a constraint violation sufficiently large.
            xpt(:, 1) = consts_obj.ZERO;
            xpt(:, 2:n + 1) = rhobeg * linalg_obj.eye1(n);
            xpt(:, n + 2:npt) = -rhobeg * linalg_obj.eye2(n, npt - n - 1); % XPT(:, 2*N+2 : NPT) = ZERO if it is nonempty.

            % Set IJ.
            % In general, when NPT = (N+1)*(N+2)/2, we can set IJ(:, 1 : NPT - (2*N+1)) to ANY permutation
            % of {{I, J} : 1 <= I /= J <= N}; when NPT < (N+1)*(N+2)/2, we can set it to the first NPT - (2*N+1)
            % elements of such a permutation. The following IJ is defined according to Powell's code. See also
            % Section 3 of the NEWUOA paper and (2.4) of the BOBYQA paper.
            ij(:, :) = powalg_obj.setij(n, npt);

            % Set XPT(:, 2*N + 2 : NPT).
            % Indeed, XPT(:, K) has only two nonzeros for each K >= 2*N + 2,
            % N.B.: The 1 in IJ + 1 comes from the fact that XPT(:, 1) corresponds to XBASE.
            xpt(:, 2 * n + 2:npt) = xpt(:, ij(1, :) + 1) + xpt(:, ij(2, :) + 1);

            % Update the constraint right-hand sides to allow for the shift XBASE.
            b(:) = b - linalg_obj.matprod12(xbase, amat);

            % Define FEASIBLE, which will be used when defining KOPT.
            for k = 1:npt
                % Internally, we use AMAT and B to evaluate the constraints.
                cval(k) = linalg_obj.maximum1([consts_obj.ZERO; reshape(linalg_obj.matprod12(xpt(:, k), amat) - b, [], 1)]);
                if infnan_obj.is_nan_sp(cval(k))
                    cval(k) = consts_obj.REALMAX;
                end
                % Powell's implementation contains the following procedure that shifts every infeasible point if
                % necessary so that its constraint violation is at least 0.2*RHOBEG. According to a test on
                % 20230209, it does not evidently improve the performance of LINCOA. Indeed, it worsens a bit the
                % performance in the early stage. Thus we decided to remove it.
                %----------------------------------------------------------------------------------------------%
                %mincv = 0.2_RP * rhobeg
                %constr(1:m) = matprod(xpt(:, k), amat) - b
                %if (cval(k) < mincv .and. cval(k) > 0) then
                %    j = int(maxloc(constr(1:m), dim=1), kind(j))
                %    xpt(:, k) = xpt(:, k) + (mincv - constr(j)) * amat(:, j)
                %end if
                %----------------------------------------------------------------------------------------------%
            end
            feasible(:) = (cval <= 0);

            % Set FVAL by evaluating F. Totally parallelizable except for FMSG.
            % IXL and IXU are the indices of the nontrivial lower and upper bounds, respectively.
            ixl = memory_obj.alloc_ivector(ixl, fix(nnz(xl > -consts_obj.BOUNDMAX))); % Removable in F2003.
            ixu = memory_obj.alloc_ivector(ixu, fix(nnz(xu < consts_obj.BOUNDMAX))); % Removable in F2003.
            ixl = linalg_obj.trueloc(xl > -consts_obj.BOUNDMAX);
            ixu = linalg_obj.trueloc(xu < consts_obj.BOUNDMAX);
            for k = 1:npt
                x(:) = xbase + xpt(:, k);
                f = evaluate_obj.evaluatef(calfun, x);
                % Evaluate the constraints.
                constr_leq(:) = linalg_obj.matprod21(Aeq, x) - beq;
                constr(:) = [reshape(xl(ixl) - x(ixl), [], 1); reshape(x(ixu) - xu(ixu), [], 1); reshape(-constr_leq, [], 1); reshape(constr_leq, [], 1); reshape(linalg_obj.matprod21(Aineq, x) - bineq, [], 1)];
                cstrv = linalg_obj.maximum1([consts_obj.ZERO; reshape(constr, [], 1)]);

                % Print a message about the function evaluation according to IPRINT.
                message_obj.fmsg(solver, "Initialization", iprint, k, rhobeg, f, x, 'cstrv', cstrv, 'constr', constr);
                % Save X, F, CSTRV into the history.
                [xhist, fhist, chist] = history_obj.savehist(k, x, xhist, f, fhist, 'cstrv', cstrv, 'chist', chist);

                evaluated(k) = true;
                cval(k) = cstrv; % CVAL will be used to initialize CFILT.
                fval(k) = f;

                % Check whether to exit.
                subinfo = checkexit_obj.checkexit_con(maxfun, k, cstrv, ctol, f, ftarget, x);
                if subinfo ~= infos_obj.INFO_DFT
                    info = subinfo;
                    break
                end
            end

            % Deallocate IXL and IXU as they have finished their job.
            ixl = []; ixu = [];

            nf = fix(nnz(evaluated));
            % Since the starting point is supposed to be feasible, there should be at least one feasible point.
            % We set feasible to TRUE for the evaluated point with the smallest constraint violation. This is
            % necessary, or KOPT defined below may become 0 if EVALUATED .AND. FEASIBLE is all FALSE.
            feasible(fortran.minloc(cval, 'mask', evaluated, 'dim', 1)) = true;
            kopt = fix(fortran.minloc(fval, 'mask', (evaluated & feasible), 'dim', 1));
            %%MATLAB:
            %%fopt = min(fval(evaluated & feasible));
            %%kopt = find(evaluated & feasible & ~(fval > fopt), 1, 'first');

            %====================%
            %  Calculation ends  %
            %====================%

            % Postconditions
            if consts_obj.DEBUGGING
                debug_obj.assert(size(ij, 1) == 2 && size(ij, 2) == max(0, npt - 2 * n - 1), "SIZE(IJ) == [2, NPT - 2*N - 1]", srname);
                debug_obj.assert(all(ij >= 1 & ij <= 2 * n, 'all'), "1 <= IJ <= 2*N", srname);
                debug_obj.assert(all(ij(1, :) ~= ij(2, :), 'all'), "IJ(1, :) /= IJ(:, 2)", srname);
                debug_obj.assert(nf <= npt, "NF <= NPT", srname);
                debug_obj.assert(kopt >= 1 && kopt <= nf, "1 <= KOPT <= NF", srname);
                debug_obj.assert(numel(xbase) == n && all(infnan_obj.is_finite(xbase), 'all'), "SIZE(XBASE) == N, XBASE is finite", srname);
                debug_obj.assert(size(xpt, 1) == n && size(xpt, 2) == npt, "SIZE(XPT) == [N, NPT]", srname);
                debug_obj.assert(all(infnan_obj.is_finite(xpt), 'all'), "XPT is finite", srname);
                debug_obj.assert(numel(cval) == npt && ~any(evaluated & (infnan_obj.is_nan(cval) | infnan_obj.is_posinf(cval)), 'all'), "SIZE(CVAL) == NPT and CVAL is not NaN or +Inf", srname);
                debug_obj.assert(numel(fval) == npt && ~any(evaluated & (infnan_obj.is_nan(fval) | infnan_obj.is_posinf(fval)), 'all'), "SIZE(FVAL) == NPT and FVAL is not NaN or +Inf", srname);
                debug_obj.assert(~any(evaluated & feasible & fval < fval(kopt), 'all'), "FVAL(KOPT) = MINVAL(FVAL)", srname);
                debug_obj.assert(numel(fhist) == maxfhist, "SIZE(FHIST) == MAXFHIST", srname);
                debug_obj.assert(numel(chist) == maxchist, "SIZE(CHIST) == MAXCHIST", srname);
                debug_obj.assert(size(xhist, 1) == n && size(xhist, 2) == maxxhist, "SIZE(XHIST) == [N, MAXXHIST]", srname);
                % LINCOA always starts with a feasible point.
                if m > 0
                    debug_obj.assert(all(linalg_obj.matprod12(xpt(:, 1), amat) - b <= max(fortran.power(consts_obj.TEN, max(-12, -consts_obj.MAXPOW10)), 100.0 * consts_obj.EPS) * (consts_obj.ONE + sum(abs(xpt(:, 1)), 'all') + sum(abs(b), 'all')), 'all'), "The starting point is feasible", srname);
                end
            end

        end
        function [idz, bmat, zmat, info] = inith(~, ij, xpt, bmat, zmat, varargin)
            %--------------------------------------------------------------------------------------------------%
            % This subroutine initializes [IDZ, BMAT, ZMAT] which represents the matrix H in (3.12) of the
            % NEWUOA paper (see also (2.7) of the BOBYQA paper).
            %--------------------------------------------------------------------------------------------------%
            % Common modules
            consts_obj = consts_mod();
            debug_obj = debug_mod();
            infnan_obj = infnan_mod();
            infos_obj = infos_mod();
            linalg_obj = linalg_mod();
            %use, non_intrinsic :: powalg_mod, only : errh


            % Inputs
            % IJ(2, MAX(0_IK, NPT - 2_IK * N - 1_IK))
            % XPT(N, NPT)
            % N.B.: XPT is essentially only used for debugging, to test the error in the initial H. The initial
            % ZMAT and BMAT are completely defined by RHOBEG and IJ.

            % Outputs


            % BMAT(N, NPT + N)
            % ZMAT(NPT, NPT - N - 1)

            % Local variables
            srname = "INITH";
            k = NaN;
            n = NaN;
            npt = NaN;
            recip = NaN;
            reciq = NaN;
            rhobeg = NaN;
            rhosq = NaN;

            % Sizes
            n = size(xpt, 1);
            npt = size(xpt, 2);

            % Preconditions
            if consts_obj.DEBUGGING
                debug_obj.assert(n >= 1 && npt >= n + 2, "N >= 1, NPT >= N + 2", srname);
                debug_obj.assert(size(ij, 1) == 2 && size(ij, 2) == max(0, npt - 2 * n - 1), "SIZE(IJ) == [2, NPT - 2*N - 1]", srname);
                debug_obj.assert(all(ij >= 1 & ij <= 2 * n, 'all'), "1 <= IJ <= 2*N", srname);
                debug_obj.assert(all(ij(1, :) ~= ij(2, :), 'all'), "IJ(1, :) /= IJ(2, :)", srname);
                debug_obj.assert(all(infnan_obj.is_finite(xpt), 'all'), "XPT is finite", srname);
                debug_obj.assert(size(bmat, 1) == n && size(bmat, 2) == npt + n, "SIZE(BMAT)==[N, NPT+N]", srname);
                debug_obj.assert(size(zmat, 1) == npt && size(zmat, 2) == npt - n - 1, "SIZE(ZMAT) == [NPT, NPT - N - 1]", srname);
            end

            %====================%
            % Calculation starts %
            %====================%

            rhobeg = max(abs(xpt(:, 2)), [], 'all'); % Read RHOBEG from XPT.
            rhosq = fortran.power(rhobeg, 2);

            % Set BMAT.
            recip = consts_obj.ONE / rhobeg;
            reciq = consts_obj.HALF / rhobeg;
            bmat = repmat(consts_obj.ZERO, size(bmat));
            if npt <= 2 * n + 1
                % Set BMAT(1 : NPT-N-1, :)
                bmat(1:npt - n - 1, 2:npt - n) = reciq * linalg_obj.eye1(npt - n - 1);
                bmat(1:npt - n - 1, n + 2:npt) = -reciq * linalg_obj.eye1(npt - n - 1);
                % Set BMAT(NPT-N : N, :)
                bmat(npt - n:n, 1) = -recip;
                bmat(npt - n:n, npt - n + 1:n + 1) = recip * linalg_obj.eye1(2 * n - npt + 1);
                bmat(npt - n:n, 2 * npt - n:npt + n) = -(consts_obj.HALF * rhosq) * linalg_obj.eye1(2 * n - npt + 1);
            else
                bmat(:, 2:n + 1) = reciq * linalg_obj.eye1(n);
                bmat(:, n + 2:2 * n + 1) = -reciq * linalg_obj.eye1(n);
            end

            % Set ZMAT.
            recip = consts_obj.ONE / rhosq;
            reciq = fortran.sqrt(consts_obj.HALF) / rhosq;
            zmat = repmat(consts_obj.ZERO, size(zmat));
            if npt <= 2 * n + 1
                zmat(1, :) = -reciq - reciq;
                zmat(2:npt - n, :) = reciq * linalg_obj.eye1(npt - n - 1);
                zmat(n + 2:npt, :) = reciq * linalg_obj.eye1(npt - n - 1);
            else
                % Set ZMAT(:, 1:N).
                zmat(1, 1:n) = -reciq - reciq;
                zmat(2:n + 1, 1:n) = reciq * linalg_obj.eye1(n);
                zmat(n + 2:2 * n + 1, 1:n) = reciq * linalg_obj.eye1(n);
                % Set ZMAT(:, N+1 : NPT-N-1).
                zmat(1, n + 1:npt - n - 1) = recip;
                zmat(2 * n + 2:npt, n + 1:npt - n - 1) = recip * linalg_obj.eye1(npt - 2 * n - 1);
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
                if any(infnan_obj.is_nan(bmat), 'all') || any(infnan_obj.is_nan(zmat), 'all')
                    info = infos_obj.NAN_INF_MODEL;
                else
                    info = infos_obj.INFO_DFT;
                end
            end

            %====================%
            %  Calculation ends  %
            %====================%

            % Postconditions
            if consts_obj.DEBUGGING
                debug_obj.assert(idz >= 1 && idz <= size(zmat, 2) + 1, "1 <= IDZ <= SIZE(ZMAT, 2) + 1", srname);
                debug_obj.assert(size(bmat, 1) == n && size(bmat, 2) == npt + n, "SIZE(BMAT)==[N, NPT+N]", srname);
                debug_obj.assert(linalg_obj.issymmetric(bmat(:, npt + 1:npt + n)), "BMAT(:, NPT+1:NPT+N) is symmetric", srname);
                debug_obj.assert(size(zmat, 1) == npt && size(zmat, 2) == npt - n - 1, "SIZE(ZMAT) == [NPT, NPT - N - 1]", srname);
                %call assert(errh(idz, bmat, zmat, xpt) <= max(1.0E-3_RP, 1.0E2_RP * real(npt, RP) * EPS), &
                %    & '[IDZ, BMA, ZMAT] represents H = W^{-1}', srname)

            end

        end

    end
end