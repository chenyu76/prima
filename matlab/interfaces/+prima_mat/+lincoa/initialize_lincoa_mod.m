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


            checkexit_obj = prima_mat.common.checkexit_mod();

            evaluate_obj = prima_mat.common.evaluate_mod();
            history_obj = prima_mat.common.history_mod();

            message_obj = prima_mat.common.message_mod();

            powalg_obj = prima_mat.common.powalg_mod();

            % AMAT(Meq, N)
            % AMAT(Mineq, N)
            % AMAT(N, M)
            % Beq(M)
            % Bineq(M)


            % XL(N)
            % XU(N)
            % X0(N)


            % B(M)


            % IJ(2, MAX(0_IK, NPT-2*N-1))


            % EVALUATED(NPT)
            % CHIST(MAXCHIST)
            % CVAL(NPT)
            % FHIST(MAXFHIST)
            % FVAL(NPT)
            % XBASE(N)
            % XHIST(N, MAXXHIST)
            % XPT(N, NPT)


            solver = "LINCOA";

            feasible = false(size(xpt, 2), 1);
            constr = NaN(nnz(xl > -(0.25 * realmax)) + nnz(xu < 0.25 * realmax) + 2 * numel(beq) + numel(bineq), 1);
            constr_leq = NaN(numel(beq), 1);

            x = NaN(numel(x0), 1);

            n = size(xpt, 1);
            npt = size(xpt, 2);

            %====================%
            % Calculation starts %
            %====================%

            % Initialize INFO to the default value. At return, an INFO different from this value will indicate
            % an abnormal return.
            info = 0;

            % Initialize XBASE to X0.
            xbase(:) = x0;

            % EVALUATED is a boolean array with EVALUATED(I) indicating whether the function value of the I-th
            % interpolation point has been evaluated. We need it for a portable counting of the number of
            % function evaluations, especially if the loop is conducted asynchronously. However, the loop here
            % is not fully parallelizable if NPT>2N+1, as the definition XPT(;, 2N+2:end) involves FVAL(1:2N+1).
            evaluated(:) = false;

            % Initialize XHIST, FHIST, CHIST, FVAL, and CVAL. Otherwise, compilers may complain that they are
            % not (completely) initialized if the initialization aborts due to abnormality (see CHECKEXIT).
            % N.B.: 1. Initializing them to NaN would be more reasonable (NaN is not available in Fortran).
            % 2. Do not initialize the models if the current initialization aborts due to abnormality. Otherwise,
            % errors or exceptions may occur, as FVAL and XPT etc are uninitialized.
            xhist = repmat(-realmax, size(xhist));
            fhist(:) = realmax;
            chist(:) = realmax;
            fval(:) = realmax;
            cval(:) = realmax;

            % Set the nonzero coordinates of XPT(K,.), K=1,2,...,min[2*N+1,NPT], but they may be altered
            % later to make a constraint violation sufficiently large.
            xpt(:, 1) = 0.0;
            xpt(:, 2:n + 1) = rhobeg * eye(n);
            xpt(:, n + 2:npt) = -rhobeg * eye(n, npt - n - 1); % XPT(:, 2*N+2 : NPT) = ZERO if it is nonempty.

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
            b(:) = b - amat.' * xbase;

            % Define FEASIBLE, which will be used when defining KOPT.
            for k = 1:npt
                % Internally, we use AMAT and B to evaluate the constraints.
                cval(k) = max([0.0; amat.' * xpt(:, k) - b], [], 'all');
                if isnan(cval(k))
                    cval(k) = realmax;
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
            % Removable in F2003.
            % Removable in F2003.
            ixl = find(xl > -(0.25 * realmax));
            ixu = find(xu < 0.25 * realmax);
            for k = 1:npt
                x(:) = xbase + xpt(:, k);
                f = evaluate_obj.evaluatef(calfun, x);
                % Evaluate the constraints.
                constr_leq(:) = Aeq * x - beq;
                constr(:) = [xl(ixl) - x(ixl); x(ixu) - xu(ixu); -constr_leq; constr_leq; Aineq * x - bineq];
                cstrv = max([0.0; constr], [], 'all');

                % Print a message about the function evaluation according to IPRINT.
                message_obj.fmsg(solver, "Initialization", iprint, k, rhobeg, f, x, 'cstrv', cstrv, 'constr', constr);
                % Save X, F, CSTRV into the history.
                [xhist, fhist, chist] = history_obj.savehist(k, x, xhist, f, fhist, 'cstrv', cstrv, 'chist', chist);

                evaluated(k) = true;
                cval(k) = cstrv; % CVAL will be used to initialize CFILT.
                fval(k) = f;

                % Check whether to exit.
                subinfo = checkexit_obj.checkexit_con(maxfun, k, cstrv, ctol, f, ftarget, x);
                if subinfo ~= 0
                    info = subinfo;
                    break
                end
            end

            % Deallocate IXL and IXU as they have finished their job.


            nf = nnz(evaluated);
            % Since the starting point is supposed to be feasible, there should be at least one feasible point.
            % We set feasible to TRUE for the evaluated point with the smallest constraint violation. This is
            % necessary, or KOPT defined below may become 0 if EVALUATED .AND. FEASIBLE is all FALSE.
            feasible(fortran.minloc(cval, 'mask', evaluated, 'dim', 1)) = true;
            kopt = fortran.minloc(fval, 'mask', (evaluated & feasible), 'dim', 1);
            %%MATLAB:
            %%fopt = min(fval(evaluated & feasible));
            %%kopt = find(evaluated & feasible & ~(fval > fopt), 1, 'first');

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