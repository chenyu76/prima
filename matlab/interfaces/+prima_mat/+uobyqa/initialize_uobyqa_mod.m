classdef initialize_uobyqa_mod
    %--------------------------------------------------------------------------------------------------%
    % This module performs the initialization of UOBYQA, described in Section 4 of the UOBYQA paper.
    %
    % Coded by Zaikun ZHANG (www.zhangzk.net) based on Powell's code and the UOBYQA paper.
    %
    % Started: July 2020
    %
    % Dedicated to the late Professor M. J. D. Powell FRS (1936--2015).
    %
    % Last Modified: Tue 10 Feb 2026 02:01:01 PM CET
    %--------------------------------------------------------------------------------------------------%

    methods
        function [kopt, nf, fhist, fval, xbase, xhist, xpt, info] = initxf(~, calfun, iprint, maxfun, ftarget, rhobeg, x0, fhist, fval, xbase, xhist, xpt)
            %--------------------------------------------------------------------------------------------------%
            % This subroutine does the initialization about the interpolation points & their function values.
            % See Section 4 of the UOBYQA paper.
            %--------------------------------------------------------------------------------------------------%

            % Common modules
            checkexit_obj = prima_mat.common.checkexit_mod();
            consts_obj = prima_mat.common.consts_mod();
            debug_obj = prima_mat.common.debug_mod();
            evaluate_obj = prima_mat.common.evaluate_mod();
            history_obj = prima_mat.common.history_mod();
            infnan_obj = prima_mat.common.infnan_mod();
            infos_obj = prima_mat.common.infos_mod();
            linalg_obj = prima_mat.common.linalg_mod();
            message_obj = prima_mat.common.message_mod();
            pintrf_obj = prima_mat.common.pintrf_mod();


            % Inputs
            % N.B.: INTENT cannot be specified if a dummy procedure is not a POINTER



            % X0(N)

            % Outputs



            % FHIST(MAXFHIST)

            % XBASE(N)
            % XHIST(N, MAXXHIST)
            % XPT(N, NPT)

            % Local variables
            solver = "UOBYQA";
            srname = "INITXF";
            ip = NaN;
            iq = NaN;
            k = NaN;
            kk = NaN(numel(x0), 1);
            maxfhist = NaN;
            maxhist = NaN;
            maxxhist = NaN;
            n = NaN;
            npt = NaN;
            subinfo = NaN;
            evaluated = false(size(xpt, 2), 1);
            f = NaN;
            x = NaN(numel(x0), 1);
            xw = NaN(numel(x0), 1);

            % Sizes
            n = size(xpt, 1);
            npt = size(xpt, 2);
            maxxhist = size(xhist, 2);
            maxfhist = fix(numel(fhist));
            maxhist = max(maxxhist, maxfhist);

            % Preconditions
            if consts_obj.DEBUGGING
                debug_obj.assert(abs(iprint) <= 3, "IPRINT is 0, 1, -1, 2, -2, 3, or -3", srname);
                debug_obj.assert(n >= 1 && npt == (n + 1) * (n + 2) / 2, "N >= 1, NPT == (N+1)*(N+2)/2", srname);
                debug_obj.assert(maxfun >= npt + 1, "MAXFUN >= NPT + 1", srname);
                debug_obj.assert(maxhist >= 0 && maxhist <= maxfun, "0 <= MAXHIST <= MAXFUN", srname);
                debug_obj.assert(maxfhist * (maxfhist - maxhist) == 0, "SIZE(FHIST) == 0 or MAXHIST", srname);
                debug_obj.assert(numel(fval) == npt, "SIZE(FVAL) == NPT", srname);
                debug_obj.assert(size(xhist, 1) == n && maxxhist * (maxxhist - maxhist) == 0, "SIZE(XHIST, 1) == N, SIZE(XHIST, 2) == 0 or MAXHIST", srname);
                debug_obj.assert(rhobeg > 0, "RHOBEG > 0", srname);
                debug_obj.assert(numel(x0) == n && all(infnan_obj.is_finite(x0), 'all'), "SIZE(X0) == N, X0 is finite", srname);
                debug_obj.assert(numel(xbase) == n, "SIZE(XBASE) == N", srname);
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
            % is not fully parallelizable if NPT>2N+1, as the definition XPT(:, 2N+2:end) involves FVAL(1:2N+1).
            evaluated(:) = false;

            % Initialize XHIST, FHIST, and FVAL. Otherwise, compilers may complain that they are not
            % (completely) initialized if the initialization aborts due to abnormality (see CHECKEXIT).
            % N.B.: 1. Initializing them to NaN would be more reasonable (NaN is not available in Fortran).
            % 2. Do not initialize the models if the current initialization aborts due to abnormality. Otherwise,
            % errors or exceptions may occur, as FVAL and XPT etc are uninitialized.
            xhist = repmat(-consts_obj.REALMAX, size(xhist));
            fhist(:) = consts_obj.REALMAX;
            fval(:) = consts_obj.REALMAX;

            % Set XPT(:, 1 : 2*N+1) and FVAL(:, 1 : 2*N+1).
            xpt = repmat(consts_obj.ZERO, size(xpt));
            kk(:) = linalg_obj.linspace_i(2, 2 * n, n);
            xpt(:, kk) = rhobeg * linalg_obj.eye1(n);
            for k = 1:2 * n + 1
                x(:) = xpt(:, k) + xbase;
                f = evaluate_obj.evaluatef(calfun, x);

                % Print a message about the function evaluation according to IPRINT.
                message_obj.fmsg(solver, "Initialization", iprint, k, rhobeg, f, x);
                % Save X and F into the history.
                [xhist, fhist] = history_obj.savehist(k, x, xhist, f, fhist);

                evaluated(k) = true;
                fval(k) = f;

                % When K is even, determine XPT(:, K+1) according to F(K).
                % N.B.: This heuristic strategy does increase the performance. In addition, it is a GOOD idea to
                % evaluate F at XBASE + 2*XPT(:, K) or XBASE - XPT(:, K) IMMEDIATELY after XBASE + XPT(:, K).
                % This increases the probability of finding a smaller function value earlier in the sampling
                % process for the first model. This process itself can be regarded as a simple direct search. It
                % is IMPORTANT for the performance of the algorithm during the early stage, even before the
                % first model is built.
                if mod(k, 2) == 0
                    if fval(k) < fval(1)
                        xpt(:, k + 1) = consts_obj.TWO * xpt(:, k); % XPT(K / 2, K + 1) = TWO * RHOBEG

                    else
                        xpt(:, k + 1) = -xpt(:, k); % XPT(K / 2, K + 1) = -RHOBEG
                    end
                end

                % Check whether to exit.
                subinfo = checkexit_obj.checkexit_unc(maxfun, k, f, ftarget, x);
                if subinfo ~= infos_obj.INFO_DFT
                    info = subinfo;
                    break
                end
            end

            if info == infos_obj.INFO_DFT
                xw(:) = -rhobeg;
                xw(linalg_obj.trueloc(fval(kk) < fval(1))) = rhobeg;
                % See (42)--(43) of the UOBYQA paper for IP and IQ.
                ip = 0;
                iq = 2;
                for k = 2 * n + 2:npt
                    % Pick the shift from XBASE to the next initial interpolation point that provides the
                    % off-diagonal second derivatives of the quadratic interpolant.
                    ip = ip + 1;
                    if ip == iq
                        iq = iq + 1;
                        ip = 1;
                    end
                    xpt([ip, iq], k) = xw([ip, iq]);
                    x(:) = xpt(:, k) + xbase;
                    f = evaluate_obj.evaluatef(calfun, x);

                    % Print a message about the function evaluation according to IPRINT.
                    message_obj.fmsg(solver, "Initialization", iprint, k, rhobeg, f, x);
                    % Save X and F into the history.
                    [xhist, fhist] = history_obj.savehist(k, x, xhist, f, fhist);

                    evaluated(k) = true;
                    fval(k) = f;

                    % Check whether to exit.
                    subinfo = checkexit_obj.checkexit_unc(maxfun, k, f, ftarget, x);
                    if subinfo ~= infos_obj.INFO_DFT
                        info = subinfo;
                        break
                    end
                end
            end

            nf = fix(nnz(evaluated)); %%MATLAB: nf = sum(evaluated);
            kopt = fix(fortran.minloc(fval, 'mask', evaluated, 'dim', 1));
            %%MATLAB: fopt = min(fval(evaluated)); kopt = find(evaluated & ~(fval > fopt), 1, 'first')

            %====================%
            %  Calculation ends  %
            %====================%

            % Postconditions
            if consts_obj.DEBUGGING
                debug_obj.assert(nf <= npt, "NF <= NPT", srname);
                debug_obj.assert(kopt >= 1 && kopt <= npt, "1 <= KOPT <= NPT", srname);
                debug_obj.assert(numel(xbase) == n && all(infnan_obj.is_finite(xbase), 'all'), "SIZE(XBASE) == N, XBASE is finite", srname);
                debug_obj.assert(size(xpt, 1) == n && size(xpt, 2) == npt, "SIZE(XPT) == [N, NPT]", srname);
                debug_obj.assert(all(infnan_obj.is_finite(xpt), 'all'), "XPT is finite", srname);
                debug_obj.assert(numel(fval) == npt && ~any(evaluated & (infnan_obj.is_nan(fval) | infnan_obj.is_posinf(fval)), 'all'), "SIZE(FVAL) == NPT and FVAL is not NaN or +Inf", srname);
                debug_obj.assert(~any(evaluated & fval < fval(kopt), 'all'), "FVAL(KOPT) = MINVAL(FVAL)", srname);
                debug_obj.assert(numel(fhist) == maxfhist, "SIZE(FHIST) == MAXFHIST", srname);
                debug_obj.assert(size(xhist, 1) == n && size(xhist, 2) == maxxhist, "SIZE(XHIST) == [N, MAXXHIST]", srname);
            end

        end
        function [pq, info] = initq(~, fval, xpt, pq, varargin)
            %--------------------------------------------------------------------------------------------------%
            % This subroutine initializes the quadratic model, whose coefficients are stored in PQ, where
            % PQ(1 : N) containing the gradient of the model at XBASE, and PQ(N+1 : NPT-1) containing the upper
            % triangular part of the Hessian, column by column. See Section 4 of the UOBYQA paper.
            %--------------------------------------------------------------------------------------------------%

            % Common modules
            consts_obj = prima_mat.common.consts_mod();
            debug_obj = prima_mat.common.debug_mod();
            infnan_obj = prima_mat.common.infnan_mod();
            infos_obj = prima_mat.common.infos_mod();


            % Inputs
            % XPT(N, NPT)
            % XPT(N, NPT)

            % Outputs

            % PQ((N + 1) * (N + 2) / 2 - 1)

            % Local variables
            srname = "INITQ";
            k1 = NaN;
            ih = NaN;
            ip = NaN;
            iq = NaN;
            k = NaN;
            k0 = NaN;
            n = NaN;
            npt = NaN;
            deriv = NaN(size(xpt, 1), 1);
            fbase = NaN;
            rhobeg = NaN;
            rhosq = NaN;

            % Sizes
            n = size(xpt, 1);
            npt = size(xpt, 2);

            % Postconditions
            if consts_obj.DEBUGGING
                debug_obj.assert(n >= 1 && npt == (n + 1) * (n + 2) / 2, "N >= 1, NPT == (N+1)*(N+2)/2", srname);
                debug_obj.assert(size(xpt, 1) == n && size(xpt, 2) == npt, "SIZE(XPT) == [N, NPT]", srname);
                debug_obj.assert(all(infnan_obj.is_finite(xpt), 'all'), "XPT is finite", srname);
                debug_obj.assert(numel(fval) == npt && ~any((infnan_obj.is_nan(fval) | infnan_obj.is_posinf(fval)), 'all'), "SIZE(FVAL) == NPT and FVAL is not NaN or +Inf", srname);
            end

            %====================%
            % Calculation starts %
            %====================%

            rhobeg = max(abs(xpt(:, 2)), [], 'all');
            rhosq = fortran.power(rhobeg, 2);
            fbase = fval(1);

            % Form the gradient and diagonal second derivatives of the quadratic model.
            for k = 1:n
                k0 = 2 * k;
                k1 = 2 * k + 1;
                % Find the (K, K) element of the Hessian.
                ih = n + k * (k + 1) / 2;
                if xpt(k, k1) > 0
                    % XPT(K, K1) = 2*RHO
                    deriv(k) = (fbase + fval(k1) - consts_obj.TWO * fval(k0)) / rhosq;
                    pq(k) = (4.0 * fval(k0) - 3.0 * fbase - fval(k1)) / (consts_obj.TWO * rhobeg);
                else                    % XPT(K, K1) = -RHO
                    deriv(k) = (fval(k0) + fval(k1) - consts_obj.TWO * fbase) / rhosq;
                    pq(k) = (fval(k0) - fval(k1)) / (consts_obj.TWO * rhobeg);
                end
                pq(ih) = deriv(k);
            end

            % Form the off-diagonal second derivatives of the initial quadratic model.
            ip = 0;
            iq = 2;
            for k = 2 * n + 2:npt
                ip = ip + 1;
                if ip == iq
                    iq = iq + 1;
                    ip = 1;
                end
                % Find the (IQ, IP) entry of the Hessian.
                ih = n + (iq - 1) * iq / 2 + ip;
                pq(ih) = (fval(k) - fbase - xpt(ip, k) * pq(ip) - xpt(iq, k) * pq(iq) - consts_obj.HALF * rhosq * (deriv(ip) + deriv(iq))) / (xpt(ip, k) * xpt(iq, k));
            end

            ipObj = inputParser();
            addParameter(ipObj, 'info', NaN);
            parse(ipObj, varargin{:});
            info = ipObj.Results.info;
            if nargout >= 2
                if any(infnan_obj.is_nan(pq), 'all')
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
                debug_obj.assert(numel(pq) == npt - 1, "SIZE(PQ) == NPT - 1", srname);
            end

        end
        function [pl, info] = initl(~, xpt, pl, varargin)
            %--------------------------------------------------------------------------------------------------%
            % This subroutine initializes the Lagrange functions. The coefficients of the K-th Lagrange function
            % is stored in PL(:, K), with PL(1 : N, K) containing the gradient of the function at XBASE, and
            % PL(N+1 : NPT-1, K) containing the upper triangular part of the Hessian, column by column.
            % See Section 4 of the UOBYQA paper.
            %--------------------------------------------------------------------------------------------------%
            % Common modules
            consts_obj = prima_mat.common.consts_mod();
            debug_obj = prima_mat.common.debug_mod();
            infnan_obj = prima_mat.common.infnan_mod();
            infos_obj = prima_mat.common.infos_mod();


            % Inputs
            % XPT(N, NPT)

            % Outputs

            % PL((N + 1) * (N + 2) / 2 - 1, (N + 1) * (N + 2) / 2)

            % Local variables
            srname = "INITL";
            ih = NaN;
            ip = NaN;
            iq = NaN;
            k = NaN;
            k0 = NaN;
            k1 = NaN;
            n = NaN;
            npt = NaN;
            rhobeg = NaN;
            rhosq = NaN;
            temp = NaN;

            % Sizes
            n = size(xpt, 1);
            npt = size(xpt, 2);

            % Postconditions
            if consts_obj.DEBUGGING
                debug_obj.assert(n >= 1 && npt == (n + 1) * (n + 2) / 2, "N >= 1, NPT == (N+1)*(N+2)/2", srname);
                debug_obj.assert(size(xpt, 1) == n && size(xpt, 2) == npt, "SIZE(XPT) == [N, NPT]", srname);
                debug_obj.assert(all(infnan_obj.is_finite(xpt), 'all'), "XPT is finite", srname);
            end

            %====================%
            % Calculation starts %
            %====================%

            rhobeg = max(abs(xpt(:, 2)), [], 'all');
            rhosq = fortran.power(rhobeg, 2);

            pl = repmat(consts_obj.ZERO, size(pl));

            % Form the gradient and diagonal second derivatives of the Lagrange functions.
            for k = 1:n
                k0 = 2 * k;
                k1 = 2 * k + 1;
                ih = n + k * (k + 1) / 2; % The (K, K) entry of the Hessian
                if xpt(k, k1) > 0
                    % XPT(K, K1) = 2*RHO
                    pl(k, 1) = -1.5 / rhobeg;
                    pl(ih, 1) = consts_obj.ONE / rhosq;
                    pl(k, k0) = consts_obj.TWO / rhobeg;
                    pl(ih, k0) = -consts_obj.TWO / rhosq;
                else                    % XPT(K, K1) = -RHO
                    pl(ih, 1) = -consts_obj.TWO / rhosq;
                    pl(k, k0) = consts_obj.HALF / rhobeg;
                    pl(ih, k0) = consts_obj.ONE / rhosq;
                end
                pl(k, k1) = -consts_obj.HALF / rhobeg;
                pl(ih, k1) = consts_obj.ONE / rhosq;
            end

            % Form the off-diagonal second derivatives of the Lagrange functions.
            ip = 0;
            iq = 2;
            for k = 2 * n + 2:npt
                ip = ip + 1;
                if ip == iq
                    iq = iq + 1;
                    ip = 1;
                end

                % Find the (IQ, IP) entry of the Hessian.
                temp = consts_obj.ONE / (xpt(ip, k) * xpt(iq, k));
                ih = n + (iq - 1) * iq / 2 + ip;

                pl(ih, 1) = temp;
                pl(ih, k) = temp;

                if xpt(ip, k) < 0
                    pl(ih, 2 * ip + 1) = -temp;
                else
                    pl(ih, 2 * ip) = -temp;
                end

                if xpt(iq, k) < 0
                    pl(ih, 2 * iq + 1) = -temp;
                else
                    pl(ih, 2 * iq) = -temp;
                end
            end

            ipObj = inputParser();
            addParameter(ipObj, 'info', NaN);
            parse(ipObj, varargin{:});
            info = ipObj.Results.info;
            if nargout >= 2
                if any(infnan_obj.is_nan(pl), 'all')
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
                debug_obj.assert(size(pl, 1) == npt - 1 && size(pl, 2) == npt, "SIZE(PL) == [NPT - 1, NPT]", srname);
            end

        end

    end
end