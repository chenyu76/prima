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


            checkexit_obj = prima_mat.common.checkexit_mod();


            evaluate_obj = prima_mat.common.evaluate_mod();
            history_obj = prima_mat.common.history_mod();


            message_obj = prima_mat.common.message_mod();


            % N.B.: INTENT cannot be specified if a dummy procedure is not a POINTER



            % X0(N)



            % FHIST(MAXFHIST)

            % XBASE(N)
            % XHIST(N, MAXXHIST)
            % XPT(N, NPT)


            solver = "UOBYQA";


            kk = NaN(numel(x0), 1);


            evaluated = false(size(xpt, 2), 1);

            x = NaN(numel(x0), 1);
            xw = NaN(numel(x0), 1);


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
            % is not fully parallelizable if NPT>2N+1, as the definition XPT(:, 2N+2:end) involves FVAL(1:2N+1).
            evaluated(:) = false;

            % Initialize XHIST, FHIST, and FVAL. Otherwise, compilers may complain that they are not
            % (completely) initialized if the initialization aborts due to abnormality (see CHECKEXIT).
            % N.B.: 1. Initializing them to NaN would be more reasonable (NaN is not available in Fortran).
            % 2. Do not initialize the models if the current initialization aborts due to abnormality. Otherwise,
            % errors or exceptions may occur, as FVAL and XPT etc are uninitialized.
            xhist = repmat(-realmax, size(xhist));
            fhist(:) = realmax;
            fval(:) = realmax;

            % Set XPT(:, 1 : 2*N+1) and FVAL(:, 1 : 2*N+1).
            xpt = zeros(size(xpt));
            kk(:) = linspace(2, 2 * n, n).';
            xpt(:, kk) = rhobeg * eye(n);
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
                        xpt(:, k + 1) = 2.0 * xpt(:, k); % XPT(K / 2, K + 1) = TWO * RHOBEG

                    else
                        xpt(:, k + 1) = -xpt(:, k); % XPT(K / 2, K + 1) = -RHOBEG
                    end
                end

                % Check whether to exit.
                subinfo = checkexit_obj.checkexit_unc(maxfun, k, f, ftarget, x);
                if subinfo ~= 0
                    info = subinfo;
                    break
                end
            end

            if info == 0
                xw(:) = -rhobeg;
                xw(fval(kk) < fval(1)) = rhobeg;
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
                    if subinfo ~= 0
                        info = subinfo;
                        break
                    end
                end
            end

            nf = nnz(evaluated); %%MATLAB: nf = sum(evaluated);
            kopt = fortran.minloc(fval, 'mask', evaluated, 'dim', 1);
            %%MATLAB: fopt = min(fval(evaluated)); kopt = find(evaluated & ~(fval > fopt), 1, 'first')

            %====================%
            %  Calculation ends  %
            %====================%



        end
        function [pq, info] = initq(~, fval, xpt, pq, varargin)
            %--------------------------------------------------------------------------------------------------%
            % This subroutine initializes the quadratic model, whose coefficients are stored in PQ, where
            % PQ(1 : N) containing the gradient of the model at XBASE, and PQ(N+1 : NPT-1) containing the upper
            % triangular part of the Hessian, column by column. See Section 4 of the UOBYQA paper.
            %--------------------------------------------------------------------------------------------------%



            % XPT(N, NPT)
            % XPT(N, NPT)



            % PQ((N + 1) * (N + 2) / 2 - 1)



            deriv = NaN(size(xpt, 1), 1);


            n = size(xpt, 1);
            npt = size(xpt, 2);


            %====================%
            % Calculation starts %
            %====================%

            rhobeg = max(abs(xpt(:, 2)), [], 'all');
            rhosq = rhobeg ^ 2;
            fbase = fval(1);

            % Form the gradient and diagonal second derivatives of the quadratic model.
            for k = 1:n
                k0 = 2 * k;
                k1 = 2 * k + 1;
                % Find the (K, K) element of the Hessian.
                ih = n + k * (k + 1) / 2;
                if xpt(k, k1) > 0
                    % XPT(K, K1) = 2*RHO
                    deriv(k) = (fbase + fval(k1) - 2.0 * fval(k0)) / rhosq;
                    pq(k) = (4.0 * fval(k0) - 3.0 * fbase - fval(k1)) / (2.0 * rhobeg);
                else                    % XPT(K, K1) = -RHO
                    deriv(k) = (fval(k0) + fval(k1) - 2.0 * fbase) / rhosq;
                    pq(k) = (fval(k0) - fval(k1)) / (2.0 * rhobeg);
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
                pq(ih) = (fval(k) - fbase - xpt(ip, k) * pq(ip) - xpt(iq, k) * pq(iq) - 0.5 * rhosq * (deriv(ip) + deriv(iq))) / (xpt(ip, k) * xpt(iq, k));
            end

            ipObj = inputParser();
            addParameter(ipObj, 'info', NaN);
            parse(ipObj, varargin{:});
            info = ipObj.Results.info;
            if nargout >= 2
                if any(isnan(pq), 'all')
                    info = -3;
                else
                    info = 0;
                end
            end

            %====================%
            %  Calculation ends  %
            %====================%



        end
        function [pl, info] = initl(~, xpt, pl, varargin)
            %--------------------------------------------------------------------------------------------------%
            % This subroutine initializes the Lagrange functions. The coefficients of the K-th Lagrange function
            % is stored in PL(:, K), with PL(1 : N, K) containing the gradient of the function at XBASE, and
            % PL(N+1 : NPT-1, K) containing the upper triangular part of the Hessian, column by column.
            % See Section 4 of the UOBYQA paper.
            %--------------------------------------------------------------------------------------------------%



            % XPT(N, NPT)



            % PL((N + 1) * (N + 2) / 2 - 1, (N + 1) * (N + 2) / 2)



            n = size(xpt, 1);
            npt = size(xpt, 2);


            %====================%
            % Calculation starts %
            %====================%

            rhobeg = max(abs(xpt(:, 2)), [], 'all');
            rhosq = rhobeg ^ 2;

            pl = zeros(size(pl));

            % Form the gradient and diagonal second derivatives of the Lagrange functions.
            for k = 1:n
                k0 = 2 * k;
                k1 = 2 * k + 1;
                ih = n + k * (k + 1) / 2; % The (K, K) entry of the Hessian
                if xpt(k, k1) > 0
                    % XPT(K, K1) = 2*RHO
                    pl(k, 1) = -1.5 / rhobeg;
                    pl(ih, 1) = 1.0 / rhosq;
                    pl(k, k0) = 2.0 / rhobeg;
                    pl(ih, k0) = -2.0 / rhosq;
                else                    % XPT(K, K1) = -RHO
                    pl(ih, 1) = -2.0 / rhosq;
                    pl(k, k0) = 0.5 / rhobeg;
                    pl(ih, k0) = 1.0 / rhosq;
                end
                pl(k, k1) = -0.5 / rhobeg;
                pl(ih, k1) = 1.0 / rhosq;
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
                temp = 1.0 / (xpt(ip, k) * xpt(iq, k));
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
                if any(isnan(pl), 'all')
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