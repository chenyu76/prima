classdef initialize_cobyla_mod
    %--------------------------------------------------------------------------------------------------%
    % This module contains subroutines for initialization.
    %
    % Coded by Zaikun ZHANG (www.zhangzk.net) based on Powell's code and the COBYLA paper.
    %
    % Dedicated to the late Professor M. J. D. Powell FRS (1936--2015).
    %
    % Started: July 2021
    %
    % Last Modified: Tue 16 Sep 2025 12:35:23 PM CST
    %--------------------------------------------------------------------------------------------------%

    methods
        function [nf, chist, conhist, conmat, cval, fhist, fval, sim, simi, xhist, evaluated, info] = initxfc(~, calcfc, iprint, maxfun, amat, bvec, constr0, ctol, f0, ftarget, rhobeg, x0, chist, conhist, conmat, cval, fhist, fval, sim, simi, xhist, evaluated)
            %--------------------------------------------------------------------------------------------------%
            % This subroutine does the initialization concerning X, function values, and constraints.
            %--------------------------------------------------------------------------------------------------%

            % Common modules
            checkexit_obj = prima_mat.common.checkexit_mod();


            evaluate_obj = prima_mat.common.evaluate_mod();
            history_obj = prima_mat.common.history_mod();


            message_obj = prima_mat.common.message_mod();


            % Inputs
            % N.B.: INTENT cannot be specified if a dummy procedure is not a POINTER


            % AMAT(N, M_LCON)
            % BVEC(M_LCON)
            % CONSTR0(M)



            % X0(N)

            % Outputs


            % EVALUATED(N+1)
            % CHIST(MAXCHIST)
            % CONHIST(M, MAXCONHIST)
            % CONMAT(M, N+1)
            % CVAL(N+1)
            % FHIST(MAXFHIST)
            % FVAL(N+1)
            % SIM(N, N+1)
            % SIMI(N, N)
            % XHIST(N, MAXXHIST)

            % Local variables
            solver = "COBYLA";


            constr = NaN(size(conmat, 1), 1);


            x = NaN(numel(x0), 1);


            % Sizes
            m_lcon = numel(bvec);
            m = size(conmat, 1);
            n = size(sim, 1);


            % Preconditions


            %====================%
            % Calculation starts %
            %====================%

            % Initialize INFO to the default value. At return, an INFO different from this value will indicate
            % an abnormal return.
            info = 0;

            % Initialize the simplex. It will be revised during the initialization.
            sim(:, :) = eye(n, n + 1) * rhobeg;
            sim(:, n + 1) = x0;

            % Initialize the matrix SIMI. This initial value will be discarded at the end of the initialization.
            % If we do not do this, compilers may complain if we return due to CHECKEXIT before SIMI is set.
            simi(:, :) = eye(n) ./ rhobeg;

            % EVALUATED(J) = TRUE iff the function/constraint of SIM(:, J) has been evaluated.
            evaluated(:) = false;

            % Initialize XHIST, FHIST, CHIST, CONHIST, FVAL, CVAL, and CONMAT. Otherwise, compilers may complain
            %that they are not (completely) initialized if the initialization aborts due to abnormality (see
            %CHECKEXIT).
            % N.B.: 1. Initializing them to NaN would be more reasonable (NaN is not available in Fortran).
            % 2. Do not initialize the models if the current initialization aborts due to abnormality. Otherwise,
            % errors or exceptions may occur, as FVAL and XPT etc are uninitialized.
            xhist = repmat(-realmax, size(xhist));
            fhist(:) = realmax;
            chist(:) = realmax;
            conhist = repmat(realmax, size(conhist));
            fval(:) = realmax;
            cval(:) = realmax;
            conmat = repmat(realmax, size(conmat));

            for k = 1:n + 1
                x(:) = sim(:, n + 1);
                % We will evaluate F corresponding to SIM(:, J).
                if k == 1
                    j = n + 1;
                    f = f0;
                    constr(:) = constr0;
                else
                    j = k - 1;
                    x(j) = x(j) + rhobeg;
                    constr(1:m_lcon) = evaluate_obj.moderatec(amat.' * x - bvec); % Linear constraints.
                    [f, constr_slice] = evaluate_obj.evaluatefc(calcfc, x, constr(m_lcon + 1:m)); constr(m_lcon + 1:m) = constr_slice; % Nonlinear constraints.
                    % Note that EVALUATE moderates the nonlinear constraint values. Thus we also moderate the
                    % linear constraint values here to make CSTRV consistent.
                end
                cstrv = max([0.0; constr], [], 'all');

                % Print a message about the function/constraint evaluation according to IPRINT.
                message_obj.fmsg(solver, "Initialization", iprint, k, rhobeg, f, x, 'cstrv', cstrv, 'constr', constr);
                % Save X, F, CONSTR, CSTRV into the history.
                [xhist, fhist, chist, conhist] = history_obj.savehist(k, x, xhist, f, fhist, 'cstrv', cstrv, 'chist', chist, 'constr', constr, 'conhist', conhist);

                % Save F, CONSTR, and CSTRV to FVAL, CONMAT, and CVAL respectively. This must be done before
                % checking whether to exit. If exit, FVAL, CONMAT, and CVAL will define FFILT, CONFILT, and
                % CFILT, which will define the returned X, F, CONSTR, and CSTRV.
                evaluated(j) = true;
                fval(j) = f;
                conmat(:, j) = constr;
                cval(j) = cstrv;

                % Check whether to exit.
                subinfo = checkexit_obj.checkexit_con(maxfun, k, cstrv, ctol, f, ftarget, x);
                if subinfo ~= 0
                    info = subinfo;
                    break
                end

                % Exchange the new vertex of the initial simplex with the optimal vertex if necessary.
                % This is the ONLY part that is essentially non-parallel.
                if j <= n && fval(j) < fval(n + 1)
                    fval([j, n + 1]) = fval([n + 1, j]);
                    cval([j, n + 1]) = cval([n + 1, j]);
                    conmat(:, [j, n + 1]) = conmat(:, [n + 1, j]);
                    sim(:, n + 1) = x;
                    sim(j, 1:j) = -rhobeg; % SIM(:, 1:N) is lower triangular.

                end
            end

            nf = nnz(evaluated);

            if all(evaluated, 'all')
                % Initialize SIMI to the inverse of SIM(:, 1:N).
                simi(:, :) = inv(sim(:, 1:n));
            end

            %====================%
            %  Calculation ends  %
            %====================%

            % Postconditions


        end
        function [nfilt, cfilt, confilt, ffilt, xfilt] = initfilt(~, conmat, ctol, cweight, cval, fval, sim, evaluated, cfilt, confilt, ffilt, xfilt)
            %--------------------------------------------------------------------------------------------------%
            % This subroutine initializes the filters (XFILT, etc) that will be used when selecting X at the
            % end of the solver.
            % N.B.:
            % 1. Why not initialize the filters using XHIST, etc? Because the history is empty if the user
            % chooses not to output it.
            % 2. We decouple INITXFC and INITFILT so that it is easier to parallelize the former if needed.
            %--------------------------------------------------------------------------------------------------%

            % Common modules



            selectx_obj = prima_mat.common.selectx_mod();

            % Inputs



            % In-outputs



            % Local variables



            x = NaN(size(sim, 1), 1);

            % Sizes

            n = size(sim, 1);


            % Preconditions


            %====================%
            % Calculation starts %
            %====================%

            nfilt = 0;
            for i = 1:n + 1
                if evaluated(i)
                    if i <= n
                        x(:) = sim(:, i) + sim(:, n + 1);
                    else
                        x(:) = sim(:, i); % I == N+1
                    end
                    [nfilt, cfilt, ffilt, xfilt, confilt] = selectx_obj.savefilt(cval(i), ctol, cweight, fval(i), x, nfilt, cfilt, ffilt, xfilt, 'constr', conmat(:, i), 'confilt', confilt);
                end
            end

            %====================%
            %  Calculation ends  %
            %====================%

            % Postconditions

        end

    end
end