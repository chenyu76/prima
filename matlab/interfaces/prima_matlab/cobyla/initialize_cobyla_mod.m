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
            checkexit_obj = checkexit_mod();
            consts_obj = consts_mod();
            debug_obj = debug_mod();
            evaluate_obj = evaluate_mod();
            history_obj = history_mod();
            infnan_obj = infnan_mod();
            infos_obj = infos_mod();
            linalg_obj = linalg_mod();
            message_obj = message_mod();
            pintrf_obj = pintrf_mod();


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
            srname = "INITIALIZE";
            j = NaN;
            k = NaN;
            m = NaN;
            m_lcon = NaN;
            maxchist = NaN;
            maxconhist = NaN;
            maxfhist = NaN;
            maxhist = NaN;
            maxxhist = NaN;
            n = NaN;
            subinfo = NaN;
            constr = NaN(size(conmat, 1), 1);
            cstrv = NaN;
            f = NaN;
            x = NaN(numel(x0), 1);
            itol = consts_obj.TENTH;

            % Sizes
            m_lcon = fix(numel(bvec));
            m = size(conmat, 1);
            n = size(sim, 1);
            maxchist = fix(numel(chist));
            maxconhist = size(conhist, 2);
            maxfhist = fix(numel(fhist));
            maxxhist = size(xhist, 2);
            maxhist = fix(max(maxchist, max(maxconhist, max(maxfhist, maxxhist))));

            % Preconditions
            if consts_obj.DEBUGGING
                debug_obj.assert(m >= 0, "M >= 0", srname);
                debug_obj.assert(n >= 1, "N >= 1", srname);
                debug_obj.assert(abs(iprint) <= 3, "IPRINT is 0, 1, -1, 2, -2, 3, or -3", srname);
                debug_obj.assert(size(amat, 1) == n && size(amat, 2) == numel(bvec), "SIZE(AMAT) == [N, SIZE(BVEC)]", srname);
                debug_obj.assert(size(conmat, 1) == m && size(conmat, 2) == n + 1, "SIZE(CONMAT) = [M, N+1]", srname);
                debug_obj.assert(numel(cval) == n + 1, "SIZE(CVAL) == N+1", srname);
                debug_obj.assert(numel(fval) == n + 1, "SIZE(FVAL) == N+1", srname);
                debug_obj.assert(size(sim, 1) == n && size(sim, 2) == n + 1, "SIZE(SIM) == [N, N+1]", srname);
                debug_obj.assert(size(simi, 1) == n && size(simi, 2) == n, "SIZE(SIMI) == [N, N]", srname);
                debug_obj.assert(numel(evaluated) == n + 1, "SIZE(EVALUATED) == N + 1", srname);
                debug_obj.assert(maxchist * (maxchist - maxhist) == 0, "SIZE(CHIST) == 0 or MAXHIST", srname);
                debug_obj.assert(size(conhist, 1) == m && maxconhist * (maxconhist - maxhist) == 0, "SIZE(CONHIST, 1) == M, SIZE(CONHIST, 2) == 0 or MAXHIST", srname);
                debug_obj.assert(maxfhist * (maxfhist - maxhist) == 0, "SIZE(FHIST) == 0 or MAXHIST", srname);
                debug_obj.assert(size(xhist, 1) == n && maxxhist * (maxxhist - maxhist) == 0, "SIZE(XHIST, 1) == N, SIZE(XHIST, 2) == 0 or MAXHIST", srname);
                debug_obj.assert(numel(x0) == n && all(infnan_obj.is_finite(x0), 'all'), "SIZE(X0) == N, X0 is finite", srname);
                debug_obj.assert(rhobeg > 0, "RHOBEG > 0", srname);
            end

            %====================%
            % Calculation starts %
            %====================%

            % Initialize INFO to the default value. At return, an INFO different from this value will indicate
            % an abnormal return.
            info = infos_obj.INFO_DFT;

            % Initialize the simplex. It will be revised during the initialization.
            sim(:, :) = linalg_obj.eye2(n, n + 1) * rhobeg;
            sim(:, n + 1) = x0;

            % Initialize the matrix SIMI. This initial value will be discarded at the end of the initialization.
            % If we do not do this, compilers may complain if we return due to CHECKEXIT before SIMI is set.
            simi(:, :) = linalg_obj.eye1(n) ./ rhobeg;

            % EVALUATED(J) = TRUE iff the function/constraint of SIM(:, J) has been evaluated.
            evaluated = repmat(false, size(evaluated));

            % Initialize XHIST, FHIST, CHIST, CONHIST, FVAL, CVAL, and CONMAT. Otherwise, compilers may complain
            %that they are not (completely) initialized if the initialization aborts due to abnormality (see
            %CHECKEXIT).
            % N.B.: 1. Initializing them to NaN would be more reasonable (NaN is not available in Fortran).
            % 2. Do not initialize the models if the current initialization aborts due to abnormality. Otherwise,
            % errors or exceptions may occur, as FVAL and XPT etc are uninitialized.
            xhist = repmat(-consts_obj.REALMAX, size(xhist));
            fhist = repmat(consts_obj.REALMAX, size(fhist));
            chist = repmat(consts_obj.REALMAX, size(chist));
            conhist = repmat(consts_obj.REALMAX, size(conhist));
            fval = repmat(consts_obj.REALMAX, size(fval));
            cval = repmat(consts_obj.REALMAX, size(cval));
            conmat = repmat(consts_obj.REALMAX, size(conmat));

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
                    constr(1:m_lcon) = evaluate_obj.moderatec(linalg_obj.matprod12(x, amat) - bvec); % Linear constraints.
                    [f, constr_slice] = evaluate_obj.evaluatefc(calcfc, x, constr(m_lcon + 1:m)); constr(m_lcon + 1:m) = constr_slice; % Nonlinear constraints.
                    % Note that EVALUATE moderates the nonlinear constraint values. Thus we also moderate the
                    % linear constraint values here to make CSTRV consistent.
                end
                cstrv = linalg_obj.maximum1([consts_obj.ZERO; reshape(constr, [], 1)]);

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
                if subinfo ~= infos_obj.INFO_DFT
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

            nf = fix(nnz(evaluated));

            if all(evaluated, 'all')
                % Initialize SIMI to the inverse of SIM(:, 1:N).
                simi(:, :) = linalg_obj.inv(sim(:, 1:n));
            end

            %====================%
            %  Calculation ends  %
            %====================%

            % Postconditions
            if consts_obj.DEBUGGING
                debug_obj.assert(nf <= maxfun, "NF <= MAXFUN", srname);
                debug_obj.assert(numel(evaluated) == n + 1, "SIZE(EVALUATED) == N + 1", srname);
                debug_obj.assert(numel(chist) == maxchist, "SIZE(CHIST) == MAXCHIST", srname);
                debug_obj.assert(~any(chist(1:min(nf, maxchist)) < 0 | infnan_obj.is_nan(chist(1:min(nf, maxchist))) | infnan_obj.is_posinf(chist(1:min(nf, maxchist))), 'all'), "CHIST does not contain negative values or NaN/+Inf", srname);
                debug_obj.assert(size(conhist, 1) == m && size(conhist, 2) == maxconhist, "SIZE(CONHIST) == [M, MAXCONHIST]", srname);
                debug_obj.assert(~any(infnan_obj.is_nan(conhist(:, 1:min(nf, maxconhist))) | infnan_obj.is_posinf(conhist(:, 1:min(nf, maxconhist))), 'all'), "CONHIST does not contain NaN/+Inf", srname);
                debug_obj.assert(size(conmat, 1) == m && size(conmat, 2) == n + 1, "SIZE(CONMAT) = [M, N+1]", srname);
                debug_obj.assert(~any(infnan_obj.is_nan(conmat) | infnan_obj.is_posinf(conmat), 'all'), "CONMAT does not contain NaN/+Inf", srname);
                debug_obj.assert(numel(cval) == n + 1 && ~any(cval < 0 | infnan_obj.is_nan(cval) | infnan_obj.is_posinf(cval), 'all'), "SIZE(CVAL) == N+1 and CVAL does not contain negative values or NaN/+Inf", srname);
                debug_obj.assert(numel(fhist) == maxfhist, "SIZE(FHIST) == MAXFHIST", srname);
                debug_obj.assert(maxfhist * (maxfhist - maxhist) == 0, "SIZE(FHIST) == 0 or MAXHIST", srname);
                debug_obj.assert(~any(infnan_obj.is_nan(fhist(1:min(nf, maxfhist))) | infnan_obj.is_posinf(fhist(1:min(nf, maxfhist))), 'all'), "FHIST does not contain NaN/+Inf", srname);
                debug_obj.assert(numel(fval) == n + 1 && ~any(infnan_obj.is_nan(fval) | infnan_obj.is_posinf(fval), 'all'), "SIZE(FVAL) == N+1 and FVAL does not contain NaN/+Inf", srname);
                debug_obj.assert(size(xhist, 1) == n && size(xhist, 2) == maxxhist, "SIZE(XHIST) == [N, MAXXHIST]", srname);
                debug_obj.assert(~any(infnan_obj.is_nan(xhist(:, 1:min(nf, maxxhist))), 'all'), "XHIST does not contain NaN", srname);
                debug_obj.assert(size(sim, 1) == n && size(sim, 2) == n + 1, "SIZE(SIM) == [N, N+1]", srname);
                debug_obj.assert(all(infnan_obj.is_finite(sim), 'all'), "SIM is finite", srname);
                debug_obj.assert(all(sum(abs(sim(:, 1:n)), 1).' > 0, 'all'), "SIM(:, 1:N) has no zero column", srname);
                debug_obj.assert(size(simi, 1) == n && size(simi, 2) == n, "SIZE(SIMI) == [N, N]", srname);
                debug_obj.assert(all(infnan_obj.is_finite(simi), 'all'), "SIMI is finite", srname);
                debug_obj.assert(linalg_obj.isinv(sim(:, 1:n), simi, 'tol', itol) || any(~evaluated, 'all'), "SIMI = SIM(:, 1:N)^{-1}", srname);
            end

        end
        function [nfilt, cfilt, confilt, ffilt, xfilt] = initfilt(~, conmat, ctol, cweight, cval, fval, sim, evaluated, nfilt, cfilt, confilt, ffilt, xfilt)
            %--------------------------------------------------------------------------------------------------%
            % This subroutine initializes the filters (XFILT, etc) that will be used when selecting X at the
            % end of the solver.
            % N.B.:
            % 1. Why not initialize the filters using XHIST, etc? Because the history is empty if the user
            % chooses not to output it.
            % 2. We decouple INITXFC and INITFILT so that it is easier to parallelize the former if needed.
            %--------------------------------------------------------------------------------------------------%

            % Common modules
            consts_obj = consts_mod();
            debug_obj = debug_mod();
            infnan_obj = infnan_mod();
            selectx_obj = selectx_mod();

            % Inputs








            % In-outputs






            % Local variables
            srname = "INITFILT";
            i = NaN;
            m = NaN;
            maxfilt = NaN;
            n = NaN;
            x = NaN(size(sim, 1), 1);

            % Sizes
            m = size(conmat, 1);
            n = size(sim, 1);
            maxfilt = fix(numel(ffilt));

            % Preconditions
            if consts_obj.DEBUGGING
                debug_obj.assert(m >= 0, "M >= 0", srname);
                debug_obj.assert(n >= 1, "N >= 1", srname);
                debug_obj.assert(maxfilt >= 1, "MAXFILT >= 1", srname);
                debug_obj.assert(size(confilt, 1) == m && size(confilt, 2) == maxfilt, "SIZE(CONFILT) == [M, MAXFILT]", srname);
                debug_obj.assert(numel(cfilt) == maxfilt, "SIZE(CFILT) == MAXFILT", srname);
                debug_obj.assert(size(xfilt, 1) == n && size(xfilt, 2) == maxfilt, "SIZE(XFILT) == [N, MAXFILT]", srname);
                debug_obj.assert(numel(ffilt) == maxfilt, "SIZE(FFILT) == MAXFILT", srname);
                debug_obj.assert(size(conmat, 1) == m && size(conmat, 2) == n + 1, "SIZE(CONMAT) = [M, N+1]", srname);
                debug_obj.assert(~any(infnan_obj.is_nan(conmat) | infnan_obj.is_posinf(conmat), 'all'), "CONMAT does not contain NaN/+Inf", srname);
                debug_obj.assert(numel(cval) == n + 1 && ~any(cval < 0 | infnan_obj.is_nan(cval) | infnan_obj.is_posinf(cval), 'all'), "SIZE(CVAL) == N+1 and CVAL does not contain negative values or NaN/+Inf", srname);
                debug_obj.assert(numel(fval) == n + 1 && ~any(infnan_obj.is_nan(fval) | infnan_obj.is_posinf(fval), 'all'), "SIZE(FVAL) == N+1 and FVAL does not contain NaN/+Inf", srname);
                debug_obj.assert(size(sim, 1) == n && size(sim, 2) == n + 1, "SIZE(SIM) == [N, N+1]", srname);
                debug_obj.assert(all(infnan_obj.is_finite(sim), 'all'), "SIM is finite", srname);
                debug_obj.assert(all(sum(abs(sim(:, 1:n)), 1).' > 0, 'all'), "SIM(:, 1:N) has no zero column", srname);
                debug_obj.assert(numel(evaluated) == n + 1, "SIZE(EVALUATED) == N + 1", srname);
            end

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
            if consts_obj.DEBUGGING
                debug_obj.assert(nfilt <= maxfilt, "NFILT <= MAXFILT", srname);
                debug_obj.assert(size(confilt, 1) == m && size(confilt, 2) == maxfilt, "SIZE(CONFILT) == [M, MAXFILT]", srname);
                debug_obj.assert(~any(infnan_obj.is_nan(confilt(:, 1:nfilt)) | infnan_obj.is_posinf(confilt(:, 1:nfilt)), 'all'), "CONFILT does not contain NaN/+Inf", srname);
                debug_obj.assert(numel(cfilt) == maxfilt, "SIZE(CFILT) == MAXFILT", srname);
                debug_obj.assert(~any(cfilt(1:nfilt) < 0 | infnan_obj.is_nan(cfilt(1:nfilt)) | infnan_obj.is_posinf(cfilt(1:nfilt)), 'all'), "CFILT does not contain negative values or NaN/Inf", srname);
                debug_obj.assert(size(xfilt, 1) == n && size(xfilt, 2) == maxfilt, "SIZE(XFILT) == [N, MAXFILT]", srname);
                debug_obj.assert(~any(infnan_obj.is_nan(xfilt(:, 1:nfilt)), 'all'), "XFILT does not contain NaN", srname);
                % The last calculated X can be Inf (finite + finite can be Inf numerically).
                debug_obj.assert(numel(ffilt) == maxfilt, "SIZE(FFILT) == MAXFILT", srname);
                debug_obj.assert(~any(infnan_obj.is_nan(ffilt(1:nfilt)) | infnan_obj.is_posinf(ffilt(1:nfilt)), 'all'), "FFILT does not contain NaN/+Inf", srname);
            end
        end

    end
end