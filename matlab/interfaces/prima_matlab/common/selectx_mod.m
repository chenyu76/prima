classdef selectx_mod
    %--------------------------------------------------------------------------------------------------%
    % This module provides subroutines that ensure the returned X is optimal among all the calculated
    % points in the sense that no other point achieves both lower function value and lower constraint
    % violation at the same time. The module is needed only in the constrained case.
    %
    % Coded by Zaikun ZHANG (www.zhangzk.net).
    %
    % Started: September 2021
    %
    % Last Modified: Saturday, March 02, 2024 AM12:42:28
    %--------------------------------------------------------------------------------------------------%

    methods
        function varargout = isbetter(obj, varargin)
            if numel(varargin) == 5 && isscalar(varargin{1}) && isscalar(varargin{2}) && isvector(varargin{3}) && isvector(varargin{4})
                [varargout{1:nargout}] = obj.isbetter01(varargin{:});
            elseif numel(varargin) == 5 && isvector(varargin{1}) && isvector(varargin{2}) && isscalar(varargin{3}) && isscalar(varargin{4})
                [varargout{1:nargout}] = obj.isbetter10(varargin{:});
            end
        end
        function [nfilt, cfilt, ffilt, xfilt, confilt] = savefilt(obj, cstrv, ctol, cweight, f, x, nfilt, cfilt, ffilt, xfilt, varargin)
            %--------------------------------------------------------------------------------------------------%
            % This subroutine saves X, F, and CSTRV in XFILT, FFILT, and CFILT (and CONSTR in CONFILT if they
            % are present), unless a vector in XFILT(:, 1:NFILT) is better than X. If X is better than some
            % vectors in XFILT(:, 1:NFILT), then these vectors will be removed. If X is not better than any of X
            % FILT(:, 1:NFILT) but NFILT=MAXFILT, then we remove a column from XFILT according to the merit
            % function PHI = FFILT + CWEIGHT * MAX(CFILT - CTOL, ZERO).
            % N.B.:
            % 1. Only XFILT(:, 1:NFILT) and FFILT(:, 1:NFILT) etc contains valid information, while
            % XFILT(:, NFILT+1:MAXFILT) and FFILT(:, NFILT+1:MAXFILT) etc are not initialized yet.
            % 2. We decide whether an X is better than another by the ISBETTER function.
            %--------------------------------------------------------------------------------------------------%
            consts_obj = consts_mod();
            debug_obj = debug_mod();
            infnan_obj = infnan_mod();
            linalg_obj = linalg_mod();


            % Inputs




            % N
            % M

            % In-outputs

            % MAXFILT
            % MAXFILT
            % (N, MAXFILT)
            % (M, MAXFILT)

            % Local variables
            srname = "SAVEFILT";
            index_to_keep = NaN(numel(ffilt), 1);
            kworst = NaN;
            m = NaN;
            maxfilt = NaN;
            n = NaN;
            keep = false(nfilt, 1);
            cfilt_shifted = NaN(numel(ffilt), 1);
            cref = NaN;
            fref = NaN;
            phi = NaN(numel(ffilt), 1);
            phimax = NaN;

            % Sizes
            ipObj = inputParser();
            addParameter(ipObj, 'constr', NaN);
            addParameter(ipObj, 'confilt', NaN);
            parse(ipObj, varargin{:});
            constr = ipObj.Results.constr;
            confilt = ipObj.Results.confilt;
            if ismember('constr', ipObj.UsingDefaults)
                m = 0;
            else
                m = fix(numel(constr));
            end
            n = fix(numel(x));
            maxfilt = fix(numel(ffilt));

            % Preconditions
            if consts_obj.DEBUGGING
                % Check the size of X.
                debug_obj.assert(n >= 1, "N >= 1", srname);
                % Check CWEIGHT and CTOL
                debug_obj.assert(cweight >= 0, "CWEIGHT >= 0", srname);
                debug_obj.assert(ctol >= 0, "CTOL >= 0", srname);
                % Check NFILT
                debug_obj.assert(nfilt >= 0 && nfilt <= maxfilt, "0 <= NFILT <= MAXFILT", srname);
                % Check the sizes of XFILT, FFILT, CFILT.
                debug_obj.assert(maxfilt >= 1, "MAXFILT >= 1", srname);
                debug_obj.assert(size(xfilt, 1) == n && size(xfilt, 2) == maxfilt, "SIZE(XFILT) == [N, MAXFILT]", srname);
                debug_obj.assert(numel(cfilt) == maxfilt, "SIZE(CFILT) == MAXFILT", srname);
                % Check the values of XFILT, FFILT, CFILT.
                debug_obj.assert(~any(infnan_obj.is_nan(xfilt(:, 1:nfilt)), 'all'), "XFILT does not contain NaN", srname);
                debug_obj.assert(~any(infnan_obj.is_nan(ffilt(1:nfilt)) | infnan_obj.is_posinf(ffilt(1:nfilt)), 'all'), "FFILT does not contain NaN/+Inf", srname);
                debug_obj.assert(~any(cfilt(1:nfilt) < 0 | infnan_obj.is_nan(cfilt(1:nfilt)) | infnan_obj.is_posinf(cfilt(1:nfilt)), 'all'), "CFILT does not contain nonnegative values of NaN/+Inf", srname);
                % Check the values of X, F, CSTRV.
                % X does not contain NaN if X0 does not and the trust-region/geometry steps are proper.
                debug_obj.assert(~any(infnan_obj.is_nan(x), 'all'), "X does not contain NaN", srname);
                % F cannot be NaN/+Inf due to the moderated extreme barrier.
                debug_obj.assert(~(infnan_obj.is_nan_sp(f) || infnan_obj.is_posinf(f)), "F is not NaN/+Inf", srname);
                % CSTRV cannot be NaN/+Inf due to the moderated extreme barrier.
                debug_obj.assert(~(cstrv < 0 || infnan_obj.is_nan_sp(cstrv) || infnan_obj.is_posinf(cstrv)), "CSTRV is nonnegative and not NaN/+Inf", srname);
                % Check CONSTR and CONFILT.
                debug_obj.assert(~ismember('constr', ipObj.UsingDefaults) == ((~ismember('confilt', ipObj.UsingDefaults)) || (nargout >= 5)), "CONSTR and CONFILT are both present or both absent", srname);
                if ~ismember('constr', ipObj.UsingDefaults)
                    % CONSTR cannot contain NaN/+Inf due to the moderated extreme barrier.
                    debug_obj.assert(~any(infnan_obj.is_nan(constr) | infnan_obj.is_posinf(constr), 'all'), "CONSTR does not contain NaN/+Inf", srname);
                    debug_obj.assert(size(confilt, 1) == m && size(confilt, 2) == maxfilt, "SIZE(CONFILT) == [M, MAXFILT]", srname);
                    debug_obj.assert(~any(infnan_obj.is_nan(confilt(:, 1:nfilt)) | infnan_obj.is_posinf(confilt(:, 1:nfilt)), 'all'), "CONFILT does not contain NaN/+Inf", srname);
                end
            end

            %====================%
            % Calculation starts %
            %====================%

            % Return immediately if any column of XFILT is better than X. Note that ISBETTER checks "strictly
            % better", handling NaN/Inf properly, but we need "non-strictly better" here, allowing equality.
            % This is why we need to supplement ISBETTER with (FFILT <= F .AND. CFILT <= CSTRV).
            if any(obj.isbetter10(ffilt(1:nfilt), cfilt(1:nfilt), f, cstrv, ctol), 'all') || any(ffilt(1:nfilt) <= f & cfilt(1:nfilt) <= cstrv, 'all')
                return
            end

            % Decide which columns of XFILT to keep.
            keep(:) = (~obj.isbetter01(f, cstrv, ffilt(1:nfilt), cfilt(1:nfilt), ctol));

            % If NFILT == MAXFILT and X is not better than any column of XFILT, then we remove the worst column
            % of XFILT according to the merit function PHI = FFILT + CWEIGHT * MAX(CFILT - CTOL, ZERO).
            if nnz(keep) == maxfilt
                % In this case, NFILT = SIZE(KEEP) = COUNT(KEEP) = MAXFILT > 0.
                cfilt_shifted(:) = max(cfilt - ctol, consts_obj.ZERO);
                if cweight <= 0
                    phi(:) = ffilt;
                elseif infnan_obj.is_posinf(cweight)
                    phi(:) = cfilt_shifted;
                    % We should not use CFILT here; if MAX(CFILT_SHIFTED) is attained at multiple indices, then
                    % we will check FFILT to exhaust the remaining degree of freedom.

                else
                    phi(:) = max(ffilt, -consts_obj.REALMAX) + cweight * cfilt_shifted;
                    % MAX(FFILT, -REALMAX) makes sure that PHI will not contain NaN (unless there is a bug).
                end
                % We select X to maximize PHI. In case there are multiple maximizers, we take the one with the
                % largest CSTRV_SHIFTED; if there are more than one choices, we take the one with the largest F;
                % if there are several candidates, we take the one with the largest CSTRV; if the last comparison
                % still leads to more than one possibilities, then they are equally bad and we choose the first.
                % N.B.:
                % 1. This process is the opposite of selecting KOPT in SELECTX.
                % 2. In finite-precision arithmetic, PHI_1 == PHI_2 and CSTRV_SHIFTED_1 == CSTRV_SHIFTED_2 do
                % not ensure that F_1 == F_2!
                phimax = max(phi, [], 'all');
                cref = max(fortran.merge('tsource', cfilt_shifted, 'fsource', -realmax, 'mask', (phi >= phimax)), [], 'all');
                fref = max(fortran.merge('tsource', ffilt, 'fsource', -realmax, 'mask', (cfilt_shifted >= cref)), [], 'all');
                kworst = fix(fortran.maxloc(cfilt, 'mask', (ffilt >= fref), 'dim', 1));
                %%MATLAB: cmax = max(cfilt(ffilt >= fref)); kworst = find(ffilt >= fref & ~(cfilt < cmax), 1,'first');
                if kworst < 1 || kworst > numel(keep)
                    % For security. Should not happen.
                    kworst = 1;
                end
                keep(kworst) = false;
            end

            nfilt = fix(nnz(keep));
            index_to_keep(1:nfilt) = linalg_obj.trueloc(keep);
            xfilt(:, 1:nfilt) = xfilt(:, index_to_keep(1:nfilt));
            ffilt(1:nfilt) = ffilt(index_to_keep(1:nfilt));
            cfilt(1:nfilt) = cfilt(index_to_keep(1:nfilt));
            if ((~ismember('confilt', ipObj.UsingDefaults)) || (nargout >= 5)) && ~ismember('constr', ipObj.UsingDefaults)
                confilt(:, 1:nfilt) = confilt(:, index_to_keep(1:nfilt));
            end

            nfilt = nfilt + 1;
            xfilt(:, nfilt) = x;
            ffilt(nfilt) = f;
            cfilt(nfilt) = cstrv;
            if ((~ismember('confilt', ipObj.UsingDefaults)) || (nargout >= 5)) && ~ismember('constr', ipObj.UsingDefaults)
                confilt(:, nfilt) = constr;
            end

            %====================%
            %  Calculation ends  %
            %====================%

            % Postconditions
            if consts_obj.DEBUGGING
                % Check NFILT and the sizes of XFILT, FFILT, CFILT.
                debug_obj.assert(nfilt >= 1 && nfilt <= maxfilt, "1 <= NFILT <= MAXFILT", srname);
                debug_obj.assert(size(xfilt, 1) == n && size(xfilt, 2) == maxfilt, "SIZE(XFILT) == [N, MAXFILT]", srname);
                debug_obj.assert(numel(ffilt) == maxfilt, "SIZE(FFILT) = MAXFILT", srname);
                debug_obj.assert(numel(cfilt) == maxfilt, "SIZE(CFILT) = MAXFILT", srname);
                % Check the values of XFILT, FFILT, CFILT.
                debug_obj.assert(~any(infnan_obj.is_nan(xfilt(:, 1:nfilt)), 'all'), "XFILT does not contain NaN", srname);
                debug_obj.assert(~any(infnan_obj.is_nan(ffilt(1:nfilt)) | infnan_obj.is_posinf(ffilt(1:nfilt)), 'all'), "FFILT does not contain NaN/+Inf", srname);
                debug_obj.assert(~any(cfilt(1:nfilt) < 0 | infnan_obj.is_nan(cfilt(1:nfilt)) | infnan_obj.is_posinf(cfilt(1:nfilt)), 'all'), "CFILT does not contain nonnegative values of NaN/+Inf", srname);
                % Check that no point in the filter is better than X, and X is better than no point.
                debug_obj.assert(~any(obj.isbetter10(ffilt(1:nfilt), cfilt(1:nfilt), f, cstrv, ctol), 'all'), "No point in the filter is better than X", srname);
                debug_obj.assert(~any(obj.isbetter01(f, cstrv, ffilt(1:nfilt), cfilt(1:nfilt), ctol), 'all'), "X is better than no point in the filter", srname);
                % Check CONFILT.
                if (~ismember('confilt', ipObj.UsingDefaults)) || (nargout >= 5)
                    debug_obj.assert(size(confilt, 1) == m && size(confilt, 2) == maxfilt, "SIZE(CONFILT) == [M, MAXFILT]", srname);
                    debug_obj.assert(~any(infnan_obj.is_nan(confilt(:, 1:nfilt)) | infnan_obj.is_posinf(confilt(:, 1:nfilt)), 'all'), "CONFILT does not contain NaN/+Inf", srname);
                end
            end

        end
        function kopt = selectx(obj, fhist, chist, cweight, ctol)
            %--------------------------------------------------------------------------------------------------%
            % This subroutine selects X according to the FHIST and CHIST, which represents (a part of) history
            % of F and CSTRV. Normally, FHIST and CHIST are not the full history but only a filter, e.g., FFILT
            % and CFILT generated by SAVEFILT. However, we name them as FHIST and CHIST because the [F, CSTRV]
            % in a filter should not dominate each other, but this subroutine does NOT assume such a property.
            % N.B.: CTOL is the tolerance of constraint violation (CSTRV). A point is considered feasible if
            % its constraint violation is at most CTOL. Note that CTOL is absolute, not relative.
            %--------------------------------------------------------------------------------------------------%

            consts_obj = consts_mod();
            infnan_obj = infnan_mod();
            debug_obj = debug_mod();


            % Inputs





            % Outputs
            kopt = NaN;

            % Local variables
            srname = "SELECTX";
            nhist = NaN;
            chist_shifted = NaN(numel(fhist), 1);
            cmin = NaN;
            cref = NaN;
            fref = NaN;
            phi = NaN(numel(fhist), 1);
            phimin = NaN;

            % Sizes
            nhist = fix(numel(fhist));

            % Preconditions
            if consts_obj.DEBUGGING
                debug_obj.assert(nhist >= 1, "SIZE(FHIST) >= 1", srname);
                debug_obj.assert(numel(chist) == nhist, "SIZE(FHIST) == SIZE(CHIST)", srname);
                debug_obj.assert(~any(infnan_obj.is_nan(fhist) | infnan_obj.is_posinf(fhist), 'all'), "FHIST does not contain NaN/+Inf", srname);
                debug_obj.assert(~any(chist < 0 | infnan_obj.is_nan(chist) | infnan_obj.is_posinf(chist), 'all'), "CHIST does not contain nonnegative values or NaN/+Inf", srname);
                debug_obj.assert(cweight >= 0, "CWEIGHT >= 0", srname);
                debug_obj.assert(ctol >= 0, "CTOL >= 0", srname);
            end

            %====================%
            % Calculation starts %
            %====================%

            % We select X among the points with F < FREF and CSTRV < CREF.
            % Do NOT use F <= FREF, because F == FREF (FUNCMAX or REALMAX) may mean F == INF in practice!
            if any(fhist < consts_obj.FUNCMAX & chist < consts_obj.CONSTRMAX, 'all')
                fref = consts_obj.FUNCMAX;
                cref = consts_obj.CONSTRMAX;
            elseif any(fhist < consts_obj.REALMAX & chist < consts_obj.CONSTRMAX, 'all')
                fref = consts_obj.REALMAX;
                cref = consts_obj.CONSTRMAX;
            elseif any(fhist < consts_obj.FUNCMAX & chist < consts_obj.REALMAX, 'all')
                fref = consts_obj.FUNCMAX;
                cref = consts_obj.REALMAX;
            else
                fref = consts_obj.REALMAX;
                cref = consts_obj.REALMAX;
            end

            if any(fhist < fref & chist < cref, 'all')
                % Shift the constraint violations by CTOL, so that CSTRV <= CTOL is regarded as no violation.
                chist_shifted(:) = max(chist - ctol, consts_obj.ZERO);
                % CMIN is the minimal shifted constraint violation attained in the history.
                cmin = min(fortran.merge('tsource', chist_shifted, 'fsource', realmax, 'mask', (fhist < fref)), [], 'all');
                % We consider only the points whose shifted constraint violations are at most the CREF below.
                % N.B.: Without taking MAX(EPS, .), CREF would be 0 if CMIN = 0. In that case, asking for
                % CSTRV_SHIFTED < CREF would be WRONG!
                cref = max(consts_obj.EPS, consts_obj.TWO * cmin);
                % We use the following PHI as our merit function to select X.
                if cweight <= 0
                    phi(:) = fhist;
                elseif infnan_obj.is_posinf(cweight)
                    phi(:) = chist_shifted;
                    % We should not use CHIST here; if MIN(CHIST_SHIFTED) is attained at multiple indices, then
                    % we will check FHIST to exhaust the remaining degree of freedom.

                else
                    phi(:) = max(fhist, -consts_obj.REALMAX) + cweight * chist_shifted;
                    % MAX(FHIST, -REALMAX) makes sure that PHI will not contain NaN (unless there is a bug).
                end
                % We select X to minimize PHI subject to F < FREF and CSTRV_SHIFTED <= CREF (see the comments
                % above for the reason of taking "<" and "<=" in these two constraints). In case there are
                % multiple minimizers, we take the one with the least CSTRV_SHIFTED; if there are more than one
                % choices, we take the one with the least F; if there are several candidates, we take the one
                % with the least CSTRV; if the last comparison still leads to more than one possibilities, then
                % they are equally good and we choose the first.
                % N.B.:
                % 1. This process is the opposite of selecting KWORST in SAVEFILT.
                % 2. In finite-precision arithmetic, PHI_1 == PHI_2 and CSTRV_SHIFTED_1 == CSTRV_SHIFTED_2 do
                % not ensure that F_1 == F_2!
                phimin = min(fortran.merge('tsource', phi, 'fsource', realmax, 'mask', (fhist < fref & chist_shifted <= cref)), [], 'all');
                cref = min(fortran.merge('tsource', chist_shifted, 'fsource', realmax, 'mask', (fhist < fref & phi <= phimin)), [], 'all');
                fref = min(fortran.merge('tsource', fhist, 'fsource', realmax, 'mask', (chist_shifted <= cref)), [], 'all');
                kopt = fix(fortran.minloc(chist, 'mask', (fhist <= fref), 'dim', 1));
                %%MATLAB: cmin = min(chist(fhist <= fref)); kopt = find(fhist <= fref & ~(chist > cmin), 1,'first');
            else
                kopt = nhist;
            end

            %====================%
            %  Calculation ends  %
            %====================%

            % Postconditions
            if consts_obj.DEBUGGING
                debug_obj.assert(kopt >= 1 && kopt <= nhist, "1 <= KOPT <= SIZE(FHIST)", srname);
                debug_obj.assert(~any(obj.isbetter10(fhist(1:nhist), chist(1:nhist), fhist(kopt), chist(kopt), ctol), 'all'), "No point in the history is better than X", srname);
            end

        end
        function is_better = isbetter00(~, f1, c1, f2, c2, ctol)
            %--------------------------------------------------------------------------------------------------%
            % This function compares whether FC1 = (F1, C1) is (strictly) better than FC2 = (F2, C2), which
            % basically means that (F1 < F2 and C1 <= C2) or (F1 <= F2 and C1 < C2).
            % It takes care of the cases where some of these values are NaN or Inf, even though some cases
            % should never happen due to the moderated extreme barrier.
            % At return, BETTER = TRUE if and only if (F1, C1) is better than (F2, C2).
            % Here, C means constraint violation, which is a nonnegative number.
            %--------------------------------------------------------------------------------------------------%

            consts_obj = consts_mod();
            infnan_obj = infnan_mod();
            debug_obj = debug_mod();

            % Inputs






            % Outputs
            is_better = false;

            % Local variables
            srname = "ISBETTER";
            cref = NaN;

            % Preconditions
            if consts_obj.DEBUGGING
                debug_obj.assert(~any(infnan_obj.is_nan([f1, c1]) | infnan_obj.is_posinf([f2, c2]), 'all'), "FC1 does not contain NaN/+Inf", srname);
                debug_obj.assert(~any(infnan_obj.is_nan([f2, c2]) | infnan_obj.is_posinf([f2, c2]), 'all'), "FC2 does not contain NaN/+Inf", srname);
                debug_obj.assert(c1 >= 0 && c2 >= 0, "C1 >= 0, C2 >= 0", srname);
                debug_obj.assert(ctol >= 0, "CTOL >= 0", srname);
            end

            %====================%
            % Calculation starts %
            %====================%

            is_better = false;
            % Even though NaN/+Inf should not occur in FC1 or FC2 due to the moderated extreme barrier, for
            % security and robustness, the code below does not make this assumption.
            is_better = is_better || (any(infnan_obj.is_nan([f2, c2]) | infnan_obj.is_posinf([f2, c2]), 'all') && ~any(infnan_obj.is_nan([f1, c1]) | infnan_obj.is_posinf([f1, c1]), 'all'));
            is_better = is_better || (f1 < f2 && c1 <= c2);
            is_better = is_better || (f1 <= f2 && c1 < c2);
            % If C1 <= CTOL and C2 is significantly larger/worse than CTOL, i.e., C2 > MAX(CTOL, CREF),
            % then FC1 is better than FC2 as long as F1 < REALMAX. Normally CREF >= CTOL so MAX(CTOL, CREF)
            % is indeed CREF. However, this may not be true if CTOL > 1E-1*CONSTRMAX.
            cref = consts_obj.TEN * max(consts_obj.EPS, min(ctol, 1.0e-2 * consts_obj.CONSTRMAX)); % The MIN avoids overflow.
            is_better = is_better || (f1 < consts_obj.REALMAX && c1 <= ctol && (c2 > max(ctol, cref) || infnan_obj.is_nan_sp(c2)));

            %====================%
            %  Calculation ends  %
            %====================%

            % Postconditions
            if consts_obj.DEBUGGING
                % Even though NaN/+Inf should not occur in FC1 due to moderated extreme barrier, for security
                % and robustness, the code below does not make this assumption.
                debug_obj.assert(~(is_better && any(infnan_obj.is_nan([f1, c1]) | infnan_obj.is_posinf([f1, c1]), 'all')), "IS_BETTER cannot be true if [F1, C1] contains NaN/+Inf", srname);
                debug_obj.assert(is_better || any(infnan_obj.is_nan([f1, c1]) | infnan_obj.is_posinf([f1, c1]), 'all') || ~any(infnan_obj.is_nan([f2, c2]) | infnan_obj.is_posinf([f2, c2]), 'all'), "if [F2, C2] contains NaN/+Inf, then either IS_BETTER is true or [F1, C1] contains NaN/+Inf", srname);
                debug_obj.assert(~(is_better && f1 >= f2 && c1 >= c2), "[F1, C1] >= [F2, C2] and IS_BETTER cannot be both true", srname);
                debug_obj.assert(is_better || ~(f1 <= f2 && c1 < c2), "if [F1, C1] <= [F2, C2] but not equal, then IS_BETTER must be true", srname);
                debug_obj.assert(is_better || ~(f1 < f2 && c1 <= c2), "if [F1, C1] <= [F2, C2] but not equal, then IS_BETTER must be true", srname);
            end

        end
        function is_better = isbetter10(obj, f1, c1, f2, c2, ctol)
            consts_obj = consts_mod();
            debug_obj = debug_mod();
            infnan_obj = infnan_mod();
            memory_obj = memory_mod();

            % Inputs






            % Outputs
            is_better = false(1, 1);

            % Local variables
            srname = "ISBETTER10";
            i = NaN;
            nfc = NaN;

            % Sizes
            nfc = fix(numel(f1));

            % Preconditions
            if consts_obj.DEBUGGING
                debug_obj.assert(nfc >= 0, "NFC >= 0", srname);
                debug_obj.assert(numel(f1) == numel(c1), "SIZE(F1) == SIZE(C1)", srname);
                debug_obj.assert(~any(infnan_obj.is_nan(f1) | infnan_obj.is_posinf(f1), 'all'), "F1 does not contain NaN/+Inf", srname);
                debug_obj.assert(~any(infnan_obj.is_nan(c1) | infnan_obj.is_posinf(c1), 'all'), "C1 does not contain NaN/+Inf", srname);
                debug_obj.assert(~(infnan_obj.is_nan_sp(f2) || infnan_obj.is_posinf(f2)), "F2 is not NaN/+Inf", srname);
                debug_obj.assert(~(infnan_obj.is_nan_sp(c2) || infnan_obj.is_posinf(c2)), "C2 is not NaN/+Inf", srname);
                debug_obj.assert(all(c1 >= 0, 'all') && c2 >= 0, "C1 >= 0, C2 >= 0", srname);
                debug_obj.assert(ctol >= 0, "CTOL >= 0", srname);
            end

            %====================%
            % Calculation starts %
            %====================%

            is_better = memory_obj.alloc_lvector(is_better, nfc);
            is_better = reshape(cell2mat(arrayfun(@(i) obj.isbetter00(f1(i), c1(i), f2, c2, ctol), (1:nfc), "UniformOutput", false)), [], 1);

            %====================%
            %  Calculation ends  %
            %====================%

            % Postconditions
            if consts_obj.DEBUGGING
                debug_obj.assert(numel(is_better) == numel(f1), "SIZE(IS_BETTER) == SIZE(F1)", srname);
            end

        end
        function is_better = isbetter01(obj, f1, c1, f2, c2, ctol)
            consts_obj = consts_mod();
            debug_obj = debug_mod();
            infnan_obj = infnan_mod();
            memory_obj = memory_mod();

            % Inputs






            % Outputs
            is_better = false(1, 1);

            % Local variables
            srname = "ISBETTER01";
            i = NaN;
            nfc = NaN;

            % Sizes
            nfc = fix(numel(f2));

            % Preconditions
            if consts_obj.DEBUGGING
                debug_obj.assert(nfc >= 0, "NFC >= 0", srname);
                debug_obj.assert(~(infnan_obj.is_nan_sp(f1) || infnan_obj.is_posinf(f1)), "F1 is not NaN/+Inf", srname);
                debug_obj.assert(~(infnan_obj.is_nan_sp(c1) || infnan_obj.is_posinf(c1)), "C1 is not NaN/+Inf", srname);
                debug_obj.assert(numel(f2) == numel(c2), "SIZE(F2) == SIZE(C2)", srname);
                debug_obj.assert(~any(infnan_obj.is_nan(f2) | infnan_obj.is_posinf(f2), 'all'), "F2 does not contain NaN/+Inf", srname);
                debug_obj.assert(~any(infnan_obj.is_nan(c2) | infnan_obj.is_posinf(c2), 'all'), "C2 does not contain NaN/+Inf", srname);
                debug_obj.assert(c1 >= 0 && all(c2 >= 0, 'all'), "C1 >= 0, C2 >= 0", srname);
                debug_obj.assert(ctol >= 0, "CTOL >= 0", srname);
            end

            %====================%
            % Calculation starts %
            %====================%

            is_better = memory_obj.alloc_lvector(is_better, nfc);
            is_better = reshape(cell2mat(arrayfun(@(i) obj.isbetter00(f1, c1, f2(i), c2(i), ctol), (1:nfc), "UniformOutput", false)), [], 1);

            %====================%
            %  Calculation ends  %
            %====================%

            % Postconditions
            if consts_obj.DEBUGGING
                debug_obj.assert(numel(is_better) == numel(f2), "SIZE(IS_BETTER) == SIZE(F2)", srname);
            end

        end

    end
end