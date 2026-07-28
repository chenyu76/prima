classdef history_mod
    %--------------------------------------------------------------------------------------------------%
    % This module provides subroutines that handle the X/F/C histories of the solver, taking into
    % account that MAXHIST may be smaller than NF.
    %
    % Coded by Zaikun ZHANG (www.zhangzk.net).
    %
    % Started: July 2020
    %
    % Last Modified: Thursday, April 04, 2024 PM10:15:40
    %--------------------------------------------------------------------------------------------------%

    methods
        function [maxhist, xhist, fhist, chist, conhist] = prehist(~, maxhist, n, output_xhist, output_fhist, varargin)
            %--------------------------------------------------------------------------------------------------%
            % This subroutine revises MAXHIST according to MAXHISTMEM, and allocates memory for the history.
            % In MATLAB/Python/Julia/R implementation, we should simply set MAXHIST = MAXFUN and initialize
            % XHIST = NaN(N, MAXFUN), FHIST = NaN(1, MAXFUN), CHIST = NaN(1, MAXFUN), CONHIST = NaN(M, MAXFUN),
            % if they are requested; replace MAXFUN with 0 for the history that is not requested.
            %--------------------------------------------------------------------------------------------------%
            consts_obj = prima_mat.common.consts_mod();
            debug_obj = prima_mat.common.debug_mod();
            linalg_obj = prima_mat.common.linalg_mod();
            memory_obj = prima_mat.common.memory_mod();

            % Inputs



            % In-outputs


            % Outputs



            % Local variables
            srname = "PREHIST";
            % INTEGER(IK) may overflow if IK corresponds to the 16-bit integer.


            % Preconditions
            ipObj = inputParser();
            addParameter(ipObj, 'output_chist', false);
            addParameter(ipObj, 'chist', NaN);
            addParameter(ipObj, 'm', NaN);
            addParameter(ipObj, 'output_conhist', false);
            addParameter(ipObj, 'conhist', NaN);
            parse(ipObj, varargin{:});
            output_chist = ipObj.Results.output_chist;
            chist = ipObj.Results.chist;
            m = ipObj.Results.m;
            output_conhist = ipObj.Results.output_conhist;
            conhist = ipObj.Results.conhist;
            if consts_obj.DEBUGGING
                debug_obj.assert(maxhist >= 0, "MAXHIST >= 0", srname);
                debug_obj.assert(n >= 1, "N >= 1", srname);
                if ~ismember('m', ipObj.UsingDefaults)
                    debug_obj.assert(m >= 0, "M >= 0", srname);
                end
                debug_obj.assert(~ismember('output_chist', ipObj.UsingDefaults) == (nargout >= 4), "OUTPUT_CHIST and CHIST are both present or both absent", srname);
                debug_obj.assert((~ismember('m', ipObj.UsingDefaults) == (nargout >= 5)) && (~ismember('output_conhist', ipObj.UsingDefaults) == (nargout >= 5)), "M, OUTPUT_CONHIST, and CONHIST are all present or all absent", srname);
            end

            %====================%
            % Calculation starts %
            %====================%

            % Save the input value of MAXHIST for debugging.
            maxhist_in = maxhist;

            % Revise MAXHIST according to MAXHISTMEM, i.e., the maximal memory allowed for the history.
            % N.B.: The `UNIT_MEMO = INT(*)` below converts integers to the default integer kind, which is the
            % kind of UNIT_MEMO. Fortran compilers may complain without the conversion. It is not needed in
            % Python/MATLAB/Julia/R. Meanwhile, INT(OUTPUT_*HIST) converts booleans to integers.
            unit_memo = linalg_obj.logical_to_int(output_xhist) * n + linalg_obj.logical_to_int(output_fhist);
            if ~ismember('output_chist', ipObj.UsingDefaults) && nargout >= 4
                unit_memo = unit_memo + linalg_obj.logical_to_int(output_chist);
            end
            if ~ismember('m', ipObj.UsingDefaults) && ~ismember('output_conhist', ipObj.UsingDefaults) && nargout >= 5
                unit_memo = unit_memo + linalg_obj.logical_to_int(output_conhist) * m;
            end
            unit_memo = unit_memo * memory_obj.size_of_sp(0.0); % INT(*) avoids overflow when IK is 16-bit.
            if unit_memo <= 0
                % No output of history is requested
                maxhist = 0;
            elseif maxhist > consts_obj.MAXHISTMEM / unit_memo
                maxhist = fix(consts_obj.MAXHISTMEM / unit_memo); % Integer division.
                % We cannot simply set MAXHIST = MIN(MAXHIST, MAXHISTMEM/UNIT_MEMO), as they may not have
                % the same kind, and compilers may complain. We may convert them, but overflow may occur.

            end

            xhist = memory_obj.alloc_rmatrix_sp(n, maxhist * linalg_obj.logical_to_int(output_xhist));
            fhist = memory_obj.alloc_rvector_sp(maxhist * linalg_obj.logical_to_int(output_fhist));
            % Even if OUTPUT_CHIST is FALSE, CHIST still needs to be allocated.
            if ~ismember('output_chist', ipObj.UsingDefaults) && nargout >= 4
                chist = memory_obj.alloc_rvector_sp(maxhist * linalg_obj.logical_to_int(output_chist));
            end
            % Even if OUTPUT_CONHIST is FALSE, CONHIST still needs to be allocated.
            if ~ismember('m', ipObj.UsingDefaults) && ~ismember('output_conhist', ipObj.UsingDefaults) && nargout >= 5
                conhist = memory_obj.alloc_rmatrix_sp(m, maxhist * linalg_obj.logical_to_int(output_conhist));
            end

            %====================%
            %  Calculation ends  %
            %====================%

            % Postconditions
            if consts_obj.DEBUGGING
                debug_obj.assert(maxhist >= 0 && maxhist <= maxhist_in, "0 <= MAXHIST <= MAXHIST_IN", srname);
                debug_obj.assert(maxhist * unit_memo <= consts_obj.MAXHISTMEM, "The history will not take more memory than MAXHISTMEM", srname);
                debug_obj.assert(exist('xhist', 'var'), "XHIST is allocated", srname);
                debug_obj.assert(size(xhist, 1) == n && size(xhist, 2) == maxhist * linalg_obj.logical_to_int(output_xhist), "if XHIST is requested, then SIZE(XHIST) == [N, MAXHIST]; otherwise, SIZE(XHIST) == [N, 0]", srname);
                debug_obj.assert(exist('fhist', 'var'), "FHIST is allocated", srname);
                debug_obj.assert(numel(fhist) == maxhist * linalg_obj.logical_to_int(output_fhist), "if FHIST is requested, then SIZE(FHIST) == MAXHIST; otherwise, SIZE(FHIST) == 0", srname);
                if ~ismember('output_chist', ipObj.UsingDefaults) && nargout >= 4
                    debug_obj.assert(exist('chist', 'var'), "CHIST is allocated", srname);
                    debug_obj.assert(numel(chist) == maxhist * linalg_obj.logical_to_int(output_chist), "if CHIST is requested, then SIZE(CHIST) == MAXHIST; otherwise, SIZE(CHIST) == 0", srname);
                end
                if ~ismember('m', ipObj.UsingDefaults) && ~ismember('output_conhist', ipObj.UsingDefaults) && nargout >= 5
                    debug_obj.assert(exist('conhist', 'var'), "CONHIST is allocated", srname);
                    debug_obj.assert(size(conhist, 1) == m && size(conhist, 2) == maxhist * linalg_obj.logical_to_int(output_conhist), "if CONHIST is requested, then SIZE(CONHIST) == [M, MAXHIST]; otherwise, SIZE(CONHIST) == [M, 0]", srname);
                end
            end
        end
        function [xhist, fhist, chist, conhist] = savehist(~, nf, x, xhist, f, fhist, varargin)
            %--------------------------------------------------------------------------------------------------%
            % This subroutine saves X, F, CSTRV, and CONSTR into XHIST, FHIST, CHIST, and CONHIST respectively.
            %--------------------------------------------------------------------------------------------------%
            consts_obj = prima_mat.common.consts_mod();
            debug_obj = prima_mat.common.debug_mod();
            infnan_obj = prima_mat.common.infnan_mod();
            string_obj = prima_mat.common.string_mod();

            % Inputs



            % In-outputs



            % Local variables



            srname = "SAVEHIST";

            % Sizes
            maxxhist = size(xhist, 2);
            maxfhist = numel(fhist);
            ipObj = inputParser();
            addParameter(ipObj, 'cstrv', NaN);
            addParameter(ipObj, 'chist', NaN);
            addParameter(ipObj, 'constr', NaN);
            addParameter(ipObj, 'conhist', NaN);
            parse(ipObj, varargin{:});
            cstrv = ipObj.Results.cstrv;
            chist = ipObj.Results.chist;
            constr = ipObj.Results.constr;
            conhist = ipObj.Results.conhist;
            if (~ismember('chist', ipObj.UsingDefaults) || nargout >= 3) && ~ismember('cstrv', ipObj.UsingDefaults)
                maxchist = numel(chist);
            else
                maxchist = 0;
            end
            if (~ismember('conhist', ipObj.UsingDefaults) || nargout >= 4) && ~ismember('constr', ipObj.UsingDefaults)
                maxconhist = size(conhist, 2);
            else
                maxconhist = 0;
            end
            maxhist = max(maxxhist, max(maxfhist, max(maxchist, maxconhist)));

            % Preconditions
            if consts_obj.DEBUGGING
                % Called after each function evaluation when debugging; can be expensive.
                % Check the presence of CSTRV, CHIST, CONSTR, CONHIST.
                debug_obj.assert(~ismember('cstrv', ipObj.UsingDefaults) == (~ismember('chist', ipObj.UsingDefaults) || nargout >= 3), "CSTRV and CHIST are both present or both absent", srname);
                debug_obj.assert(~ismember('constr', ipObj.UsingDefaults) == (~ismember('conhist', ipObj.UsingDefaults) || nargout >= 4), "CONSTR and CONHIST are both present or both absent", srname);
                % Check the size of X.
                debug_obj.assert(numel(x) >= 1, "SIZE(X) >= 1", srname);
                % Check the sizes of XHIST, FHIST, CONHIST, CHIST.
                debug_obj.assert(size(xhist, 1) == numel(x) && maxxhist * (maxxhist - maxhist) == 0, "SIZE(XHIST, 1) == SIZE(X), SIZE(XHIST, 2) == 0 or MAXHIST", srname);
                debug_obj.assert(maxfhist * (maxfhist - maxhist) == 0, "SIZE(FHIST) == 0 or MAXHIST", srname);
                debug_obj.assert(maxchist * (maxchist - maxhist) == 0, "SIZE(CHIST) == 0 or MAXHIST", srname);
                if ~ismember('constr', ipObj.UsingDefaults) && (~ismember('conhist', ipObj.UsingDefaults) || nargout >= 4)
                    debug_obj.assert(size(conhist, 1) == numel(constr) && maxconhist * (maxconhist - maxhist) == 0, "SIZE(CONHIST, 1) == SIZE(CONSTR), SIZE(CONHIST, 2) == 0 or MAXNHIST", srname);
                end
                % Check the values of XHIST, FHIST, CHIST, CONHIST, up to the (NF - 1)th position.
                % As long as this subroutine is called, XHIST contains only finite values.
                debug_obj.assert(all(infnan_obj.is_finite(xhist(:, 1:min(nf - 1, maxxhist))), 'all'), "XHIST is finite", srname);
                debug_obj.assert(~any(infnan_obj.is_nan_sp(fhist(1:min(nf - 1, maxfhist))) | infnan_obj.is_posinf(fhist(1:min(nf - 1, maxfhist))), 'all'), "FHIST does not contain NaN/+Inf", srname);
                if ~ismember('chist', ipObj.UsingDefaults) || nargout >= 3
                    debug_obj.assert(~any(chist(1:min(nf - 1, maxchist)) < 0, 'all'), "CHIST does not contain negative values", srname);
                    %------------------------------------------------------------------------------------------%
                    % The following test is not applicable to LINCOA.
                    % %call assert(.not. any(is_nan(chist(1:min(nf - 1_IK, maxchist))) .or. &
                    % %    & is_posinf(chist(1:min(nf - 1_IK, maxchist)))), 'CHIST does not contain NaN/+Inf', srname)
                    %------------------------------------------------------------------------------------------%

                end
                if ~ismember('conhist', ipObj.UsingDefaults) || nargout >= 4
                    debug_obj.assert(~any(infnan_obj.is_nan_sp(conhist(:, 1:min(nf - 1, maxconhist))) | infnan_obj.is_posinf(conhist(:, 1:min(nf - 1, maxconhist))), 'all'), "CONHIST does not contain NaN/Inf", srname);
                end
                % Check the values of X, F, CSTRV, CONSTR.
                % X does not contain NaN if X0 does not and the trust-region/geometry steps are proper.
                debug_obj.assert(~any(infnan_obj.is_nan_sp(x), 'all'), "X does not contain NaN", srname);
                % F cannot be NaN/+Inf due to the moderated extreme barrier.
                debug_obj.assert(~(infnan_obj.is_nan_sp(f) || infnan_obj.is_posinf(f)), "F is not NaN/+Inf", srname);
                if ~ismember('cstrv', ipObj.UsingDefaults)
                    debug_obj.assert(~(cstrv < 0), "CSTRV is not negative", srname);
                    %------------------------------------------------------------------------------------------%
                    % The following test is not applicable to LINCOA.
                    % %call assert(.not. (is_nan(cstrv) .or. is_posinf(cstrv)), 'CSTRV is NaN/+Inf', srname)
                    %------------------------------------------------------------------------------------------%

                end
                if ~ismember('constr', ipObj.UsingDefaults)
                    % CONSTR cannot contain NaN/+Inf due to the moderated extreme barrier.
                    debug_obj.assert(~any(infnan_obj.is_nan_sp(constr) | infnan_obj.is_posinf(constr), 'all'), "CONSTR does not contain NaN/+Inf", srname);
                end
            end

            %====================%
            % Calculation starts %
            %====================%

            % Save the history. Note that NF may exceed the maximal amount of history to save. We save X and F
            % at the position indexed by MODULO(NF - 1, MAXHIST) + 1. When the solver terminates, the history
            % will be reordered so that the information is in the chronological order. Similar for CONSTR, CSTRV.
            if maxxhist > 0
                % We could replace MODULO(NF - 1_IK, MAXXHIST) + 1_IK) with MODULO(NF - 1_IK, MAXHIST) + 1_IK)
                % based on the assumption that MAXXHIST == 0 or MAXHIST. For robustness, we do not do that.
                xhist(:, mod(nf - 1, maxxhist) + 1) = x;
            end
            if maxfhist > 0
                fhist(mod(nf - 1, maxfhist) + 1) = f;
            end
            if maxchist > 0
                % MAXCHIST > 0 implies PRESENT(CHIST) and PRESENT(CSTRV)
                chist(mod(nf - 1, maxchist) + 1) = cstrv;
            end
            if maxconhist > 0
                % MAXCONHIST > 0 implies PRESENT(CONHIST) and PRESENT (CONSTR)
                conhist(:, mod(nf - 1, maxconhist) + 1) = constr;
            end

            %====================%
            %  Calculation ends  %
            %====================%

            % Postconditions
            if consts_obj.DEBUGGING
                % Called after each function evaluation when debugging; can be expensive.
                debug_obj.assert(size(xhist, 1) == numel(x) && size(xhist, 2) == maxxhist, "SIZE(XHIST) == [SIZE(X), MAXXHIST]", srname);
                debug_obj.assert(~any(infnan_obj.is_nan_sp(xhist(:, 1:min(nf, maxxhist))), 'all'), "XHIST does not contain NaN", srname);
                % The last calculated X can be Inf (finite + finite can be Inf numerically).
                debug_obj.assert(numel(fhist) == maxfhist, "SIZE(FHIST) == MAXFHIST", srname);
                debug_obj.assert(~any(infnan_obj.is_nan_sp(fhist(1:min(nf, maxfhist))) | infnan_obj.is_posinf(fhist(1:min(nf, maxfhist))), 'all'), "FHIST does not contain NaN/+Inf", srname);
                if ~ismember('chist', ipObj.UsingDefaults) || nargout >= 3
                    debug_obj.assert(numel(chist) == maxchist, "SIZE(CHIST) == MAXCHIST", srname);
                    debug_obj.assert(~any(chist(1:min(nf, maxchist)) < 0, 'all'), "CHIST does not contain negative values", srname);
                    %------------------------------------------------------------------------------------------%
                    % The following test is not applicable to LINCOA.
                    % %call assert(.not. any(is_nan(chist(1:min(nf, maxchist))) .or. is_posinf(chist(1:min(nf, maxchist)))), &
                    % %    & 'CHIST does not contain NaN/+Inf', srname)
                    %------------------------------------------------------------------------------------------%

                end
                if (~ismember('conhist', ipObj.UsingDefaults) || nargout >= 4) && ~ismember('constr', ipObj.UsingDefaults)
                    debug_obj.assert(size(conhist, 1) == numel(constr) && size(conhist, 2) == maxconhist, "SIZE(CONHIST) == [SIZE(CONSTR), MAXCONHIST]", srname);
                    debug_obj.assert(~any(infnan_obj.is_nan_sp(conhist(:, 1:min(nf, maxconhist))) | infnan_obj.is_posinf(conhist(:, 1:min(nf, maxconhist))), 'all'), "CONHIST does not contain NaN/+Inf", srname);
                end

                % The following code checks that XHIST does not contain a segment that repeats. If such a segment
                % is found, we believe that the solver has encountered an infinite cycle, which would be a bug.
                % N.B.:
                % 1. We check this only if NF > (N+1)*(N+2)/2, when the initialization has surely finished. This
                % is because XHIST may contain repeating segments during the initialization if X0 + RHOBEG = X0,
                % which can happen if all entries of X0 are excessively large compared with RHOBEG. It is
                % possible to revise the initialization subroutine to avoid repetition, but we choose not to,
                % a motivation being to keep the initialization parallelizable.
                % 2. We skip the test if N = 1, as false positive may occur (also possible when N > 1, but rare).
                % 3. For segments of length 1, we check whether it repeats three times. For segments of length
                % i > 1, we check whether it repeats twice. Due to rounding errors, it may happen that the same
                % point is repeated twice, but the solver is not in an infinite cycle, which was observed in an
                % experiment of NEWUOA on 20240404.
                nhist = min(nf, maxxhist);
                n = numel(x);
                if n > 1 && nf > (n + 1) * (n + 2) / 2
                    if nhist >= 3
                        debug_obj.wassert(~(all(abs(xhist(:, nhist) - xhist(:, nhist - 1)) <= 0, 'all') && all(abs(xhist(:, nhist - 1) - xhist(:, nhist - 2)) <= 0, 'all')), "XHIST does not contain a repeating segment of length 1", srname);
                    end
                    for i = 2:min(100, nhist / 2)
                        debug_obj.wassert(~all(abs(xhist(:, nhist - i + 1:nhist) - xhist(:, nhist - 2 * i + 1:nhist - i)) <= 0, 'all'), "XHIST does not contain a repeating segment of length " + string_obj.int2str(i), srname);
                    end
                end
            end

        end
        function [xhist, fhist, chist, conhist] = rangehist(~, nf, xhist, fhist, varargin)
            %--------------------------------------------------------------------------------------------------%
            % This subroutine arranges FHIST, XHIST, CHIST, and CONHIST in the chronological order.
            %--------------------------------------------------------------------------------------------------%
            consts_obj = prima_mat.common.consts_mod();
            infnan_obj = prima_mat.common.infnan_mod();
            debug_obj = prima_mat.common.debug_mod();

            % Inputs


            % In-outputs



            % Local variables



            srname = "RANGEHIST";

            % Sizes
            n = size(xhist, 1);
            maxxhist = size(xhist, 2);
            maxfhist = numel(fhist);
            ipObj = inputParser();
            addParameter(ipObj, 'chist', NaN);
            addParameter(ipObj, 'conhist', NaN);
            parse(ipObj, varargin{:});
            chist = ipObj.Results.chist;
            conhist = ipObj.Results.conhist;
            if ~ismember('chist', ipObj.UsingDefaults) || nargout >= 3
                maxchist = numel(chist);
            else
                maxchist = 0;
            end
            if ~ismember('conhist', ipObj.UsingDefaults) || nargout >= 4
                m = size(conhist, 1);
                maxconhist = size(conhist, 2);
            else
                m = 0;
                maxconhist = 0;
            end
            maxhist = max(maxxhist, max(maxfhist, max(maxconhist, maxchist)));

            % Preconditions
            if consts_obj.DEBUGGING
                % Check the sizes of XHIST, FHIST, CHIST, CONHIST.
                debug_obj.assert(n >= 1, "SIZE(XHIST, 1) >= 1", srname);
                debug_obj.assert(maxxhist * (maxxhist - maxhist) == 0, "SIZE(XHIST, 2) == 0 or MAXHIST", srname);
                debug_obj.assert(maxfhist * (maxfhist - maxhist) == 0, "SIZE(FHIST) == 0 or MAXHIST", srname);
                debug_obj.assert(maxchist * (maxchist - maxhist) == 0, "SIZE(CHIST) == 0 or MAXHIST", srname);
                debug_obj.assert(maxconhist * (maxconhist - maxhist) == 0, "SIZE(CONHIST, 2) == 0 or MAXHIST", srname);
                % Check the values of XHIST, FHIST, CHIST, CONHIST.
                debug_obj.assert(~any(infnan_obj.is_nan_sp(xhist(:, 1:min(nf, maxxhist))), 'all'), "XHIST does not contain NaN", srname);
                % The last calculated X can be Inf (finite + finite can be Inf numerically).
                debug_obj.assert(~any(infnan_obj.is_nan_sp(fhist(1:min(nf, maxfhist))) | infnan_obj.is_posinf(fhist(1:min(nf, maxfhist))), 'all'), "FHIST does not contain NaN/+Inf", srname);
                if ~ismember('chist', ipObj.UsingDefaults) || nargout >= 3
                    debug_obj.assert(~any(chist(1:min(nf, maxchist)) < 0, 'all'), "CHIST does not contain negative values", srname);
                    %------------------------------------------------------------------------------------------%
                    % The following test is not applicable to LINCOA.
                    % %call assert(.not. any(is_nan(chist(1:min(nf, maxchist))) .or. is_posinf(chist(1:min(nf, maxchist)))), &
                    % %    & 'CHIST does not contain NaN/+Inf', srname)
                    %------------------------------------------------------------------------------------------%

                end
                if ~ismember('conhist', ipObj.UsingDefaults) || nargout >= 4
                    debug_obj.assert(~any(infnan_obj.is_nan_sp(conhist(:, 1:min(nf, maxconhist))) | infnan_obj.is_posinf(conhist(:, 1:min(nf, maxconhist))), 'all'), "CONHIST does not contain NaN/+Inf", srname);
                end
            end

            %====================%
            % Calculation starts %
            %====================%

            % The ranging should be done only if 0 < MAXXHIST < NF. Otherwise, it leads to errors/wrong results.
            if maxxhist > 0 && maxxhist < nf
                % We could replace MODULO(NF - 1_IK, MAXXHIST) + 1_IK) with MODULO(NF - 1_IK, MAXHIST) + 1_IK)
                % based on the assumption that MAXXHIST == 0 or MAXHIST. For robustness, we do not do that.
                khist = mod(nf - 1, maxxhist) + 1;
                xhist(:, :) = reshape([reshape(xhist(:, khist + 1:maxxhist), 1, []), reshape(xhist(:, 1:khist), 1, [])], size(xhist));
                % N.B.:
                % 1. The result of the array constructor is always a rank-1 array (e.g., vector), no matter what
                % elements are used for the construction.
                % 2. The above combination of SHAPE and RESHAPE fulfills our desire thanks to the COLUMN-MAJOR
                % order of Fortran arrays.
                % 3. In MATLAB, `xhist = [xhist(:, khist + 1:maxxhist), xhist(:, 1:khist)]` does the same thing.

            end
            % The ranging should be done only if 0 < MAXFHIST < NF. Otherwise, it leads to errors/wrong results.
            if maxfhist > 0 && maxfhist < nf
                khist = mod(nf - 1, maxfhist) + 1;
                fhist(:) = [fhist(khist + 1:maxfhist); fhist(1:khist)];
            end
            % The ranging should be done only if 0 < MAXCONHIST < NF. Otherwise, it leads to errors/wrong results.
            if maxconhist > 0 && maxconhist < nf
                khist = mod(nf - 1, maxconhist) + 1;
                conhist(:, :) = reshape([reshape(conhist(:, khist + 1:maxconhist), 1, []), reshape(conhist(:, 1:khist), 1, [])], size(conhist));
            end
            % The ranging should be done only if 0 < MAXCHIST < NF. Otherwise, it leads to errors/wrong results.
            if maxchist > 0 && maxchist < nf
                khist = mod(nf - 1, maxchist) + 1;
                chist(:) = [chist(khist + 1:maxchist); chist(1:khist)];
            end

            %====================%
            %  Calculation ends  %
            %====================%

            % Postconditions
            if consts_obj.DEBUGGING
                debug_obj.assert(size(xhist, 1) == n && size(xhist, 2) == maxxhist, "SIZE(XHIST) == [N, MAXXHIST]", srname);
                debug_obj.assert(~any(infnan_obj.is_nan_sp(xhist(:, 1:min(nf, maxxhist))), 'all'), "XHIST does not contain NaN", srname);
                % The last calculated X can be Inf (finite + finite can be Inf numerically).
                debug_obj.assert(numel(fhist) == maxfhist, "SIZE(FHIST) == MAXFHIST", srname);
                debug_obj.assert(~any(infnan_obj.is_nan_sp(fhist(1:min(nf, maxfhist))) | infnan_obj.is_posinf(fhist(1:min(nf, maxfhist))), 'all'), "FHIST does not contain NaN/+Inf", srname);
                if ~ismember('chist', ipObj.UsingDefaults) || nargout >= 3
                    debug_obj.assert(numel(chist) == maxchist, "SIZE(CHIST) == MAXCHIST", srname);
                    debug_obj.assert(~any(chist(1:min(nf, maxchist)) < 0, 'all'), "CHIST does not contain negative values", srname);
                    %------------------------------------------------------------------------------------------%
                    % The following test is not applicable to LINCOA.
                    % %call assert(.not. any(is_nan(chist(1:min(nf, maxchist))) .or. is_posinf(chist(1:min(nf, maxchist)))), &
                    % %    & 'CHIST does not contain NaN/+Inf', srname)
                    %------------------------------------------------------------------------------------------%

                end
                if ~ismember('conhist', ipObj.UsingDefaults) || nargout >= 4
                    debug_obj.assert(size(conhist, 1) == m && size(conhist, 2) == maxconhist, "SIZE(CONHIST) == [M, MAXCONHIST]", srname);
                    debug_obj.assert(~any(infnan_obj.is_nan_sp(conhist(:, 1:min(nf, maxconhist))) | infnan_obj.is_posinf(conhist(:, 1:min(nf, maxconhist))), 'all'), "CONHIST does not contain NaN/+Inf", srname);
                end
            end

        end

    end
end