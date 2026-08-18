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



            % Inputs



            % In-outputs


            % Outputs



            % Local variables

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


            %====================%
            % Calculation starts %
            %====================%

            % Save the input value of MAXHIST for debugging.


            % Revise MAXHIST according to MAXHISTMEM, i.e., the maximal memory allowed for the history.
            % N.B.: The `UNIT_MEMO = INT(*)` below converts integers to the default integer kind, which is the
            % kind of UNIT_MEMO. Fortran compilers may complain without the conversion. It is not needed in
            % Python/MATLAB/Julia/R. Meanwhile, INT(OUTPUT_*HIST) converts booleans to integers.
            unit_memo = fix(output_xhist) * n + fix(output_fhist);
            if ~ismember('output_chist', ipObj.UsingDefaults) && nargout >= 4
                unit_memo = unit_memo + fix(output_chist);
            end
            if ~ismember('m', ipObj.UsingDefaults) && ~ismember('output_conhist', ipObj.UsingDefaults) && nargout >= 5
                unit_memo = unit_memo + fix(output_conhist) * m;
            end
            unit_memo = unit_memo * fix(8); % INT(*) avoids overflow when IK is 16-bit.
            if unit_memo <= 0
                % No output of history is requested
                maxhist = 0;
            elseif maxhist > 100000000 / unit_memo
                maxhist = fix(100000000 / unit_memo); % Integer division.
                % We cannot simply set MAXHIST = MIN(MAXHIST, MAXHISTMEM/UNIT_MEMO), as they may not have
                % the same kind, and compilers may complain. We may convert them, but overflow may occur.

            end

            xhist = NaN(n, maxhist * fix(output_xhist));
            fhist = NaN(maxhist * fix(output_fhist), 1);
            % Even if OUTPUT_CHIST is FALSE, CHIST still needs to be allocated.
            if ~ismember('output_chist', ipObj.UsingDefaults) && nargout >= 4
                chist = NaN(maxhist * fix(output_chist), 1);
            end
            % Even if OUTPUT_CONHIST is FALSE, CONHIST still needs to be allocated.
            if ~ismember('m', ipObj.UsingDefaults) && ~ismember('output_conhist', ipObj.UsingDefaults) && nargout >= 5
                conhist = NaN(m, maxhist * fix(output_conhist));
            end

            %====================%
            %  Calculation ends  %
            %====================%

            % Postconditions

        end
        function [xhist, fhist, chist, conhist] = savehist(~, nf, x, xhist, f, fhist, varargin)
            %--------------------------------------------------------------------------------------------------%
            % This subroutine saves X, F, CSTRV, and CONSTR into XHIST, FHIST, CHIST, and CONHIST respectively.
            %--------------------------------------------------------------------------------------------------%



            % Inputs



            % In-outputs



            % Local variables



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


            % Preconditions


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


        end
        function [xhist, fhist, chist, conhist] = rangehist(~, nf, xhist, fhist, varargin)
            %--------------------------------------------------------------------------------------------------%
            % This subroutine arranges FHIST, XHIST, CHIST, and CONHIST in the chronological order.
            %--------------------------------------------------------------------------------------------------%



            % Inputs


            % In-outputs



            % Local variables



            % Sizes

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

                maxconhist = size(conhist, 2);
            else

                maxconhist = 0;
            end


            % Preconditions


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


        end

    end
end