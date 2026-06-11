classdef string_mod
    %--------------------------------------------------------------------------------------------------%
    % This module provides some procedures for manipulating strings.
    %
    % Coded by Zaikun ZHANG (www.zhangzk.net).
    %
    % Started: September 2021
    %
    % Last Modified: Sunday, March 31, 2024 PM09:44:38
    %--------------------------------------------------------------------------------------------------%
    properties (Access = private)
        MAX_NUM_STR_LEN;
        MAX_WIDTH;
    end

    methods
        function obj = string_mod()
            % MAX_NUM_STR_LEN is the maximum length of a string that is needed to represent a real or integer
            % number. Assuming that such a number is represented by at most 128 bits, it is safe to set this
            % maximum length to 128. We set this number to 1024 to be on the safe side.
            obj.MAX_NUM_STR_LEN = 1024;
            % MAX_WIDTH is the maximum number of characters printed in each row when printing arrays.
            obj.MAX_WIDTH = 100;
        end
        function varargout = num2str(obj, varargin)
            if numel(varargin) == 1 && isinteger(varargin{1}) && isscalar(varargin{1})
                [varargout{1:nargout}] = obj.int2str(varargin{:});
            elseif numel(varargin) >= 1 && numel(varargin) <= 3 && isfloat(varargin{1}) && isscalar(varargin{1})
                [varargout{1:nargout}] = obj.real2str_scalar(varargin{:});
            else
                [varargout{1:nargout}] = obj.real2str_vector(varargin{:});
            end
        end
        function y = lower(~, x)
            %--------------------------------------------------------------------------------------------------%
            % This function maps the characters of a string to the lower case, if applicable.
            %--------------------------------------------------------------------------------------------------%



            y = pad(" ", strlength(x));

            dist = unicode2native("A") - unicode2native("a");
            i = NaN;

            y = x;
            for i = 1:strlength(y)
                if extract(y, i) >= "A" && extract(y, i) <= "Z"
                    y = replaceBetween(y, i, i, char(unicode2native(extract(y, i)) - dist));
                end
            end
        end
        function y = upper(~, x)
            %--------------------------------------------------------------------------------------------------%
            % This function maps the characters of a string to the upper case, if applicable.
            %--------------------------------------------------------------------------------------------------%



            y = pad(" ", strlength(x));

            dist = unicode2native("A") - unicode2native("a");
            i = NaN;

            y = x;
            for i = 1:strlength(y)
                if extract(y, i) >= "a" && extract(y, i) <= "z"
                    y = replaceBetween(y, i, i, char(unicode2native(extract(y, i)) + dist));
                end
            end
        end
        function y = strip(~, x)
            %--------------------------------------------------------------------------------------------------%
            % This function removes the leading and trailing spaces of a string.
            %--------------------------------------------------------------------------------------------------%



            y = pad(" ", strlength(strtrim(strjust(x, 'left'))));

            y = strtrim(strjust(x, 'left'));
        end
        function y = istr(~, x)
            %--------------------------------------------------------------------------------------------------%
            % This function converts a string to an integer array.
            %--------------------------------------------------------------------------------------------------%
            consts_obj = consts_mod();

            y = NaN(strlength(x), 1);

            i = NaN;

            y(:) = reshape(cell2mat(arrayfun(@(i) fix(unicode2native(extract(x, i))), (1:fix(strlength(x))), "UniformOutput", false)), [], 1);

        end
        function s = real2str_scalar(obj, x, varargin)
            %--------------------------------------------------------------------------------------------------%
            % This function converts a real scalar to a string. Optionally, NDGT is the number of decimal
            % digits to print, and NEXP is the number of digits in the exponent; they may be reduced if needed.
            %--------------------------------------------------------------------------------------------------%
            consts_obj = consts_mod();
            debug_obj = debug_mod();
            infnan_obj = infnan_mod();
            % Inputs



            % Outputs
            s = "";
            % Local variables
            srname = "REAL2STR_SCALAR";
            sformat = "";
            str = pad(" ", obj.MAX_NUM_STR_LEN);
            ndgt_loc = NaN; % The number of decimal digits to print
            nexp_loc = NaN; % The number of digits in the exponent
            wx = NaN; % The width of the printed X

            % Preconditions
            ipObj = inputParser();
            addParameter(ipObj, 'ndgt', NaN);
            addParameter(ipObj, 'nexp', NaN);
            parse(ipObj, varargin{:});
            ndgt = ipObj.Results.ndgt;
            nexp = ipObj.Results.nexp;
            if consts_obj.DEBUGGING
                if ~ismember('ndgt', ipObj.UsingDefaults)
                    debug_obj.assert(ndgt >= 0 && 2 * ndgt <= obj.MAX_NUM_STR_LEN - 5, "0 <= NDGT <= " + obj.int2str(floor(double(obj.MAX_NUM_STR_LEN - 5) / 2.0)), srname);
                end
                if ~ismember('nexp', ipObj.UsingDefaults)
                    debug_obj.assert(nexp >= 0 && 2 * nexp <= obj.MAX_NUM_STR_LEN - 5, "0 <= NEXP <= " + obj.int2str(floor(double(obj.MAX_NUM_STR_LEN - 5) / 2.0)), srname);
                end
            end

            %====================%
            % Calculation starts %
            %====================%

            if ismember('ndgt', ipObj.UsingDefaults)
                % By default, we print at most the same number of decimal digits as the double precision.
                ndgt_loc = min(floor(-log10(eps(class(x)))), floor(-log10(eps(class(0.0))))) + 1;
            else
                ndgt_loc = ndgt;
            end
            ndgt_loc = min(ndgt_loc, floor(double(obj.MAX_NUM_STR_LEN - 5) / 2.0)); % Safeguard
            if ismember('nexp', ipObj.UsingDefaults)
                nexp_loc = ceil(log10(double(floor(log10(realmax(class(x)))) + 0.1))); % Use + 0.1 in case RANGE(X) = 10^k.
            else
                nexp_loc = nexp;
            end
            nexp_loc = min(nexp_loc, floor(double(obj.MAX_NUM_STR_LEN - 5) / 2.0));

            if infnan_obj.is_finite(x)
                wx = ndgt_loc + nexp_loc + 5;
                debug_obj.validate(wx <= obj.MAX_NUM_STR_LEN, "The width of the printed number is at most " + obj.int2str(fix(obj.MAX_NUM_STR_LEN)), srname);
                sformat = "(1PE" + obj.int2str(fix(wx)) + "." + obj.int2str(fix(ndgt_loc)) + "E" + obj.int2str(fix(nexp_loc)) + ")";
                str = sprintf('%s \n', obj.num2str(x));
                s = strtrim(str); % Remove the trailing spaces, but keep the leading ones, if any.
            else
                str = sprintf('%s \n', obj.num2str(x));
                s = obj.strip(str); % Remove the leading and trailing spaces, if any.

            end

            %====================%
            %  Calculation ends  %
            %====================%

            % Postconditions
            if consts_obj.DEBUGGING
                debug_obj.assert(strlength(s) > 0 && strlength(s) <= obj.MAX_NUM_STR_LEN, "0 < LEN(S) <= MAX_NUM_STR_LEN", srname);
                debug_obj.assert(infnan_obj.is_nan_sp(x) == infnan_obj.is_nan_sp(obj.str2real(s)), "IS_NAN(X) .EQV. IS_NAN(STR2REAL(S))", srname);
                % The assertions concerning the infiniteness of X may fail due to the limited precision of
                % printing. Thus we relax the assertions as below.
                %call assert(is_posinf(x) .eqv. is_posinf(str2real(s)), 'IS_POSINF(X) .EQV. IS_POSINF(STR2REAL(S))', srname)
                %call assert(is_neginf(x) .eqv. is_neginf(str2real(s)), 'IS_NEGINF(X) .EQV. IS_NEGINF(STR2REAL(S))', srname)
                debug_obj.assert((x >= consts_obj.REALMAX * (1.0 - fortran.power(10.0, (-ndgt_loc)))) == (obj.str2real(s) >= consts_obj.REALMAX * (1.0 - fortran.power(10.0, (-ndgt_loc)))), "IS_POSINF(X) .EQV. IS_POSINF(STR2REAL(S))", srname);
                debug_obj.assert((x <= -consts_obj.REALMAX * (1.0 - fortran.power(10.0, (-ndgt_loc)))) == (obj.str2real(s) <= -consts_obj.REALMAX * (1.0 - fortran.power(10.0, (-ndgt_loc)))), "IS_NEGINF(X) .EQV. IS_NEGINF(STR2REAL(S))", srname);
                if abs(x) < consts_obj.REALMAX
                    debug_obj.assert(abs(x - obj.str2real(s)) <= abs(x) * fortran.power(10.0, (-ndgt_loc)), "STR2REAL(S) == X", srname);
                end
            end
        end
        function s = real2str_vector(obj, x, varargin)
            %--------------------------------------------------------------------------------------------------%
            % This function converts a real vector to a string. Optionally, NDGT is the number of decimal
            % digits to print, NEXP is the number of digits in the exponent, and NX is the number of entries
            % printed per row; they may be reduced if needed.
            %--------------------------------------------------------------------------------------------------%
            consts_obj = consts_mod();
            debug_obj = debug_mod();
            memory_obj = memory_mod();
            % Inputs



            % Outputs
            s = "";
            % Local variables
            srname = "REAL2STR_VECTOR";
            spaces = "  "; % The spaces between two entries in a row
            i = NaN;
            j = NaN;
            m = NaN; % The number of rows
            n = NaN; % N = SIZE(X)
            ndgt_loc = NaN; % The number of decimal digits to print
            nexp_loc = NaN; % The number of digits in the exponent
            nx_loc = NaN; % The number of entries printed per row
            slen = NaN; % The length of the string
            wx = NaN; % The width of each entry in X

            % Preconditions
            ipObj = inputParser();
            addParameter(ipObj, 'ndgt', NaN);
            addParameter(ipObj, 'nexp', NaN);
            addParameter(ipObj, 'nx', NaN);
            parse(ipObj, varargin{:});
            ndgt = ipObj.Results.ndgt;
            nexp = ipObj.Results.nexp;
            nx = ipObj.Results.nx;
            if consts_obj.DEBUGGING
                if ~ismember('ndgt', ipObj.UsingDefaults)
                    debug_obj.assert(ndgt >= 0 && 2 * ndgt <= obj.MAX_NUM_STR_LEN - 5, "0 <= NDGT <= " + obj.int2str(floor(double(obj.MAX_NUM_STR_LEN - 5) / 2.0)), srname);
                end
                if ~ismember('nexp', ipObj.UsingDefaults)
                    debug_obj.assert(nexp >= 0 && 2 * nexp <= obj.MAX_NUM_STR_LEN - 5, "0 <= NEXP <= " + obj.int2str(floor(double(obj.MAX_NUM_STR_LEN - 5) / 2.0)), srname);
                end
                if ~ismember('nx', ipObj.UsingDefaults)
                    debug_obj.assert(nx >= 1, "NX >= 1", srname);
                end
            end

            %====================%
            % Calculation starts %
            %====================%

            % Quick return if X is empty.
            if numel(x) <= 0
                s = "";
                return
            end

            if ismember('ndgt', ipObj.UsingDefaults)
                % By default, we print at most the same number of decimal digits as the double precision.
                ndgt_loc = min(floor(-log10(eps(class(x)))), floor(-log10(eps(class(0.0))))) + 1;
            else
                ndgt_loc = ndgt;
            end
            ndgt_loc = min(ndgt_loc, floor(double(obj.MAX_NUM_STR_LEN - 5) / 2.0)); % Safeguard

            if ismember('nexp', ipObj.UsingDefaults)
                nexp_loc = ceil(log10(double(floor(log10(realmax(class(x))))) + 0.1)); % Use + 0.1 in case RANGE(X) = 10^k.
            else
                nexp_loc = nexp;
            end
            nexp_loc = min(nexp_loc, floor(double(obj.MAX_NUM_STR_LEN - 5) / 2.0));

            wx = strlength(obj.real2str_scalar(0.0, 'ndgt', ndgt, 'nexp', nexp));
            n = numel(x);
            if ismember('nx', ipObj.UsingDefaults)
                nx_loc = max(1, min(floor(double(obj.MAX_WIDTH + strlength(spaces)) / (double(wx) + strlength(spaces))), numel(x)));
            else
                nx_loc = max(1, min(nx, n));
            end

            % Calculate the length of the printed string S.
            % N.B.: Here, SLEN should not be an INT16 integer, because a double-precision vector of length
            % ~3500 would be printed as a string longer than 65536. On most modern platforms, the default
            % integer kind is INT32, which is enough for printing double-precision vectors of size ~ 10^8,
            % being sufficient for this project.
            m = ceil(double(n) / double(nx_loc)); % The number of rows
            slen = wx * n + strlength(spaces) * (n - 1) + (1 - strlength(spaces)) * (m - 1);
            s = memory_obj.alloc_character(slen);

            j = 0; % J is the index of the last up-to-date character in S.
            for i = 1:n
                s = replaceBetween(s, j + 1, j + wx, obj.real2str_scalar(x(i), 'ndgt', ndgt_loc, 'nexp', nexp_loc));
                if i == n
                    break
                end
                j = j + wx;
                if mod(i, nx_loc) == 0
                    s = replaceBetween(s, j + 1, j + 1, compose('\n'));
                    j = j + 1;
                else
                    s = replaceBetween(s, j + 1, j + strlength(spaces), spaces);
                    j = j + strlength(spaces);
                end
            end

            %====================%
            %  Calculation ends  %
            %====================%
        end
        function x = str2real(~, s)
            %--------------------------------------------------------------------------------------------------%
            % This function converts a string to a real scalar.
            %--------------------------------------------------------------------------------------------------%
            consts_obj = consts_mod();
            debug_obj = debug_mod();

            x = NaN;
            srname = "STR2REAL";
            if consts_obj.DEBUGGING
                debug_obj.assert(strlength(s) > 0, "LEN(S) > 0", srname);
            end
            x = sscanf(s, '%s');
        end
        function s = int2str(obj, x)
            %--------------------------------------------------------------------------------------------------%
            % This function converts an integer scalar to a string.
            %--------------------------------------------------------------------------------------------------%
            consts_obj = consts_mod();
            debug_obj = debug_mod();

            srname = "INT2STR";
            s = "";
            str = pad(" ", obj.MAX_NUM_STR_LEN);
            % In the following, 'I0' means to use the minimum number of digits needed to print.
            % It should work also if we use * instead of I0. However, this sometimes lead to a segmentation
            % fault on Windows Server 2022 with gcc/gfortran 13.
            str = sprintf('%d\n', x);
            s = obj.strip(str);
            if consts_obj.DEBUGGING
                debug_obj.assert(strlength(s) > 0 && strlength(s) <= obj.MAX_NUM_STR_LEN, "0 < LEN(S) <= MAX_NUM_STR_LEN", srname);
                debug_obj.assert(obj.str2int(s) == x, "STR2INT(S) == X", srname);
            end
        end
        function x = str2int(~, s)
            %--------------------------------------------------------------------------------------------------%
            % This function converts a string to an integer scalar.
            %--------------------------------------------------------------------------------------------------%
            consts_obj = consts_mod();
            debug_obj = debug_mod();

            x = NaN;
            srname = "STR2INT";
            if consts_obj.DEBUGGING
                debug_obj.assert(strlength(s) > 0, "LEN(S) > 0", srname);
            end
            x = sscanf(s, '%s');
        end

    end
end