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
        function varargout = num2str_custom(obj, varargin)
            if numel(varargin) == 1 && (isinteger(varargin{1}) || isnumeric(varargin{1}) && (isreal(varargin{1}) && all(fix(varargin{1}) == varargin{1}, 'all'))) && isscalar(varargin{1})
                [varargout{1:nargout}] = int2str(varargin{:});
            elseif numel(varargin) >= 1 && numel(varargin) <= 3 && isfloat(varargin{1}) && isscalar(varargin{1})
                [varargout{1:nargout}] = obj.real2str_scalar(varargin{:});
            else
                [varargout{1:nargout}] = obj.real2str_vector(varargin{:});
            end
        end
        function y = upper(~, x)
            %--------------------------------------------------------------------------------------------------%
            % This function maps the characters of a string to the upper case, if applicable.
            %--------------------------------------------------------------------------------------------------%


            y = pad(" ", strlength(x));

            dist = 'A' - 'a';

            y = x;
            for i = 1:strlength(y)
                if extractBetween(y, i, i) >= "a" && extractBetween(y, i, i) <= "z"
                    y = replaceBetween(y, i, i, char(double(unicode2native(extractBetween(y, i, i))) + dist));
                end
            end
        end
        function y = istr(~, x)
            %--------------------------------------------------------------------------------------------------%
            % This function converts a string to an integer array.
            %--------------------------------------------------------------------------------------------------%


            y = NaN(strlength(x), 1);

            y(:) = arrayfun(@(i) fix(double(unicode2native(extractBetween(x, i, i)))), (1:fix(strlength(x)))');

        end
        function s = real2str_scalar(obj, x, varargin)
            %--------------------------------------------------------------------------------------------------%
            % This function converts a real scalar to a string. Optionally, NDGT is the number of decimal
            % digits to print, and NEXP is the number of digits in the exponent; they may be reduced if needed.
            %--------------------------------------------------------------------------------------------------%


            s = "";

            % The number of decimal digits to print
            % The number of digits in the exponent
            % The width of the printed X


            ipObj = inputParser();
            addParameter(ipObj, 'ndgt', NaN);
            addParameter(ipObj, 'nexp', NaN);
            parse(ipObj, varargin{:});
            ndgt = ipObj.Results.ndgt;
            nexp = ipObj.Results.nexp;

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

            if isfinite(x)
                wx = ndgt_loc + nexp_loc + 5;
                if ~(wx <= obj.MAX_NUM_STR_LEN)
                    error("The width of the printed number is at most " + int2str(obj.MAX_NUM_STR_LEN));
                end

                str = sprintf('%s \n', num2str(x));
                s = strtrim(str); % Remove the trailing spaces, but keep the leading ones, if any.
            else
                str = sprintf('%s \n', num2str(x));
                s = strip(str); % Remove the leading and trailing spaces, if any.

            end

            %====================%
            %  Calculation ends  %
            %====================%


        end
        function s = real2str_vector(obj, x, varargin)
            %--------------------------------------------------------------------------------------------------%
            % This function converts a real vector to a string. Optionally, NDGT is the number of decimal
            % digits to print, NEXP is the number of digits in the exponent, and NX is the number of entries
            % printed per row; they may be reduced if needed.
            %--------------------------------------------------------------------------------------------------%


            s = "";

            spaces = "  "; % The spaces between two entries in a row


            % The number of rows
            % N = SIZE(X)
            % The number of decimal digits to print
            % The number of digits in the exponent
            % The number of entries printed per row
            % The length of the string
            % The width of each entry in X


            ipObj = inputParser();
            addParameter(ipObj, 'ndgt', NaN);
            addParameter(ipObj, 'nexp', NaN);
            addParameter(ipObj, 'nx', NaN);
            parse(ipObj, varargin{:});
            ndgt = ipObj.Results.ndgt;
            nexp = ipObj.Results.nexp;
            nx = ipObj.Results.nx;

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
            s = repmat("", [slen, 1]);

            j = 0; % J is the index of the last up-to-date character in S.
            for i = 1:n
                s = replaceBetween(s, j + 1, j + wx, obj.real2str_scalar(x(i), 'ndgt', ndgt_loc, 'nexp', nexp_loc));
                if i == n
                    break
                end
                j = j + wx;
                if mod(i, nx_loc) == 0
                    s = replaceBetween(s, j + 1, j + 1, newline);
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


            x = NaN;

            x = sscanf(s, '%s');
        end
        function x = str2int(~, s)
            %--------------------------------------------------------------------------------------------------%
            % This function converts a string to an integer scalar.
            %--------------------------------------------------------------------------------------------------%


            x = NaN;

            x = sscanf(s, '%s');
        end

    end
end