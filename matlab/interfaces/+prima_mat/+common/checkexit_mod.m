%TODO: merge CHECKEXIT_UNC and CHECKEXIT_CON, using optional CSTRV and CTOL.
classdef checkexit_mod
    %--------------------------------------------------------------------------------------------------%
    % This module checks whether to exit the solver.
    %
    % Coded by Zaikun ZHANG (www.zhangzk.net).
    %
    % Started: September 2021
    %
    % Last Modified: Tuesday, September 26, 2023 AM10:51:16
    %--------------------------------------------------------------------------------------------------%

    methods
        function varargout = checkexit(obj, varargin)
            if numel(varargin) == 7 && isscalar(varargin{5})
                [varargout{1:nargout}] = obj.checkexit_con(varargin{:});
            else
                [varargout{1:nargout}] = obj.checkexit_unc(varargin{:});
            end
        end
        function info = checkexit_unc(~, maxfun, nf, f, ftarget, x)
            %--------------------------------------------------------------------------------------------------%
            % This module checks whether to exit the solver in the unconstrained case.
            %--------------------------------------------------------------------------------------------------%


            info = NaN;

            %====================%
            % Calculation starts %
            %====================%

            info = 0; % Default info, indicating that the solver should not exit.

            % Although X should not contain NaN unless there is a bug, we include the following for security.
            % X can be Inf, as finite + finite can be Inf numerically.
            if any(isnan(x) | isinf(x), 'all')
                info = -1;
            end

            % Although NAN_INF_F should not happen unless there is a bug, we include the following for security.
            if isnan(f) | isinf(f) & f > 0
                info = -2;
            end

            if f <= ftarget
                info = 1;
            end

            if nf >= maxfun
                info = 3;
            end

            %====================%
            %  Calculation ends  %
            %====================%


        end
        function info = checkexit_con(~, maxfun, nf, cstrv, ctol, f, ftarget, x)
            %--------------------------------------------------------------------------------------------------%
            % This module checks whether to exit the solver in the constrained case.
            %--------------------------------------------------------------------------------------------------%


            info = NaN;

            %====================%
            % Calculation starts %
            %====================%

            info = 0; % Default info, indicating that the solver should not exit.

            % Although X should not contain NaN unless there is a bug, we include the following for security.
            % X can be Inf, as finite + finite can be Inf numerically.
            if any(isnan(x) | isinf(x), 'all')
                info = -1;
            end

            % Although NAN_INF_F should not happen unless there is a bug, we include the following for security.
            if (isnan(f) | isinf(f) & f > 0 || isnan(cstrv)) | isinf(cstrv) & cstrv > 0
                info = -2;
            end

            if cstrv <= ctol && f <= ftarget
                info = 1;
            end

            if nf >= maxfun
                info = 3;
            end

            %====================%
            %  Calculation ends  %
            %====================%


        end

    end
end