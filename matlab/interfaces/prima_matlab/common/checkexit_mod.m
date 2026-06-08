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

            % Common modules
            consts_obj = consts_mod();
            debug_obj = debug_mod();
            infnan_obj = infnan_mod();
            infos_obj = infos_mod();


            % Inputs






            % Outputs
            info = NaN;

            % Local variables
            srname = "CHECKEXIT_UNC";

            % Preconditions
            if consts_obj.DEBUGGING
                debug_obj.assert(~any([infos_obj.NAN_INF_X, infos_obj.NAN_INF_F, infos_obj.FTARGET_ACHIEVED, infos_obj.MAXFUN_REACHED] == infos_obj.INFO_DFT, 'all'), "NAN_INF_X, NAN_INF_F, FTARGET_ACHIEVED, and MAXFUN_REACHED differ from INFO_DFT", srname);
                % X does not contain NaN if the initial X does not contain NaN and the subroutines generating
                % trust-region/geometry steps work properly so that they never produce a step containing NaN/Inf.
                debug_obj.assert(~any(infnan_obj.is_nan(x), 'all'), "X does not contain NaN", srname);
                % With the moderated extreme barrier, F cannot be NaN/+Inf.
                debug_obj.assert(~(infnan_obj.is_nan_sp(f) || infnan_obj.is_posinf(f)), "F is not NaN/+Inf", srname);
            end

            %====================%
            % Calculation starts %
            %====================%

            info = infos_obj.INFO_DFT; % Default info, indicating that the solver should not exit.

            % Although X should not contain NaN unless there is a bug, we include the following for security.
            % X can be Inf, as finite + finite can be Inf numerically.
            if any(infnan_obj.is_nan(x) | infnan_obj.is_inf(x), 'all')
                info = infos_obj.NAN_INF_X;
            end

            % Although NAN_INF_F should not happen unless there is a bug, we include the following for security.
            if infnan_obj.is_nan_sp(f) || infnan_obj.is_posinf(f)
                info = infos_obj.NAN_INF_F;
            end

            if f <= ftarget
                info = infos_obj.FTARGET_ACHIEVED;
            end

            if nf >= maxfun
                info = infos_obj.MAXFUN_REACHED;
            end

            %====================%
            %  Calculation ends  %
            %====================%

            % Postconditions
            if consts_obj.DEBUGGING
                debug_obj.assert(any([infos_obj.INFO_DFT, infos_obj.NAN_INF_X, infos_obj.FTARGET_ACHIEVED, infos_obj.MAXFUN_REACHED] == info, 'all'), "INFO is NAN_INF_X, FTARGET_ACHIEVED, MAXFUN_REACHED, or INFO_DFT", srname);
            end

        end
        function info = checkexit_con(~, maxfun, nf, cstrv, ctol, f, ftarget, x)
            %--------------------------------------------------------------------------------------------------%
            % This module checks whether to exit the solver in the constrained case.
            %--------------------------------------------------------------------------------------------------%

            % Common modules
            consts_obj = consts_mod();
            debug_obj = debug_mod();
            infnan_obj = infnan_mod();
            infos_obj = infos_mod();


            % Inputs








            % Outputs
            info = NaN;

            % Local variables
            srname = "CHECKEXIT_CON";

            % Preconditions
            if consts_obj.DEBUGGING
                debug_obj.assert(~any([infos_obj.NAN_INF_X, infos_obj.NAN_INF_F, infos_obj.FTARGET_ACHIEVED, infos_obj.MAXFUN_REACHED] == infos_obj.INFO_DFT, 'all'), "NAN_INF_X, NAN_INF_F, FTARGET_ACHIEVED, and MAXFUN_REACHED differ from INFO_DFT", srname);
                % X does not contain NaN if the initial X does not contain NaN and the subroutines generating
                % trust-region/geometry steps work properly so that they never produce a step containing NaN/Inf.
                debug_obj.assert(~any(infnan_obj.is_nan(x), 'all'), "X does not contain NaN", srname);
                % With the moderated extreme barrier, F or CSTRV cannot be NaN/+Inf.
                debug_obj.assert(~(infnan_obj.is_nan_sp(f) || infnan_obj.is_posinf(f) || infnan_obj.is_nan_sp(cstrv) || infnan_obj.is_posinf(cstrv)), "F or CSTRV is not NaN/+Inf", srname);
            end

            %====================%
            % Calculation starts %
            %====================%

            info = infos_obj.INFO_DFT; % Default info, indicating that the solver should not exit.

            % Although X should not contain NaN unless there is a bug, we include the following for security.
            % X can be Inf, as finite + finite can be Inf numerically.
            if any(infnan_obj.is_nan(x) | infnan_obj.is_inf(x), 'all')
                info = infos_obj.NAN_INF_X;
            end

            % Although NAN_INF_F should not happen unless there is a bug, we include the following for security.
            if infnan_obj.is_nan_sp(f) || infnan_obj.is_posinf(f) || infnan_obj.is_nan_sp(cstrv) || infnan_obj.is_posinf(cstrv)
                info = infos_obj.NAN_INF_F;
            end

            if cstrv <= ctol && f <= ftarget
                info = infos_obj.FTARGET_ACHIEVED;
            end

            if nf >= maxfun
                info = infos_obj.MAXFUN_REACHED;
            end

            %====================%
            %  Calculation ends  %
            %====================%

            % Postconditions
            if consts_obj.DEBUGGING
                debug_obj.assert(any([infos_obj.INFO_DFT, infos_obj.NAN_INF_F, infos_obj.FTARGET_ACHIEVED, infos_obj.MAXFUN_REACHED] == info, 'all'), "INFO is NAN_INF_X, FTARGET_ACHIEVED, MAXFUN_REACHED, or INFO_DFT", srname);
            end

        end

    end
end