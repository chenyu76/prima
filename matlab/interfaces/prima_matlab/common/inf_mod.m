classdef inf_mod
    %--------------------------------------------------------------------------------------------------%
    % This module provides functions that check whether a real number X is infinite or finite.
    % See infnan.f90 for more comments.
    %
    % Coded by Zaikun ZHANG (www.zhangzk.net).
    %
    % Started: July 2020.
    %
    % Last Modified: Tuesday, February 27, 2024 PM10:57:47
    %--------------------------------------------------------------------------------------------------%
    properties
        huge_obj;
    end

    methods
        function obj = inf_mod()
            obj.huge_obj = huge_mod();
        end
        function varargout = is_finite(obj, varargin)
            if numel(varargin) == 1
                [varargout{1:nargout}] = obj.is_finite_sp(varargin{:});
            elseif numel(varargin) == 1
                [varargout{1:nargout}] = obj.is_finite_dp(varargin{:});
            elseif numel(varargin) == 1
                [varargout{1:nargout}] = obj.is_finite_qp(varargin{:});
            end
        end
        function varargout = is_posinf(obj, varargin)
            if numel(varargin) == 1
                [varargout{1:nargout}] = obj.is_posinf_sp(varargin{:});
            elseif numel(varargin) == 1
                [varargout{1:nargout}] = obj.is_posinf_dp(varargin{:});
            elseif numel(varargin) == 1
                [varargout{1:nargout}] = obj.is_posinf_qp(varargin{:});
            end
        end
        function varargout = is_neginf(obj, varargin)
            if numel(varargin) == 1
                [varargout{1:nargout}] = obj.is_neginf_sp(varargin{:});
            elseif numel(varargin) == 1
                [varargout{1:nargout}] = obj.is_neginf_dp(varargin{:});
            elseif numel(varargin) == 1
                [varargout{1:nargout}] = obj.is_neginf_qp(varargin{:});
            end
        end
        function varargout = is_inf(obj, varargin)
            if numel(varargin) == 1
                [varargout{1:nargout}] = obj.is_inf_sp(varargin{:});
            elseif numel(varargin) == 1
                [varargout{1:nargout}] = obj.is_inf_dp(varargin{:});
            elseif numel(varargin) == 1
                [varargout{1:nargout}] = obj.is_inf_qp(varargin{:});
            end
        end
        function y = is_finite_sp(obj, x)
            consts_obj = consts_mod();


            y = (x <= obj.huge_obj.huge_value(x) & x >= -obj.huge_obj.huge_value(x));
        end
        function y = is_finite_dp(obj, x)
            consts_obj = consts_mod();


            y = (x <= obj.huge_obj.huge_value(x) & x >= -obj.huge_obj.huge_value(x));
        end
        function y = is_posinf_sp(obj, x)
            consts_obj = consts_mod();


            y = (abs(x) > obj.huge_obj.huge_value(x)) & (x > 0);
        end
        function y = is_posinf_dp(obj, x)
            consts_obj = consts_mod();


            y = (abs(x) > obj.huge_obj.huge_value(x)) & (x > 0);
        end
        function y = is_neginf_sp(obj, x)
            consts_obj = consts_mod();


            y = (abs(x) > obj.huge_obj.huge_value(x)) & (x < 0);
        end
        function y = is_neginf_dp(obj, x)
            consts_obj = consts_mod();


            y = (abs(x) > obj.huge_obj.huge_value(x)) & (x < 0);
        end
        function y = is_inf_sp(obj, x)
            consts_obj = consts_mod();


            y = (abs(x) > obj.huge_obj.huge_value(x));
        end
        function y = is_inf_dp(obj, x)
            consts_obj = consts_mod();


            y = (abs(x) > obj.huge_obj.huge_value(x));
        end
        function y = is_finite_qp(obj, x)
            consts_obj = consts_mod();


            y = (x <= obj.huge_obj.huge_value(x) & x >= -obj.huge_obj.huge_value(x));
        end
        function y = is_posinf_qp(obj, x)
            consts_obj = consts_mod();


            y = (abs(x) > obj.huge_obj.huge_value(x)) & (x > 0);
        end
        function y = is_neginf_qp(obj, x)
            consts_obj = consts_mod();


            y = (abs(x) > obj.huge_obj.huge_value(x)) & (x < 0);
        end
        function y = is_inf_qp(obj, x)
            consts_obj = consts_mod();


            y = (abs(x) > obj.huge_obj.huge_value(x));
        end

    end
end