classdef huge_mod
    %--------------------------------------------------------------------------------------------------%
    % This module provides a function that returns HUGE(X). See infnan.f90 for more comments.
    %
    % Coded by Zaikun ZHANG (www.zhangzk.net).
    %
    % Started: July 2020.
    %
    % Last Modified: Tuesday, February 27, 2024 PM11:02:29
    %--------------------------------------------------------------------------------------------------%

    methods
        function varargout = huge_value(obj, varargin)
            if numel(varargin) == 1
                [varargout{1:nargout}] = obj.huge_value_sp(varargin{:});
            elseif numel(varargin) == 1
                [varargout{1:nargout}] = obj.huge_value_dp(varargin{:});
            else
                [varargout{1:nargout}] = obj.huge_value_qp(varargin{:});
            end
        end
        function y = huge_value_sp(~, x)
            consts_obj = prima_mat.common.consts_mod();


            y = realmax;
        end
        function y = huge_value_dp(~, x)
            consts_obj = prima_mat.common.consts_mod();


            y = realmax;
        end
        function y = huge_value_qp(~, x)
            consts_obj = prima_mat.common.consts_mod();


            y = realmax;
        end

    end
end