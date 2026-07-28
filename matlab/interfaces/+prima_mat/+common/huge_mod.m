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
                [varargout{1:nargout}] = obj.huge_value_sp();
            elseif numel(varargin) == 1
                [varargout{1:nargout}] = obj.huge_value_dp();
            else
                [varargout{1:nargout}] = obj.huge_value_qp();
            end
        end
        function y = huge_value_sp(~)


            y = realmax;
        end
        function y = huge_value_dp(~)


            y = realmax;
        end
        function y = huge_value_qp(~)


            y = realmax;
        end

    end
end