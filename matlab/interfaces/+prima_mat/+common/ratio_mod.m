classdef ratio_mod
    %--------------------------------------------------------------------------------------------------%
    % This module calculates the reduction ratio for trust-region methods.
    %
    % Coded by Zaikun ZHANG (www.zhangzk.net).
    %
    % Started: September 2021
    %
    % Last Modified: Sunday, December 11, 2022 AM01:29:35
    %--------------------------------------------------------------------------------------------------%

    methods
        function ratio = redrat(~, ared, pred, rshrink)
            %--------------------------------------------------------------------------------------------------%
            % This function evaluates the reduction ratio of a trust-region step, handling Inf/NaN properly.
            %--------------------------------------------------------------------------------------------------%



            % When RATIO <= RSHRINK, DELTA will be shrunk.


            ratio = NaN;


            %====================%
            % Calculation starts %
            %====================%

            if isnan(ared)
                % This should not happen in unconstrained problems due to the moderated extreme barrier.
                ratio = -realmax;
            elseif isnan(pred) || pred <= 0
                % The trust-region subproblem solver fails in this rare case. Instead of terminating as Powell's
                % original code does, we set RATIO as follows so that the solver may continue to progress.
                if ared > 0
                    % The trial point will be accepted, but the trust-region radius will be shrunk if RSHRINK>0.
                    ratio = 0.5 * rshrink;
                else
                    % Set ratio to a large negative number to signify a bad trust-region step, so that the
                    % solver will check whether to take a geometry step or reduce RHO.
                    ratio = -realmax;
                end
            elseif isinf(pred) & pred > 0 & (isinf(ared) & ared > 0)
                ratio = 1.0; % ARED/PRED = NaN if calculated directly.

            elseif isinf(pred) & pred > 0 & (isinf(ared) & ared < 0)
                ratio = -realmax; % ARED/PRED = NaN if calculated directly.

            else
                ratio = ared / pred;
            end

            %====================%
            %  Calculation ends  %
            %====================%



        end

    end
end