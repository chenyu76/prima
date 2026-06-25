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
            consts_obj = consts_mod();
            infnan_obj = infnan_mod();
            debug_obj = debug_mod();

            % Inputs


            % When RATIO <= RSHRINK, DELTA will be shrunk.

            % Outputs
            ratio = NaN;

            % Local variables
            srname = "REDRAT";

            % Preconditions
            if consts_obj.DEBUGGING
                debug_obj.assert(rshrink >= 0, "RSHRINK >= 0", srname);
            end

            %====================%
            % Calculation starts %
            %====================%

            if infnan_obj.is_nan_sp(ared)
                % This should not happen in unconstrained problems due to the moderated extreme barrier.
                ratio = -consts_obj.REALMAX;
            elseif infnan_obj.is_nan_sp(pred) || pred <= 0
                % The trust-region subproblem solver fails in this rare case. Instead of terminating as Powell's
                % original code does, we set RATIO as follows so that the solver may continue to progress.
                if ared > 0
                    % The trial point will be accepted, but the trust-region radius will be shrunk if RSHRINK>0.
                    ratio = consts_obj.HALF * rshrink;
                else
                    % Set ratio to a large negative number to signify a bad trust-region step, so that the
                    % solver will check whether to take a geometry step or reduce RHO.
                    ratio = -consts_obj.REALMAX;
                end
            elseif infnan_obj.is_posinf(pred) && infnan_obj.is_posinf(ared)
                ratio = consts_obj.ONE; % ARED/PRED = NaN if calculated directly.

            elseif infnan_obj.is_posinf(pred) && infnan_obj.is_neginf(ared)
                ratio = -consts_obj.REALMAX; % ARED/PRED = NaN if calculated directly.

            else
                ratio = ared / pred;
            end

            %====================%
            %  Calculation ends  %
            %====================%

            % Postconditions
            if consts_obj.DEBUGGING
                debug_obj.assert(~infnan_obj.is_nan_sp(ratio), "RATIO is not NaN", srname);
            end

        end

    end
end