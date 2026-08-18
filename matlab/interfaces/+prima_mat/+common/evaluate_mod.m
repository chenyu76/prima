classdef evaluate_mod
    %--------------------------------------------------------------------------------------------------%
    % This is a module evaluating the objective/constraint function with Nan/Inf handling.
    %
    % Coded by Zaikun ZHANG (www.zhangzk.net).
    %
    % Started: August 2021
    %
    % Last Modified: Monday, September 25, 2023 PM08:52:04
    %--------------------------------------------------------------------------------------------------%

    methods
        function varargout = evaluate(obj, varargin)
            if numel(varargin) == 3
                [varargout{1:nargout}] = obj.evaluatef(varargin{:});
            else
                [varargout{1:nargout}] = obj.evaluatefc(varargin{:});
            end
        end
        function y = moderatex(~, x)
            %--------------------------------------------------------------------------------------------------%
            % This function moderates a decision variable. It replaces NaN by 0 and Inf/-Inf by REALMAX/-REALMAX.
            %--------------------------------------------------------------------------------------------------%


            y = NaN(numel(x), 1);

            y(:) = x;
            y(isnan(x)) = 0.0;
            y = max(-realmax, min(realmax, y));
        end
        function y = moderatef(~, f)
            %--------------------------------------------------------------------------------------------------%
            % This function moderates the function value of a MINIMIZATION problem. It replaces NaN and any
            % value above FUNCMAX by FUNCMAX.
            %--------------------------------------------------------------------------------------------------%


            y = f;
            if isnan(y)
                y = 1.0e30;
            end
            y = max(-realmax, min(1.0e30, y));
            % We may moderate huge negative function values as follows, but we decide not to.
            %y = max(-FUNCMAX, min(FUNCMAX, y))
        end
        function y = moderatec(~, c)
            %--------------------------------------------------------------------------------------------------%
            % This function moderates the constraint value, the constraint demanding this value to be NONNEGATIVE.
            % It replaces any value below -CONSTRMAX by -CONSTRMAX, and any NaN or value above CONSTRMAX by
            % CONSTRMAX.
            %--------------------------------------------------------------------------------------------------%


            y = NaN(numel(c), 1);

            y(:) = c;
            y(isnan(c)) = 1.0e30;
            y = max(-1.0e30, min(1.0e30, y));
        end
        function f = evaluatef(obj, calfun, x)
            %--------------------------------------------------------------------------------------------------%
            % This function evaluates CALFUN at X, setting F to the objective function value. Nan/Inf are
            % handled by a moderated extreme barrier.
            %--------------------------------------------------------------------------------------------------%


            %====================%
            % Calculation starts %
            %====================%

            if any(isnan(x), 'all')
                % Although this should not happen unless there is a bug, we include this case for robustness.
                f = sum(x, 'all'); % Set F to NaN

            else
                f = calfun(obj.moderatex(x)); % Evaluate F; We moderate X before doing so.

                % Moderated extreme barrier: replace NaN/huge objective or constraint values with a large but
                % finite value. This is naive. Better approaches surely exist.
                f = obj.moderatef(f);
            end

            %====================%
            %  Calculation ends  %
            %====================%


        end
        function [f, constr] = evaluatefc(obj, calcfc, x, constr)
            %--------------------------------------------------------------------------------------------------%
            % This function evaluates CALCFC at X, setting F to the objective function value and CONSTR to the
            % constraint value. Nan/Inf are handled by a moderated extreme barrier.
            %--------------------------------------------------------------------------------------------------%


            %====================%
            % Calculation starts %
            %====================%

            if any(isnan(x), 'all')
                % Although this should not happen unless there is a bug, we include this case for robustness.
                % Set F, CONSTR, and CSTRV to NaN.
                f = sum(x, 'all');
                constr(:) = f;
            else
                [f, constr] = calcfc(obj.moderatex(x), constr); % Evaluate F and CONSTR; We moderate X before doing so.

                % Moderated extreme barrier: replace NaN/huge objective or constraint values with a large but
                % finite value. This is naive, and better approaches surely exist.
                f = obj.moderatef(f);
                constr(:) = obj.moderatec(constr);
            end

            %====================%
            %  Calculation ends  %
            %====================%


        end

    end
end