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
            consts_obj = prima_mat.common.consts_mod();
            infnan_obj = prima_mat.common.infnan_mod();
            linalg_obj = prima_mat.common.linalg_mod();

            % Inputs

            % Outputs
            y = NaN(numel(x), 1);

            y(:) = x;
            y(linalg_obj.trueloc(infnan_obj.is_nan_sp(x))) = consts_obj.ZERO;
            y(:) = max(-consts_obj.REALMAX, min(consts_obj.REALMAX, y));
        end
        function y = moderatef(~, f)
            %--------------------------------------------------------------------------------------------------%
            % This function moderates the function value of a MINIMIZATION problem. It replaces NaN and any
            % value above FUNCMAX by FUNCMAX.
            %--------------------------------------------------------------------------------------------------%
            consts_obj = prima_mat.common.consts_mod();
            infnan_obj = prima_mat.common.infnan_mod();

            % Inputs

            % Outputs


            y = f;
            if infnan_obj.is_nan_sp(y)
                y = consts_obj.FUNCMAX;
            end
            y = max(-consts_obj.REALMAX, min(consts_obj.FUNCMAX, y));
            % We may moderate huge negative function values as follows, but we decide not to.
            %y = max(-FUNCMAX, min(FUNCMAX, y))
        end
        function y = moderatec(~, c)
            %--------------------------------------------------------------------------------------------------%
            % This function moderates the constraint value, the constraint demanding this value to be NONNEGATIVE.
            % It replaces any value below -CONSTRMAX by -CONSTRMAX, and any NaN or value above CONSTRMAX by
            % CONSTRMAX.
            %--------------------------------------------------------------------------------------------------%
            consts_obj = prima_mat.common.consts_mod();
            infnan_obj = prima_mat.common.infnan_mod();
            linalg_obj = prima_mat.common.linalg_mod();

            % Inputs

            % Outputs
            y = NaN(numel(c), 1);

            y(:) = c;
            y(linalg_obj.trueloc(infnan_obj.is_nan_sp(c))) = consts_obj.CONSTRMAX;
            y(:) = max(-consts_obj.CONSTRMAX, min(consts_obj.CONSTRMAX, y));
        end
        function f = evaluatef(obj, calfun, x)
            %--------------------------------------------------------------------------------------------------%
            % This function evaluates CALFUN at X, setting F to the objective function value. Nan/Inf are
            % handled by a moderated extreme barrier.
            %--------------------------------------------------------------------------------------------------%
            % Common modules
            consts_obj = prima_mat.common.consts_mod();
            debug_obj = prima_mat.common.debug_mod();
            infnan_obj = prima_mat.common.infnan_mod();


            % Inputs
            % N.B.: INTENT cannot be specified if a dummy procedure is not a POINTER


            % Output


            % Local variables
            srname = "EVALUATEF";

            % Preconditions
            if consts_obj.DEBUGGING
                % X should not contain NaN if the initial X does not contain NaN and the subroutines generating
                % trust-region/geometry steps work properly so that they never produce a step containing NaN/Inf.
                debug_obj.assert(~any(infnan_obj.is_nan_sp(x), 'all'), "X does not contain NaN", srname);
            end

            %====================%
            % Calculation starts %
            %====================%

            if any(infnan_obj.is_nan_sp(x), 'all')
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

            % Postconditions
            if consts_obj.DEBUGGING
                % With X not containing NaN, and with the moderated extreme barrier, F cannot be NaN/+Inf.
                debug_obj.assert(~(infnan_obj.is_nan_sp(f) || infnan_obj.is_posinf(f)), "F is not NaN/+Inf", srname);
            end

        end
        function [f, constr] = evaluatefc(obj, calcfc, x, constr)
            %--------------------------------------------------------------------------------------------------%
            % This function evaluates CALCFC at X, setting F to the objective function value and CONSTR to the
            % constraint value. Nan/Inf are handled by a moderated extreme barrier.
            %--------------------------------------------------------------------------------------------------%
            % Common modules
            consts_obj = prima_mat.common.consts_mod();
            debug_obj = prima_mat.common.debug_mod();
            infnan_obj = prima_mat.common.infnan_mod();


            % Inputs
            % N.B.: INTENT cannot be specified if a dummy procedure is not a POINTER


            % Outputs



            % Local variables
            srname = "EVALUATEFC";

            % Preconditions
            if consts_obj.DEBUGGING
                % X should not contain NaN if the initial X does not contain NaN and the subroutines generating
                % trust-region/geometry steps work properly so that they never produce a step containing NaN/Inf.
                debug_obj.assert(~any(infnan_obj.is_nan_sp(x), 'all'), "X does not contain NaN", srname);
            end

            %====================%
            % Calculation starts %
            %====================%

            if any(infnan_obj.is_nan_sp(x), 'all')
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

            % Postconditions
            if consts_obj.DEBUGGING
                % With X not containing NaN, and with the moderated extreme barrier, F cannot be NaN/+Inf, and
                % CONSTR cannot be NaN/+Inf.
                debug_obj.assert(~(infnan_obj.is_nan_sp(f) || infnan_obj.is_posinf(f)), "F is not NaN/+Inf", srname);
                debug_obj.assert(~any(infnan_obj.is_nan_sp(constr) | infnan_obj.is_posinf(constr), 'all'), "CONSTR does not contain NaN/+Inf", srname);
            end

        end

    end
end