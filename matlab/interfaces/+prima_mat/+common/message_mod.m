classdef message_mod
    %--------------------------------------------------------------------------------------------------%
    % This module provides some subroutines that print messages to terminal/files.
    %
    % N.B.:
    % 1. In case parallelism is desirable (especially during initialization), the subroutines may
    % have to be modified or disabled due to the IO operations.
    % 2. IPRINT indicates the level of verbosity, which increases with the absolute value of IPRINT.
    % IPRINT = +/-3 can be expensive due to high IO operations.
    %
    % Coded by Zaikun ZHANG (www.zhangzk.net) based on Powell's code and papers.
    %
    % Started: July 2020
    %
    % Last Modified: Sunday, March 31, 2024 PM04:55:58
    %--------------------------------------------------------------------------------------------------%
    properties (Access = private)
        spaces;
    end

    methods
        function obj = message_mod()
            obj.spaces = "   ";
        end
        function retmsg(obj, solver, info, iprint, nf, f, x, varargin)
            %--------------------------------------------------------------------------------------------------%
            % This subroutine prints messages at return.
            %--------------------------------------------------------------------------------------------------%


            fprint_obj = prima_mat.common.fprint_mod();


            string_obj = prima_mat.common.string_mod();

            % Compulsory inputs



            % Optional inputs



            newline_custom = newline;


            funit = NaN; % File storage unit for the writing. Should be an integer of default kind.



            %====================%
            % Calculation starts %
            %====================%

            if abs(iprint) < 1
                % No printing
                return
            elseif iprint > 0
                % Print the message to the standard out.
                funit = 1;
                fname = "";
            else                % Print the message to a file named FNAME.
                fname = strip(solver) + "_output.txt";
            end

            % Decide whether the problem is truly constrained.
            ipObj = inputParser();
            addParameter(ipObj, 'cstrv', NaN);
            addParameter(ipObj, 'constr', NaN);
            parse(ipObj, varargin{:});
            cstrv = ipObj.Results.cstrv;
            constr = ipObj.Results.constr;
            if ismember('constr', ipObj.UsingDefaults)
                is_constrained = ~ismember('cstrv', ipObj.UsingDefaults);
            else
                is_constrained = (numel(constr) > 0);
            end

            % Decide the constraint violation.
            if ~ismember('cstrv', ipObj.UsingDefaults)
                cstrv_loc = cstrv;
            elseif ~ismember('constr', ipObj.UsingDefaults)
                cstrv_loc = max([0.0; -constr], [], 'all'); % N.B.: We assume that the constraint is CONSTR >= 0.

            else
                cstrv_loc = 0.0;
            end

            % Decide the return message.
            switch info
            case 1
                reason = "the target function value is achieved.";
            case 3
                reason = "the maximal number of function evaluations has been reached.";
            case 20
                reason = "the maximal number of trust region iterations has been reached.";
            case 0
                reason = "the trust region radius reaches its lower bound.";
            case 2
                reason = "a trust region step has failed to reduce the quadratic model.";
            case -1
                reason = "NaN or Inf occurs in x.";
            case -2
                reason = "the objective or constraint functions return NaN or +Inf.";
            case -3
                reason = "NaN or Inf occurs in the models.";
            case 7
                reason = "rounding errors are becoming damaging.";
            case 6
                reason = "there is no space between the lower and upper bounds of variable.";
            case 8
                reason = "one of the linear constraints has a zero gradient";
            case 30
                reason = "callback function requested termination of optimization";
            otherwise
                reason = "UNKNOWN EXIT FLAG";
            end
            ret_message = newline_custom + "Return from " + solver + " because " + strip(reason);

            if numel(x) <= 2
                x_message = newline_custom + "The corresponding X is: " + string_obj.real2str_vector(x); % Printed in one line

            else
                x_message = newline_custom + "The corresponding X is:" + newline_custom + string_obj.real2str_vector(x);
            end

            if is_constrained
                nf_message = newline_custom + "Number of function values = " + int2str(nf) + obj.spaces + "Least value of F = " + string_obj.real2str_scalar(f) + obj.spaces + "Constraint violation = " + string_obj.real2str_scalar(cstrv_loc);
            else
                nf_message = newline_custom + "Number of function values = " + int2str(nf) + obj.spaces + "Least value of F = " + string_obj.real2str_scalar(f);
            end

            if is_constrained && ~ismember('constr', ipObj.UsingDefaults)
                if numel(constr) <= 2
                    constr_message = newline_custom + "The constraint value is: " + string_obj.real2str_vector(constr); % Printed in one line

                else
                    constr_message = newline_custom + "The constraint value is:" + newline_custom + string_obj.real2str_vector(constr);
                end
            else
                constr_message = "";
            end

            % Print the message.
            if abs(iprint) >= 2
                message = newline_custom + ret_message + nf_message + x_message + constr_message + newline_custom;
            else
                message = ret_message + nf_message + x_message + constr_message + newline_custom;
            end
            if strlength(fname) > 0
                fprint_obj.fprint(message, 'fname', fname, 'faction', "append");
            else
                fprint_obj.fprint(message, 'funit', funit, 'faction', "append");
            end

            %====================%
            %  Calculation ends  %
            %====================%
        end
        function rhomsg(obj, solver, iprint, nf, delta, f, rho, x, varargin)
            %--------------------------------------------------------------------------------------------------%
            % This subroutine prints messages when RHO is updated.
            %--------------------------------------------------------------------------------------------------%

            fprint_obj = prima_mat.common.fprint_mod();

            string_obj = prima_mat.common.string_mod();

            % Compulsory inputs



            % Optional inputs



            newline_custom = newline;


            funit = NaN; % Logical unit for the writing. Should be an integer of default kind.



            %====================%
            % Calculation starts %
            %====================%

            if abs(iprint) < 2
                % No printing
                return
            elseif iprint > 0
                % Print the message to the standard out.
                funit = 1;
                fname = "";
            else                % Print the message to a file named FNAME.
                fname = strip(solver) + "_output.txt";
            end

            % Decide whether the problem is truly constrained.
            ipObj = inputParser();
            addParameter(ipObj, 'cstrv', NaN);
            addParameter(ipObj, 'constr', NaN);
            addParameter(ipObj, 'cpen', NaN);
            parse(ipObj, varargin{:});
            cstrv = ipObj.Results.cstrv;
            constr = ipObj.Results.constr;
            cpen = ipObj.Results.cpen;
            if ismember('constr', ipObj.UsingDefaults)
                is_constrained = ~ismember('cstrv', ipObj.UsingDefaults);
            else
                is_constrained = (numel(constr) > 0);
            end

            % Decide the constraint violation.
            if ~ismember('cstrv', ipObj.UsingDefaults)
                cstrv_loc = cstrv;
            elseif ~ismember('constr', ipObj.UsingDefaults)
                cstrv_loc = max([0.0; -constr], [], 'all'); % N.B.: We assume that the constraint is CONSTR >= 0.

            else
                cstrv_loc = 0.0;
            end

            if ismember('cpen', ipObj.UsingDefaults)
                rho_message = newline_custom + "New RHO = " + string_obj.real2str_scalar(rho) + obj.spaces + "Delta = " + string_obj.real2str_scalar(delta);
            else
                rho_message = newline_custom + "New RHO = " + string_obj.real2str_scalar(rho) + obj.spaces + "Delta = " + string_obj.real2str_scalar(delta) + obj.spaces + "CPEN = " + string_obj.real2str_scalar(cpen);
            end

            if numel(x) <= 2
                x_message = newline_custom + "The corresponding X is: " + string_obj.real2str_vector(x); % Printed in one line

            else
                x_message = newline_custom + "The corresponding X is:" + newline_custom + string_obj.real2str_vector(x);
            end

            if is_constrained
                nf_message = newline_custom + "Number of function values = " + int2str(nf) + obj.spaces + "Least value of F = " + string_obj.real2str_scalar(f) + obj.spaces + "Constraint violation = " + string_obj.real2str_scalar(cstrv_loc);
            else
                nf_message = newline_custom + "Number of function values = " + int2str(nf) + obj.spaces + "Least value of F = " + string_obj.real2str_scalar(f);
            end

            if is_constrained && ~ismember('constr', ipObj.UsingDefaults)
                if numel(constr) <= 2
                    constr_message = newline_custom + "The constraint value is: " + string_obj.real2str_vector(constr); % Printed in one line

                else
                    constr_message = newline_custom + "The constraint value is:" + newline_custom + string_obj.real2str_vector(constr);
                end
            else
                constr_message = "";
            end

            % Print the message.
            if abs(iprint) >= 3
                message = newline_custom + rho_message + nf_message + x_message + constr_message;
            else
                message = rho_message + nf_message + x_message + constr_message;
            end
            if strlength(fname) > 0
                fprint_obj.fprint(message, 'fname', fname, 'faction', "append");
            else
                fprint_obj.fprint(message, 'funit', funit, 'faction', "append");
            end

            %====================%
            %  Calculation ends  %
            %====================%
        end
        function cpenmsg(~, solver, iprint, varargin)
            %--------------------------------------------------------------------------------------------------%
            % This subroutine prints a message when CPEN is updated.
            %--------------------------------------------------------------------------------------------------%

            fprint_obj = prima_mat.common.fprint_mod();
            string_obj = prima_mat.common.string_mod();

            % Compulsory inputs



            % Optional inputs



            newline_custom = newline;


            funit = NaN; % Logical unit for the writing. Should be an integer of default kind.

            %====================%
            % Calculation starts %
            %====================%

            if abs(iprint) < 2
                % No printing
                return
            elseif iprint > 0
                % Print the message to the standard out.
                funit = 1;
                fname = "";
            else                % Print the message to a file named FNAME.
                fname = strip(solver) + "_output.txt";
            end

            % Print the message.
            ipObj = inputParser();
            addParameter(ipObj, 'cpen', NaN);
            parse(ipObj, varargin{:});
            cpen = ipObj.Results.cpen;
            if abs(iprint) >= 3
                message = newline_custom + "Set CPEN to " + string_obj.real2str_scalar(cpen);
            else
                message = newline_custom + newline_custom + "Set CPEN to " + string_obj.real2str_scalar(cpen);
            end
            if strlength(fname) > 0
                fprint_obj.fprint(message, 'fname', fname, 'faction', "append");
            else
                fprint_obj.fprint(message, 'funit', funit, 'faction', "append");
            end

            %====================%
            %  Calculation ends  %
            %====================%
        end
        function fmsg(obj, solver, state, iprint, nf, delta, f, x, varargin)
            %--------------------------------------------------------------------------------------------------%
            % This subroutine prints messages for each evaluation of the objective function.
            %--------------------------------------------------------------------------------------------------%

            fprint_obj = prima_mat.common.fprint_mod();

            string_obj = prima_mat.common.string_mod();

            % Compulsory inputs

            % `state` is a string indicating the solver's state when the function evaluation is invoked. Its
            % value can be 'Initialization', 'Trust region', 'Geometry', or 'Rescue'.



            % Optional inputs



            newline_custom = newline;


            funit = NaN; % Logical unit for the writing. Should be an integer of default kind.



            %====================%
            % Calculation starts %
            %====================%

            if abs(iprint) < 3
                % No printing
                return
            elseif iprint > 0
                % Print the message to the standard out.
                funit = 1;
                fname = "";
            else                % Print the message to a file named FNAME.
                fname = strip(solver) + "_output.txt";
            end

            % Decide whether the problem is truly constrained.
            ipObj = inputParser();
            addParameter(ipObj, 'cstrv', NaN);
            addParameter(ipObj, 'constr', NaN);
            parse(ipObj, varargin{:});
            cstrv = ipObj.Results.cstrv;
            constr = ipObj.Results.constr;
            if ismember('constr', ipObj.UsingDefaults)
                is_constrained = ~ismember('cstrv', ipObj.UsingDefaults);
            else
                is_constrained = (numel(constr) > 0);
            end

            % Decide the constraint violation.
            if ~ismember('cstrv', ipObj.UsingDefaults)
                cstrv_loc = cstrv;
            elseif ~ismember('constr', ipObj.UsingDefaults)
                cstrv_loc = max([0.0; -constr], [], 'all'); % N.B.: We assume that the constraint is CONSTR >= 0.

            else
                cstrv_loc = 0.0;
            end

            delta_message = newline_custom + state + " step with radius = " + string_obj.real2str_scalar(delta);

            if is_constrained
                nf_message = newline_custom + "Function number " + int2str(nf) + obj.spaces + "F = " + string_obj.real2str_scalar(f) + obj.spaces + "Constraint violation = " + string_obj.real2str_scalar(cstrv_loc);
            else
                nf_message = newline_custom + "Function number " + int2str(nf) + obj.spaces + "F = " + string_obj.real2str_scalar(f);
            end

            if numel(x) <= 2
                x_message = newline_custom + "The corresponding X is: " + string_obj.real2str_vector(x); % Printed in one line

            else
                x_message = newline_custom + "The corresponding X is:" + newline_custom + string_obj.real2str_vector(x);
            end

            if is_constrained && ~ismember('constr', ipObj.UsingDefaults)
                if numel(constr) <= 2
                    constr_message = newline_custom + "The constraint value is: " + string_obj.real2str_vector(constr); % Printed in one line

                else
                    constr_message = newline_custom + "The constraint value is:" + newline_custom + string_obj.real2str_vector(constr);
                end
            else
                constr_message = "";
            end

            % Print the message.
            message = delta_message + nf_message + x_message + constr_message;
            if strlength(fname) > 0
                fprint_obj.fprint(message, 'fname', fname, 'faction', "append");
            else
                fprint_obj.fprint(message, 'funit', funit, 'faction', "append");
            end

            %====================%
            %  Calculation ends  %
            %====================%
        end

    end
end