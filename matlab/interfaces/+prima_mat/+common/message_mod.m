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
            consts_obj = prima_mat.common.consts_mod();
            debug_obj = prima_mat.common.debug_mod();
            fprint_obj = prima_mat.common.fprint_mod();
            infos_obj = prima_mat.common.infos_mod();
            linalg_obj = prima_mat.common.linalg_mod();
            string_obj = prima_mat.common.string_mod();

            % Compulsory inputs



            % Optional inputs



            % Local variables
            newline_custom = newline;
            srname = "RETMSG";


            funit = NaN; % File storage unit for the writing. Should be an integer of default kind.
            valid_exit_flags = [infos_obj.FTARGET_ACHIEVED, infos_obj.MAXFUN_REACHED, infos_obj.MAXTR_REACHED, infos_obj.SMALL_TR_RADIUS, infos_obj.TRSUBP_FAILED, infos_obj.NAN_INF_F, infos_obj.NAN_INF_X, infos_obj.NAN_INF_MODEL, infos_obj.DAMAGING_ROUNDING, infos_obj.NO_SPACE_BETWEEN_BOUNDS, infos_obj.ZERO_LINEAR_CONSTRAINT];


            % Preconditions
            if consts_obj.DEBUGGING
                debug_obj.assert(any(info == valid_exit_flags, 'all'), "The exit flag is valid", srname);
            end

            %====================%
            % Calculation starts %
            %====================%

            if abs(iprint) < 1
                % No printing
                return
            elseif iprint > 0
                % Print the message to the standard out.
                funit = consts_obj.STDOUT;
                fname = "";
            else                % Print the message to a file named FNAME.
                fname = string_obj.strip(solver) + "_output.txt";
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
                cstrv_loc = linalg_obj.maximum1([consts_obj.ZERO; -constr]); % N.B.: We assume that the constraint is CONSTR >= 0.

            else
                cstrv_loc = consts_obj.ZERO;
            end

            % Decide the return message.
            switch info
            case infos_obj.FTARGET_ACHIEVED
                reason = "the target function value is achieved.";
            case infos_obj.MAXFUN_REACHED
                reason = "the maximal number of function evaluations has been reached.";
            case infos_obj.MAXTR_REACHED
                reason = "the maximal number of trust region iterations has been reached.";
            case infos_obj.SMALL_TR_RADIUS
                reason = "the trust region radius reaches its lower bound.";
            case infos_obj.TRSUBP_FAILED
                reason = "a trust region step has failed to reduce the quadratic model.";
            case infos_obj.NAN_INF_X
                reason = "NaN or Inf occurs in x.";
            case infos_obj.NAN_INF_F
                reason = "the objective or constraint functions return NaN or +Inf.";
            case infos_obj.NAN_INF_MODEL
                reason = "NaN or Inf occurs in the models.";
            case infos_obj.DAMAGING_ROUNDING
                reason = "rounding errors are becoming damaging.";
            case infos_obj.NO_SPACE_BETWEEN_BOUNDS
                reason = "there is no space between the lower and upper bounds of variable.";
            case infos_obj.ZERO_LINEAR_CONSTRAINT
                reason = "one of the linear constraints has a zero gradient";
            case infos_obj.CALLBACK_TERMINATE
                reason = "callback function requested termination of optimization";
            otherwise
                reason = "UNKNOWN EXIT FLAG";
            end
            ret_message = newline_custom + "Return from " + solver + " because " + string_obj.strip(reason);

            if numel(x) <= 2
                x_message = newline_custom + "The corresponding X is: " + string_obj.real2str_vector(x); % Printed in one line

            else
                x_message = newline_custom + "The corresponding X is:" + newline_custom + string_obj.real2str_vector(x);
            end

            if is_constrained
                nf_message = newline_custom + "Number of function values = " + string_obj.int2str(nf) + obj.spaces + "Least value of F = " + string_obj.real2str_scalar(f) + obj.spaces + "Constraint violation = " + string_obj.real2str_scalar(cstrv_loc);
            else
                nf_message = newline_custom + "Number of function values = " + string_obj.int2str(nf) + obj.spaces + "Least value of F = " + string_obj.real2str_scalar(f);
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
            consts_obj = prima_mat.common.consts_mod();
            fprint_obj = prima_mat.common.fprint_mod();
            linalg_obj = prima_mat.common.linalg_mod();
            string_obj = prima_mat.common.string_mod();

            % Compulsory inputs



            % Optional inputs



            % Local variables
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
                funit = consts_obj.STDOUT;
                fname = "";
            else                % Print the message to a file named FNAME.
                fname = string_obj.strip(solver) + "_output.txt";
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
                cstrv_loc = linalg_obj.maximum1([consts_obj.ZERO; -constr]); % N.B.: We assume that the constraint is CONSTR >= 0.

            else
                cstrv_loc = consts_obj.ZERO;
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
                nf_message = newline_custom + "Number of function values = " + string_obj.int2str(nf) + obj.spaces + "Least value of F = " + string_obj.real2str_scalar(f) + obj.spaces + "Constraint violation = " + string_obj.real2str_scalar(cstrv_loc);
            else
                nf_message = newline_custom + "Number of function values = " + string_obj.int2str(nf) + obj.spaces + "Least value of F = " + string_obj.real2str_scalar(f);
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
            consts_obj = prima_mat.common.consts_mod();
            fprint_obj = prima_mat.common.fprint_mod();
            string_obj = prima_mat.common.string_mod();

            % Compulsory inputs



            % Optional inputs


            % Local variables
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
                funit = consts_obj.STDOUT;
                fname = "";
            else                % Print the message to a file named FNAME.
                fname = string_obj.strip(solver) + "_output.txt";
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
            consts_obj = prima_mat.common.consts_mod();
            fprint_obj = prima_mat.common.fprint_mod();
            linalg_obj = prima_mat.common.linalg_mod();
            string_obj = prima_mat.common.string_mod();

            % Compulsory inputs

            % `state` is a string indicating the solver's state when the function evaluation is invoked. Its
            % value can be 'Initialization', 'Trust region', 'Geometry', or 'Rescue'.



            % Optional inputs



            % Local variables
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
                funit = consts_obj.STDOUT;
                fname = "";
            else                % Print the message to a file named FNAME.
                fname = string_obj.strip(solver) + "_output.txt";
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
                cstrv_loc = linalg_obj.maximum1([consts_obj.ZERO; -constr]); % N.B.: We assume that the constraint is CONSTR >= 0.

            else
                cstrv_loc = consts_obj.ZERO;
            end

            delta_message = newline_custom + state + " step with radius = " + string_obj.real2str_scalar(delta);

            if is_constrained
                nf_message = newline_custom + "Function number " + string_obj.int2str(nf) + obj.spaces + "F = " + string_obj.real2str_scalar(f) + obj.spaces + "Constraint violation = " + string_obj.real2str_scalar(cstrv_loc);
            else
                nf_message = newline_custom + "Function number " + string_obj.int2str(nf) + obj.spaces + "F = " + string_obj.real2str_scalar(f);
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