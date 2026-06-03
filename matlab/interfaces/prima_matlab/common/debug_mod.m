classdef debug_mod
    %--------------------------------------------------------------------------------------------------%
    % This is a module defining some procedures concerning debugging, errors, and warnings.
    %
    % Coded by Zaikun ZHANG (www.zhangzk.net).
    %
    % Started: July 2020.
    %
    % Last Modified: Tue 09 Sep 2025 11:53:28 PM CST
    %--------------------------------------------------------------------------------------------------%

    methods
        function assert(obj, condition, description, srname)
            %--------------------------------------------------------------------------------------------------%
            % This subroutine checks whether ASSERTION is true.
            % If no but DEBUGGING is true, print the following message to STDERR and then stop the program:
            % 'ERROR: ' // SRNAME // 'Assertion fails: ' // DESCRIPTION
            % MATLAB analogue: assert(condition, sprintf('%s: Assertion fails: %s', srname, description))
            % Python analogue: assert condition, srname + ': Assertion fails: ' + description
            % C analogue: assert(condition)  
            % N.B.: As in C, we design ASSERT to operate only in the debug mode, i.e., when 0 == 1;
            % when 0 == 0, ASSERT does nothing. For the checking that should take effect in both
            % the debug and release modes, use VALIDATE (see below) instead. In the optimized mode of Python
            % (python -O), the Python `assert` will also be ignored. MATLAB does not behave in this way.
            %--------------------------------------------------------------------------------------------------%
            consts_obj = consts_mod();
            infos_obj = infos_mod();
            % A condition that is expected to be true
            % Description of the condition in human language
            % Name of the subroutine that calls this procedure
            if consts_obj.DEBUGGING && ~condition
                obj.errstop(strtrim(strjust(srname, 'left')), "Assertion fails: " + strtrim(strjust(description, 'left')), 'code', infos_obj.ASSERTION_FAILS);
            end
        end
        function validate(obj, condition, description, srname)
            %--------------------------------------------------------------------------------------------------%
            % This subroutine checks whether CONDITION is true.
            % If no, print the following message to STDERR and then stop the program:
            % 'ERROR: ' // SRNAME // 'Validation fails: ' // DESCRIPTION
            % MATLAB analogue: assert(condition, sprintf('%s: Validation fails: %s', srname, description))
            % In Python or C, VALIDATE can be implemented following the Fortran implementation below.
            % N.B.: ASSERT checks the condition only when debugging, but VALIDATE does it always.
            %--------------------------------------------------------------------------------------------------%
            infos_obj = infos_mod();
            % A condition that is expected to be true
            % Description of the condition in human language
            % Name of the subroutine that calls this procedure
            if ~condition
                obj.errstop(strtrim(strjust(srname, 'left')), "Validation fails: " + strtrim(strjust(description, 'left')), 'code', infos_obj.VALIDATION_FAILS);
            end
        end
        function wassert(obj, condition, description, srname)
            %--------------------------------------------------------------------------------------------------%
            % This subroutine checks whether CONDITION is true.
            % If no but DEBUGGING is true, print the following message to STDERR (but do not stop the program):
            % 'Warning: ' // SRNAME // 'Assertion fails: ' // DESCRIPTION
            % MATLAB analogue:
            % %if ~condition
            % %    warning(sprintf('%s: Assertion fails: %s', srname, description))
            % %end
            % In Python or C, WASSERT can be implemented following the Fortran implementation below.
            % N.B.: When DEBUGGING is true, ASSERT stops the program with an error if the condition is false,
            % but WASSERT only raises a warning.
            %--------------------------------------------------------------------------------------------------%
            consts_obj = consts_mod();
            % A condition that is expected to be true
            % Description of the condition in human language
            % Name of the subroutine that calls this procedure
            if consts_obj.DEBUGGING && ~condition
                obj.backtr();
                obj.warning(strtrim(strjust(srname, 'left')), "Assertion fails: " + strtrim(strjust(description, 'left')));
            end
        end
        function errstop(obj, srname, msg, varargin)
            %--------------------------------------------------------------------------------------------------%
            % This subroutine prints 'ERROR: '//STRIP(SRNAME)//': '//STRIP(MSG)//'.' to STDERR, then stop.
            % It also calls BACKTR to print the backtrace.
            %--------------------------------------------------------------------------------------------------%
            consts_obj = consts_mod();




            % `backtr` prints a backtrace. With gfortran 12, even without calling `backtrace`, a backtrace is
            % printed when the program is stopped by an error stop.
            obj.backtr();

            fprintf(consts_obj.STDERR, '\n%s\n\n', "ERROR: " + strtrim(strjust(srname, 'left')) + ": " + strtrim(strjust(msg, 'left')) + ".");
            ipObj = inputParser();
            addParameter(ipObj, 'code', NaN);
            parse(ipObj, varargin{:});
            code = ipObj.Results.code;
            if ismember('code', ipObj.UsingDefaults)
                %Unsupported Statement: StmtErrorStop Nothing

            else
                % N.B.: In Fortran 2008, stop code must be a scalar default character or integer CONSTANT
                % expression, but Fortran 2018 lifts the requirement on constancy. gfortran is strict in this
                % aspect. Consequently, for gfortran, compile with either `-std=f2018` or no `-std` at all.
                %Unsupported Statement: StmtErrorStop (Just (Variable "code"))

            end
            % N.B.
            % 1. ERROR STOP means to stop the whole program.
            % 2. (Zaikun 20230410): We prefer ERROR STOP to STOP, as the former has been allowed in PURE
            % procedures since F2018. Later, when F2018 is better supported, we should take advantage of this
            % feature to make our subroutines PURE whenever possible.
        end
        function backtr(~)
            %--------------------------------------------------------------------------------------------------%
            % This subroutine calls a compiler-dependent intrinsic to show a backtrace if we are in the
            % debugging mode, i.e., 0 == 1.
            % N.B.:
            % 1. The intrinsic is compiler-dependent and does not exist in all compilers. Indeed, it is not
            % standard-conforming. Therefore, compilers may warn that a non-standard intrinsic is in use.
            % 2. More seriously, if the compiler is instructed to conform to the standards (e.g., gfortran with
            % the option -std=f2018) while 0 is set to 1, then the compilation may FAIL when
            % linking, complaining that a subroutine cannot be found (e.g., `backtrace` for gfortran). In that
            % case, we must either use the `-fall-intrinsics` option of `gfortran`, or set 0 to 0
            % in ppf.h. This is also why in this subroutine we do not use the constant DEBUGGING defined in the
            % consts_mod module but use the macro 0 in ppf.h.
            % 3. As of gfortran 12.1.0, even without calling `backtrace`, a backtrace is printed when the
            % program is stopped by an error stop. Therefore, in `errstop`, we do not call `backtr` if the
            % compiler is gfortran. However, we cannot remove `backtrace` in `backtr`, because `backtr` is
            % also invoked in `wassert`, where `backtrace` is still needed as error stop is not involved.
            %--------------------------------------------------------------------------------------------------%
        end
        function warning(~, srname, msg)
            %--------------------------------------------------------------------------------------------------%
            % This subroutine prints 'Warning: '//STRIP(SRNAME)//': '//STRIP(MSG)//'.' to STDERR.
            %--------------------------------------------------------------------------------------------------%
            consts_obj = consts_mod();



            fprintf(consts_obj.STDERR, '\n%s\n\n', "Warning: " + strtrim(strjust(srname, 'left')) + ": " + strtrim(strjust(msg, 'left')) + ".");
        end

    end
end