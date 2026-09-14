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
        function assert(obj)
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


            % A condition that is expected to be true
            % Description of the condition in human language
            % Name of the subroutine that calls this procedure

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

            % A condition that is expected to be true
            % Description of the condition in human language
            % Name of the subroutine that calls this procedure

        end
        function warning(~, srname, msg)
            %--------------------------------------------------------------------------------------------------%
            % This subroutine prints 'Warning: '//STRIP(SRNAME)//': '//STRIP(MSG)//'.' to STDERR.
            %--------------------------------------------------------------------------------------------------%


            fprintf(2, '\n%s\n\n', ...
                    "Warning: " + strtrim(strjust(srname, 'left')) + ": " ...
                    + strtrim(strjust(msg, 'left')) + ".");
        end

    end
end