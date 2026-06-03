classdef fprint_mod
    %--------------------------------------------------------------------------------------------------%
    % This module provides a subroutine that prints a string to STDOUT, STDERR, or a normal file.
    %
    % N.B.: When interfacing the code with MATLAB, this module needs to be revised to use the MATLAB
    % MEX function mexPrintf instead of WRITE. This is because the Fortran WRITE cannot write to the
    % STDOUT when the code is interfaced with MATLAB, since the STDOUT is hijacked by MEX. See
    % https://stackoverflow.com/questions/26271154/how-can-i-make-a-mex-function-printf-while-its-running
    % https://www.mathworks.com/matlabcentral/answers/132527-in-mex-files-where-does-output-to-stdout-and-stderr-go
    %
    % Coded by Zaikun ZHANG (www.zhangzk.net) based on Powell's code and papers.
    %
    % Started: July 2020
    %
    % Last Modified: Sunday, May 21, 2023 AM01:29:25
    %--------------------------------------------------------------------------------------------------%

    methods
        function fprint(~, string, varargin)
            consts_obj = consts_mod();
            debug_obj = debug_mod();
            string_obj = string_mod();

            % Inputs





            % Local variables
            newline = compose('\n');
            srname = "FPRINT";
            fname_loc = "";
            fstat = "";
            position = "";
            funit_loc = NaN;
            i = NaN;
            iostat = NaN;
            j = NaN;
            slen = NaN;
            fexist = false;

            % Preconditions
            ipObj = inputParser();
            addParameter(ipObj, 'funit', NaN);
            addParameter(ipObj, 'fname', "");
            addParameter(ipObj, 'faction', "");
            parse(ipObj, varargin{:});
            funit = ipObj.Results.funit;
            fname = ipObj.Results.fname;
            faction = ipObj.Results.faction;
            if consts_obj.DEBUGGING
                if ~ismember('funit', ipObj.UsingDefaults)
                    debug_obj.assert(funit ~= consts_obj.STDIN, "The file unit is not STDIN", srname);
                    if ~ismember('fname', ipObj.UsingDefaults)
                        debug_obj.assert((strlength(fname) == 0) == (funit == consts_obj.STDOUT || funit == consts_obj.STDERR), "The file name is empty if and only if the file unit is either STDOUT or STDERR", srname);
                    end
                end
                if ~ismember('faction', ipObj.UsingDefaults)
                    debug_obj.assert(faction == "write" || faction == "w" || faction == "append" || faction == "a", "FACTION is either ""write (w)"" or ""append (a)""", srname);
                end
            end

            %====================%
            % Calculation starts %
            %====================%

            % Decide the file storage unit.
            if ismember('funit', ipObj.UsingDefaults)
                if ismember('fname', ipObj.UsingDefaults)
                    funit_loc = consts_obj.STDOUT; % Print the message to the standard out.
                else
                    funit_loc = -1; % This value will not be used.

                end
            else
                funit_loc = funit;
            end

            % Decide the file name.
            if ~ismember('fname', ipObj.UsingDefaults)
                fname_loc = fname;
            elseif funit_loc ~= consts_obj.STDOUT && funit_loc ~= consts_obj.STDERR
                fname_loc = "fort." + string_obj.int2str(fix(funit_loc));
            else
                fname_loc = "";
            end

            if consts_obj.DEBUGGING
                debug_obj.assert((strlength(fname_loc) == 0) == (funit_loc == consts_obj.STDOUT || funit_loc == consts_obj.STDERR), "The file name is empty if and only if the file unit is either STDOUT or STDERR", srname);
            end

            % Open the file if necessary.
            iostat = 0;
            if strlength(fname_loc) > 0
                % Decide the position for OPEN. This is the only place where FACTION is used.
                position = "append";
                if ~ismember('faction', ipObj.UsingDefaults)
                    switch faction
                    case {"write", "w"}
                        position = "rewind";
                    case {"append", "a"}
                        position = "append";
                    otherwise
                        debug_obj.warning(srname, "Unknown file action """ + faction + """");
                    end
                end
                % Check whether the file is already existing.
                %Unsupported Statement: StmtExpr (CallFunction "inquire" [KeywordArg "file" (Variable "fname_loc"),KeywordArg "exist" (Variable "fexist")])

                fstat = fortran.merge('tsource', "old", 'fsource', "new", 'mask', fexist);
                % Open the file.

                fid_funit_loc = fopen(fname_loc, 'w');

                if iostat ~= 0
                    debug_obj.warning(srname, "Failed to open file " + fname_loc);
                    return
                end
            end

            % Print the string.
            % N.B.: `WRITE (FUNIT_LOC, '(A)') STRING` would do what we want, but it causes "Buffer overflow on
            % output" if string is long. This did occur with NAG Fortran Compiler R7.1(Hanzomon) Build 7122.
            % To avoid this problem, we print the string line by line, separated by newlines.
            i = 1;
            j = fortran.index(string, newline); % Index of the first newline in the string.
            slen = strlength(string);
            while j >= i                % J < I: No more newline in the string.
                fprintf(funit_loc, '%s\n', extractBetween(string, i, j - 1)); % Print the string before the current newline.
                i = j + 1; % Index of the character after the current newline.
                j = i + fortran.index(extractBetween(string, i, slen), newline) - 1; % Index of the next newline.
            end
            if extractBetween(string, i, slen) ~= ""
                % Print the string after the last newline.
                fprintf(funit_loc, '%s\n', extractBetween(string, i, slen));
            end

            % Close the file if necessary.
            if strlength(fname_loc) > 0 && iostat == 0
                fclose(fid_funit_loc);
            end

            %====================%
            %  Calculation ends  %
            %====================%
        end

    end
end