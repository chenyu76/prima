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


            newline_custom = newline;


            fexist = false;


            ipObj = inputParser();
            addParameter(ipObj, 'funit', NaN);
            addParameter(ipObj, 'fname', "");
            addParameter(ipObj, 'faction', "");
            parse(ipObj, varargin{:});
            funit = ipObj.Results.funit;
            fname = ipObj.Results.fname;
            faction = ipObj.Results.faction;


            %====================%
            % Calculation starts %
            %====================%

            % Decide the file storage unit.
            if ismember('funit', ipObj.UsingDefaults)
                if ismember('fname', ipObj.UsingDefaults)
                    funit_loc = 1; % Print the message to the standard out.
                else
                    funit_loc = -1; % This value will not be used.

                end
            else
                funit_loc = funit;
            end

            % Decide the file name.
            if ~ismember('fname', ipObj.UsingDefaults)
                fname_loc = fname;
            elseif funit_loc ~= 1 && funit_loc ~= 2
                fname_loc = "fort." + int2str(funit_loc);
            else
                fname_loc = "";
            end


            % Open the file if necessary.
            iostat = 0;
            if strlength(fname_loc) > 0
                % Decide the position for OPEN. This is the only place where FACTION is used.

                if ~ismember('faction', ipObj.UsingDefaults)
                    switch faction
                    case {"write", "w"}

                    case {"append", "a"}

                    otherwise

                    end
                end
                % Check whether the file is already existing.
                inquire('file', fname_loc, 'exist', fexist);
                fortran.merge('tsource', "old", 'fsource', "new", 'mask', fexist);
                % Open the file.

                funit_loc = fopen(fname_loc, 'w');

                if iostat ~= 0
                    return
                end
            end

            % Print the string.
            % N.B.: `WRITE (FUNIT_LOC, '(A)') STRING` would do what we want, but it causes "Buffer overflow on
            % output" if string is long. This did occur with NAG Fortran Compiler R7.1(Hanzomon) Build 7122.
            % To avoid this problem, we print the string line by line, separated by newlines.
            i = 1;
            j = fortran.index(string, newline_custom); % Index of the first newline in the string.
            slen = strlength(string);
            while j >= i                % J < I: No more newline in the string.
                fprintf(funit_loc, '%s\n', extractBetween(string, i, j - 1)); % Print the string before the current newline.
                i = j + 1; % Index of the character after the current newline.
                j = i + fortran.index(extractBetween(string, i, slen), newline_custom) - 1; % Index of the next newline.
            end
            if extractBetween(string, i, slen) ~= ""
                % Print the string after the last newline.
                fprintf(funit_loc, '%s\n', extractBetween(string, i, slen));
            end

            % Close the file if necessary.
            if strlength(fname_loc) > 0 && iostat == 0
                fclose(funit_loc);
            end

            %====================%
            %  Calculation ends  %
            %====================%
        end

    end
end