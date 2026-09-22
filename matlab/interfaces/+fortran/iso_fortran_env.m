classdef iso_fortran_env
    % ISO_FORTRAN_ENV Fortran intrinsic module mapped to MATLAB

    properties (Constant)
        real32      = 'single'
        real64      = 'double'
        real128     = 'double'
        int8        = 'int8'
        int16       = 'int16'
        int32       = 'int32'
        int64       = 'int64'

        % Unit numbers follow MATLAB's own file identifiers: 0 stdin, 1 stdout,
        % 2 stderr.  Every unit the Fortran module names must exist here: a USE
        % only list can name them, and the generated wrapper calls the member by
        % the spelling defined in this file.
        input_unit  = 0
        output_unit = 1
        error_unit  = 2
    end
end
