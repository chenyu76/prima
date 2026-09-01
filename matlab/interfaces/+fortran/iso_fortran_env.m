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

        output_unit = 1
        error_unit  = 2
    end
end
