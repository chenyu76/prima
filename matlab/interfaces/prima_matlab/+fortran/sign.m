% Fortran-compatible SIGN(A, B) function.
%
% Returns |A| with the sign of B, equivalent to Fortran's SIGN intrinsic.
%
% Fortran's SIGN(A, B) strictly depends the sign bit of B,
% while MATLAB's built-in sign() differs in two ways:
%   1. sign(x) returns -1, 0, or 1 based on a single argument,
%   2. sign(0.0) and sign(-0.0) both return 0, ignoring the IEEE 754 sign bit.
%
% To detect -0.0 in MATLAB, we use the trick that
%   1/(-0.0) == -Inf   and   1/0.0 == Inf.
function result = sign(a, b)
    if b < 0 || (b == 0 && 1 / b < 0)
        result = -abs(a);
    else
        result = abs(a);
    end
end
