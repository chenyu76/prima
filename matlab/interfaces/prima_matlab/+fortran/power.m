function y = power(x, n)
    if isnumeric(n) && isscalar(n) && isfinite(n)
    	if n == round(n)
            y = fortran_power_integer(x, double(n));
    	else
            y = fortran_power_real(x, double(n));
    else
        error('fortran:power:InvalidExponent', ...
            'power requires a finite scalar exponent.');
    end
end
