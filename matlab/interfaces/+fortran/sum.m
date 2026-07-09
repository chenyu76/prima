function y = sum(x, dim)
    if nargin < 2 || (ischar(dim) && strcmp(dim, 'all'))
        y = fortran_sum(x, 0);
    else
        y = fortran_sum(x, dim);
    end
end
