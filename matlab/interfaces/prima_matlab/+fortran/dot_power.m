function y = dot_power(x, n)
    if isnumeric(n) && isscalar(n) && isfinite(n) && abs(n - round(n)) < 1e-12
        y = int_pow(x, int64(round(n)));
    elseif isnumeric(n) && all(isfinite(n(:)) & abs(n(:) - round(n(:))) < 1e-12)
        y = int_pow_elems(x, int64(round(n)));
    else
        y = builtin('power', x, n);
    end
end

function y = int_pow(x, n)
    if n < 0
        y = 1.0 ./ int_pow(x, -n);
        return
    end
    if n == 0
        y = ones(size(x), 'like', x);
        return
    end
    if n == 1
        y = x;
        return
    end
    if n == 2
        y = x .* x;
        return
    end
    y = ones(size(x), 'like', x);
    base = x;
    while n > 0
        if bitand(n, 1)
            y = y .* base;
        end
        n = bitshift(n, -1);
        if n > 0
            base = base .* base;
        end
    end
end

function y = int_pow_elems(x, n)
    y = zeros(size(x), 'like', x);
    for i = 1:numel(x)
        y(i) = int_pow_elem(x(i), n(i));
    end
end

function y = int_pow_elem(x, n)
    if n < 0
        y = 1.0 / int_pow_elem(x, -n);
        return
    end
    y = 1.0;
    base = x;
    while n > 0
        if bitand(n, 1)
            y = y * base;
        end
        n = bitshift(n, -1);
        if n > 0
            base = base * base;
        end
    end
end
