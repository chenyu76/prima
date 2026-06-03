function loc = maxloc(A, varargin)
    dim = [];
    mask = [];
    back = false;
    
    i = 1;
    while i <= length(varargin)
        arg = varargin{i};
        if ischar(arg) || isstring(arg)
            key = lower(char(arg));
            switch key
                case 'dim',  dim = varargin{i+1};  i = i + 2;
                case 'mask', mask = varargin{i+1}; i = i + 2;
                case 'back', back = varargin{i+1}; i = i + 2;
                case 'kind', i = i + 2; 
                otherwise,   i = i + 1;
            end
        else
            if isscalar(arg) && isnumeric(arg) && isempty(dim) && isempty(mask)
                dim = arg;
            elseif (islogical(arg) || numel(arg) > 1 || isempty(arg)) && isempty(mask)
                mask = arg;
            elseif isscalar(arg) && isnumeric(arg) && ~isempty(mask) && isempty(dim)
                dim = arg;
            end
            i = i + 1;
        end
    end
    
    is_fortran_1d = isvector(A) && ~isscalar(A);
    
    if isempty(A)
        loc = handleEmptyReturn(A, dim, is_fortran_1d);
        return;
    end
    
    A_masked = double(A);
    if ~isempty(mask)
        if ~any(mask(:))
            loc = handleEmptyReturn(A, dim, is_fortran_1d);
            return;
        end
        A_masked(~mask) = -Inf; 
    end
    
    if isempty(dim)
        maxVal = max(A_masked(:));
        
        if back
            linearIdx = find(A_masked(:) == maxVal, 1, 'last');
        else
            linearIdx = find(A_masked(:) == maxVal, 1, 'first');
        end
        
        if is_fortran_1d
            loc = linearIdx;
        else
            numDims = ndims(A);
            locCell = cell(1, numDims);
            [locCell{:}] = ind2sub(size(A), linearIdx);
            loc = cell2mat(locCell);
            loc = loc(:).'; 
        end
    else
        actual_dim = dim;
        
        if isrow(A) && dim == 1
            actual_dim = 2;
        elseif iscolumn(A) && dim == 1
            actual_dim = 1;
        end
        
        if back
            A_flipped = flip(A_masked, actual_dim);
            [~, loc_flipped] = max(A_flipped, [], actual_dim);
            loc = size(A, actual_dim) - loc_flipped + 1;
        else
            [~, loc] = max(A_masked, [], actual_dim);
        end
        
        if ~isempty(mask)
            all_masked = all(~mask, actual_dim);
            loc(all_masked) = 0;
        end
    end

    function loc = handleEmptyReturn(A, dim, is_1d)
        if isempty(dim)
            if is_1d
                loc = 0;
            else
                loc = zeros(1, ndims(A));
            end
        else
            if is_1d
                loc = 0;
            else
                sz = size(A);
                sz(dim) = [];
                if isempty(sz), sz = [1 1]; end
                loc = zeros(sz);
            end
        end
    end
end
