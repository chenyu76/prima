function result = spread(varargin)
    if nargin == 3 && ~ischar(varargin{2}) && ~isstring(varargin{2})
        source = varargin{1};
        dim = varargin{2};
        ncopies = varargin{3};
    else
        source = []; dim = []; ncopies = [];
        i = 1;
        while i <= length(varargin)
            arg = varargin{i};
            if ischar(arg) || isstring(arg)
                key = lower(char(arg));
                switch key
                    case 'source',  source = varargin{i+1}; i = i + 2;
                    case 'dim',     dim = varargin{i+1}; i = i + 2;
                    case 'ncopies', ncopies = varargin{i+1}; i = i + 2;
                    otherwise
                        if isempty(source), source = arg;
                        elseif isempty(dim), dim = arg;
                        elseif isempty(ncopies), ncopies = arg;
                        end
                        i = i + 1;
                end
            else
                if isempty(source), source = arg;
                elseif isempty(dim), dim = arg;
                elseif isempty(ncopies), ncopies = arg;
                end
                i = i + 1;
            end
        end
    end

    sz = size(source);

    if length(sz) == 2 
        if sz(1) == 1 && sz(2) > 1 
            sz = sz(2);
            source = source(:);
        elseif sz(2) == 1
            sz = sz(1);
        elseif sz(1) == 1 && sz(2) == 1 
            sz = 1;
        end
    end
    
    if ncopies <= 0
        empty_sz = [sz(1:dim-1), 0, sz(dim:end)];
        if length(empty_sz) == 1
            empty_sz = [empty_sz, 1];
        end
        result = zeros(empty_sz, 'like', source); 
        return;
    end
    
    new_sz = [sz(1:dim-1), 1, sz(dim:end)];
    
    if length(new_sz) == 1
        new_sz = [new_sz, 1];
    end
    
    temp = reshape(source, new_sz);
    
    rep_vec = ones(1, length(new_sz));
    rep_vec(dim) = ncopies;
    
    result = repmat(temp, rep_vec);
end
