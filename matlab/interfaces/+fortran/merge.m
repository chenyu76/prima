function result = merge(varargin)
    if nargin == 3 && ~ischar(varargin{1}) && ~isstring(varargin{1})
        tsource = varargin{1};
        fsource = varargin{2};
        mask = varargin{3};
    else
        tsource = [];
        fsource = [];
        mask = [];
        i = 1;
        while i <= length(varargin)
            arg = varargin{i};
            if ischar(arg) || isstring(arg)
                key = lower(char(arg));
                switch key
                    case 'tsource'
                        tsource = varargin{i + 1};
                        i = i + 2;
                    case 'fsource'
                        fsource = varargin{i + 1};
                        i = i + 2;
                    case 'mask'
                        mask = varargin{i + 1};
                        i = i + 2;
                    otherwise
                        if isempty(tsource)
                            tsource = arg;
                        elseif isempty(fsource)
                            fsource = arg;
                        elseif isempty(mask)
                            mask = arg;
                        end
                        i = i + 1;
                end
            else
                if isempty(tsource)
                    tsource = arg;
                elseif isempty(fsource)
                    fsource = arg;
                elseif isempty(mask)
                    mask = arg;
                end
                i = i + 1;
            end
        end
    end

    if isscalar(mask)
        if mask
            if isscalar(tsource) && ~isscalar(fsource)
                result = repmat(tsource, size(fsource));
            else
                result = tsource;
            end
        else
            if isscalar(fsource) && ~isscalar(tsource)
                result = repmat(fsource, size(tsource));
            else
                result = fsource;
            end
        end
        return
    end

    if isscalar(tsource)
        tsource = repmat(tsource, size(mask));
    end
    if isscalar(fsource)
        fsource = repmat(fsource, size(mask));
    end

    result = fsource;
    result(mask) = tsource(mask);
end
