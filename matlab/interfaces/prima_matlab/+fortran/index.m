function pos = index(varargin)
    if nargin == 2 || (nargin == 3 && ~ischar(varargin{3}) && ~isstring(varargin{3}))
        s = varargin{1};
        sub = varargin{2};
        if nargin == 3
            back = varargin{3};
        else
            back = false;
        end
    else
        s = []; sub = []; back = false;
        i = 1;
        while i <= length(varargin)
            arg = varargin{i};
            if ischar(arg) || isstring(arg)
                key = lower(char(arg));
                switch key
                    case 'string',    s = varargin{i+1}; i = i + 2;
                    case 'substring', sub = varargin{i+1}; i = i + 2;
                    case 'back',      back = varargin{i+1}; i = i + 2;
                    case 'kind',      i = i + 2;
                    otherwise
                        if isempty(s), s = arg;
                        elseif isempty(sub), sub = arg;
                        else back = arg;
                        end
                        i = i + 1;
                end
            else
                if isempty(s), s = arg;
                elseif isempty(sub), sub = arg;
                else back = arg;
                end
                i = i + 1;
            end
        end
    end

    if isempty(s) || isempty(sub)
        pos = 0; return;
    end

    try
        jS = java.lang.String(s);
        jSub = java.lang.String(sub);
        
        if back
            jPos = jS.lastIndexOf(jSub);
        else
            jPos = jS.indexOf(jSub);
        end
        pos = double(jPos) + 1;
        
    catch
        if back
            res = regexp(s, regexptranslate('escape', sub), 'start', 'last');
        else
            res = regexp(s, regexptranslate('escape', sub), 'start', 'once');
        end
        
        if isempty(res)
            pos = 0;
        else
            pos = res;
        end
    end
end
