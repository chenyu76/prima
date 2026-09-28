function result = format(operation, fmt, data, types, layout)
    % FORMAT transfers one formatted I/O statement that the translator could not
    % express with a MATLAB conversion.
    %
    %   operation  'write' or 'read'.
    %   fmt        the Fortran FORMAT text, for example '(I4,F8.2)'.
    %   data       for a write, the values: a cell array with one cell per I/O item,
    %              so an array item is one cell holding its elements.  For a read,
    %              the records to read: a character array for an internal file, or
    %              the file identifier of an external one.
    %   types      one row per I/O item: {kind, kind bytes, element length}.
    %              kind is 'integer', 'real', 'complex', 'character' or 'logical'.
    %              kind bytes is a byte count, or the class name the storage model
    %              gives an intrinsic-module kind constant (real64 is 'double').
    %              element length applies to a character item, and is [] otherwise.
    %   layout     the file identifier of an external unit; [recordLength,
    %              recordCount] for a write to an internal file; the shape of each
    %              item for a read.
    %
    % Returns the records written to an internal file, [] for an external write, and
    % the values read, one cell per I/O item.
    %
    % A pass of the format is one record: a format that runs out is reused from its
    % last group, and values are transferred in array element order.  State is
    % call-local, and the nested helpers keep single-file embedding out of the
    % caller's namespace.
    %
    % An error carries a fortran:format:* identifier, so a caller can tell a bad
    % FORMAT from the end of a file.
    writing = strcmp(operation, 'write');
    if ~writing && ~strcmp(operation, 'read')
        error('fortran:format:operation', 'Expected read or write.');
    end
    source = compactFormat(char(fmt));
    cursor = 1;
    skip();
    if cursor > length(source) || source(cursor) ~= '('
        error('fortran:format:syntax', 'FORMAT must be parenthesized.');
    end
    cursor = cursor + 1;
    nodes = parseGroup();
    revert = 1;
    for n = 1:numel(nodes)
        if strcmp(nodes{n}.code, 'group')
            revert = n;
        end
    end
    scale = 0;
    plus = false;
    blankZero = false;
    item = 1;
    element = 1;
    imaginary = false;
    record = '';
    column = 1;
    recordIndex = 1;
    loaded = false;
    internal = (writing && numel(layout) == 2) || (~writing && ~isnumeric(data));
    if writing
        values = data;
        if internal
            recordLength = layout(1);
            recordCount = layout(2);
            records = repmat({repmat(' ', 1, recordLength)}, recordCount, 1);
        end
    else
        values = cell(1, size(types, 1));
        for n = 1:numel(values)
            if strcmp(types{n, 1}, 'character')
                values{n} = strings(layout{n});
            elseif strcmp(types{n, 1}, 'logical')
                values{n} = false(layout{n});
            else
                values{n} = zeros(layout{n});
            end
        end
        if internal
            records = cellstr(string(data(:)));
        end
    end
    % A kind entry is either a byte count or the class name the storage model
    % gives an intrinsic-module kind constant; the rest of this function works in
    % bytes.
    for n = 1:size(types, 1)
        types{n, 2} = byteCount(types{n, 2});
    end
    normalizeItem();
    if ~writing
        loadRecord();
    end
    stopped = execute(nodes);
    while ~stopped && item <= numel(values)
        oldItem = item;
        oldElement = element;
        oldImaginary = imaginary;
        advance();
        stopped = execute(nodes(revert:end));
        if item == oldItem && element == oldElement && imaginary == oldImaginary
            error('fortran:format:noData', 'Format reversion has no data edit descriptor.');
        end
    end
    if writing
        saveRecord();
        if internal
            result = string(records);
        else
            result = [];
        end
    else
        result = values;
    end

    function text = compactFormat(raw)
        text = '';
        quote = char(0);
        for character = raw
            if quote ~= char(0)
                text(end + 1) = character;
                if character == quote
                    quote = char(0);
                end
            elseif character == '''' || character == '"'
                quote = character;
                text(end + 1) = character;
            elseif ~isspace(character)
                text(end + 1) = character;
            end
        end
    end

    function skip()
        while cursor <= length(source) && any(source(cursor) == [',' ' ' char(9) char(10)])
            cursor = cursor + 1;
        end
    end

    function number = digits()
        start = cursor;
        while cursor <= length(source) && source(cursor) >= '0' && source(cursor) <= '9'
            cursor = cursor + 1;
        end
        if start == cursor
            number = [];
        else
            number = str2double(source(start:cursor - 1));
        end
    end

    function list = parseGroup()
        list = {};
        while true
            skip();
            if cursor > length(source)
                error('fortran:format:syntax', 'Unclosed FORMAT group.');
            end
            c = source(cursor);
            if c == ')'
                cursor = cursor + 1;
                return
            end
            node = struct('code', '', 'repeat', 1, 'width', [], 'precision', [], 'exponent', [], 'text', '', 'children', {{}});
            if c == '''' || c == '"'
                cursor = cursor + 1;
                while cursor <= length(source)
                    ch = source(cursor);
                    cursor = cursor + 1;
                    if ch == c
                        if cursor <= length(source) && source(cursor) == c
                            cursor = cursor + 1;
                        else
                            break
                        end
                    end
                    node.text(end + 1) = ch;
                end
                node.code = 'literal';
            else
                sign = 1;
                if c == '-' || c == '+'
                    if c == '-'
                        sign = -1;
                    end
                    cursor = cursor + 1;
                end
                if source(cursor) == '*'
                    repeat = Inf;
                    cursor = cursor + 1;
                else
                    repeat = digits();
                    if isempty(repeat)
                        repeat = 1;
                    end
                    repeat = sign * repeat;
                end
                node.repeat = repeat;
                c = upper(source(cursor));
                cursor = cursor + 1;
                if c == '('
                    node.code = 'group';
                    node.children = parseGroup();
                else
                    node.code = c;
                    if cursor <= length(source)
                        pair = [c upper(source(cursor))];
                        if any(strcmp(pair, {'ES', 'EN', 'SP', 'SS', 'BN', 'BZ', 'TL', 'TR'}))
                            node.code = pair;
                            cursor = cursor + 1;
                        end
                    end
                    node.width = digits();
                    if cursor <= length(source) && source(cursor) == '.'
                        cursor = cursor + 1;
                        node.precision = digits();
                    end
                    if any(strcmp(node.code, {'E', 'ES', 'EN', 'G'})) && cursor <= length(source) && upper(source(cursor)) == 'E'
                        cursor = cursor + 1;
                        node.exponent = digits();
                    end
                    if ~any(strcmp(node.code, {'I', 'B', 'O', 'Z', 'F', 'E', 'ES', 'EN', 'D', 'G', 'A', 'L', 'X', 'T', 'TL', 'TR', '/', ':', 'P', 'S', 'SP', 'SS', 'BN', 'BZ'}))
                        error('fortran:format:descriptor', 'Unsupported FORMAT descriptor %s.', node.code);
                    end
                end
            end
            list{end + 1} = node;
        end
    end

    function normalizeItem()
        while item <= numel(values) && element > numel(values{item})
            item = item + 1;
            element = 1;
        end
    end

    function stopped = execute(list)
        stopped = false;
        for j = 1:numel(list)
            node = list{j};
            code = node.code;
            switch code
                case 'group'
                    iteration = 0;
                    while iteration < node.repeat
                        before = [item element imaginary];
                        if execute(node.children)
                            stopped = true;
                            return
                        end
                        iteration = iteration + 1;
                        if isinf(node.repeat) && isequal(before, [item element imaginary])
                            error('fortran:format:noData', 'Unlimited group has no data transfer.');
                        end
                    end
                case 'literal'
                    if writing
                        put(node.text);
                    else
                        column = column + length(node.text);
                    end
                case 'P'
                    scale = node.repeat;
                case 'SP'
                    plus = true;
                case {'S', 'SS'}
                    plus = false;
                case 'BN'
                    blankZero = false;
                case 'BZ'
                    blankZero = true;
                case ':'
                    if item > numel(values)
                        stopped = true;
                        return
                    end
                case '/'
                    for k = 1:node.repeat
                        advance();
                    end
                case 'X'
                    column = column + node.repeat;
                case 'T'
                    column = node.width;
                case 'TL'
                    column = max(1, column - node.width);
                case 'TR'
                    column = column + node.width;
                otherwise
                    for k = 1:node.repeat
                        if item > numel(values)
                            stopped = true;
                            return
                        end
                        kind = types{item, 1};
                        bytes = types{item, 2};
                        len = types{item, 3};
                        if ~any(strcmp(kind, {'integer', 'real', 'complex', 'character', 'logical'}))
                            error('fortran:format:type', 'Unsupported formatted I/O type %s.', kind);
                        end
                        if any(strcmp(kind, {'real', 'complex'})) && ~any(bytes == [4 8])
                            error('fortran:format:kind', 'Unsupported real kind %d.', bytes);
                        end
                        if writing
                            value = values{item}(element);
                            if strcmp(kind, 'complex')
                                if imaginary
                                    value = imag(value);
                                else
                                    value = real(value);
                                end
                            end
                            if any(strcmp(kind, {'real', 'complex'})) && bytes == 4
                                value = double(single(value));
                            end
                            put(writeField(node, value, kind, bytes, len));
                        else
                            width = node.width;
                            if isempty(width) && strcmp(code, 'A')
                                width = len;
                            end
                            if isempty(width) || width == 0
                                error('fortran:format:width', 'Input requires a positive field width.');
                            end
                            loadRecord();
                            padded = [record repmat(' ', 1, max(0, column + width - 1 - length(record)))];
                            field = padded(column:column + width - 1);
                            column = column + width;
                            value = readField(node, field, kind, bytes, len);
                            if strcmp(kind, 'complex')
                                if imaginary
                                    values{item}(element) = complex(real(values{item}(element)), value);
                                else
                                    values{item}(element) = value;
                                end
                            else
                                values{item}(element) = value;
                            end
                        end
                        if strcmp(kind, 'complex') && ~imaginary
                            imaginary = true;
                        else
                            imaginary = false;
                            element = element + 1;
                            normalizeItem();
                        end
                    end
            end
        end
    end

    function loadRecord()
        if loaded
            return
        end
        if internal
            if recordIndex > numel(records)
                error('fortran:format:end', 'End of internal file.');
            end
            record = records{recordIndex};
        else
            record = fgetl(data);
            if ~ischar(record)
                error('fortran:format:end', 'End of file.');
            end
        end
        loaded = true;
    end

    function saveRecord()
        if internal
            if recordIndex > recordCount
                error('fortran:format:end', 'End of internal file.');
            end
            records{recordIndex} = [record repmat(' ', 1, recordLength - length(record))];
        else
            fprintf(layout, '%s\n', record);
        end
    end

    function advance()
        if writing
            saveRecord();
        end
        recordIndex = recordIndex + 1;
        column = 1;
        record = '';
        loaded = false;
        % Internal trailing slashes need not access another record unless a
        % subsequent data descriptor actually transfers a field.
        if ~writing && ~internal
            loadRecord();
        end
    end

    function put(text)
        last = column + length(text) - 1;
        if internal && last > recordLength
            error('fortran:format:record', 'Formatted output exceeds internal record length.');
        end
        if isempty(text)
            return
        end
        record(end + 1:last) = ' ';
        record(column:last) = text;
        column = last + 1;
    end

    function text = writeField(node, value, kind, bytes, len)
        code = node.code;
        w = node.width;
        d = node.precision;
        if strcmp(code, 'A')
            text = char(string(value));
            text = [text repmat(' ', 1, max(0, len - length(text)))];
            text = text(1:len);
            if ~isempty(w)
                if length(text) > w
                    text = text(1:w);
                else
                    text = [repmat(' ', 1, w - length(text)) text];
                end
            end
            return
        elseif strcmp(code, 'L')
            if value
                text = 'T';
            else
                text = 'F';
            end
        elseif any(strcmp(code, {'B', 'O', 'Z'}))
            bits = toBits(value, kind, bytes);
            base = radixBase(code);
            alphabet = '0123456789ABCDEF';
            text = '';
            while bits > 0
                digit = bitand(bits, uint64(base - 1));
                text = [alphabet(double(digit) + 1) text];
                bits = bitshift(bits, -log2(base));
            end
            if isempty(text)
                text = '0';
            end
            if ~isempty(d)
                if value == 0 && d == 0
                    text = '';
                end
                text = [repmat('0', 1, max(0, d - length(text))) text];
            end
        elseif strcmp(code, 'I')
            text = sprintf('%.0f', abs(value));
            if ~isempty(d)
                if value == 0 && d == 0
                    text = '';
                end
                text = [repmat('0', 1, max(0, d - length(text))) text];
            end
            if value < 0
                text = ['-' text];
            elseif plus && ~isempty(text)
                text = ['+' text];
            end
        else
            if strcmp(code, 'G') && isempty(d)
                if bytes == 4
                    d = 9;
                else
                    d = 17;
                end
            end
            if isempty(d)
                error('fortran:format:precision', 'Real edit descriptor requires precision.');
            end
            negative = value < 0 || (value == 0 && 1 / value < 0);
            x = abs(double(value));
            if ~isfinite(x)
                if isnan(x)
                    text = 'NaN';
                else
                    text = 'Infinity';
                end
                signWidth = double(negative || plus);
                if isinf(x) && ~isempty(w) && w > 0 && w < 8 + signWidth
                    text = 'Inf';
                end
            elseif strcmp(code, 'F')
                % Shift a decimal representation before rounding, avoiding binary
                % multiplication by 10^scale (which can overflow or round twice).
                text = fixedDecimal(x, d, scale);
            else
                text = exponential(node, x, d);
            end
            if negative
                text = ['-' text];
            elseif plus
                text = ['+' text];
            end
            % gfortran drops the leading zero of the mantissa in the F0.d form, and
            % whenever dropping it is what makes an overlong field fit ('(E7.2)' of
            % 0.5 writes ".50E+00").
            if d > 0 && ((strcmp(code, 'F') && isequal(w, 0)) || (~isempty(w) && w > 0 && length(text) > w))
                text = regexprep(text, '^([+-]?)0\.', '$1.');
            end
        end
        if isempty(text) && isequal(w, 0)
            text = ' ';
        end
        if ~isempty(w) && w > 0
            if length(text) > w
                text = repmat('*', 1, w);
            else
                text = [repmat(' ', 1, w - length(text)) text];
            end
        end
    end

    function text = fixedDecimal(x, d, shift)
        if shift == 0 || x == 0
            text = sprintf('%.*f', d, x);
        else
            % Scientific printf rounds to the requested significant digits. Moving
            % its decimal point afterwards is exact string manipulation.
            probe = sprintf('%.17e', x);
            split = strfind(probe, 'e');
            e = str2double(probe(split + 1:end)) + shift;
            significant = e + d + 1;
            if significant <= 0
                digit = '0';
                if significant == 0 && (probe(1) > '5' || (probe(1) == '5' && any(probe(3:split - 1) ~= '0')))
                    digit = '1';
                end
                if d == 0
                    text = [digit '.'];
                else
                    text = ['0.' repmat('0', 1, d - 1) digit];
                end
            else
                rounded = sprintf('%.*e', significant - 1, x);
                split = strfind(rounded, 'e');
                e = str2double(rounded(split + 1:end)) + shift;
                digitsText = strrep(rounded(1:split - 1), '.', '');
                point = e + 1;
                digitsText = [repmat('0', 1, max(0, 1 - point)) digitsText];
                point = max(1, point);
                digitsText = [digitsText repmat('0', 1, max(0, point + d - length(digitsText)))];
                text = [digitsText(1:point) '.' digitsText(point + 1:end)];
            end
        end
        if d == 0 && ~contains(text, '.')
            text = [text '.'];
        end
    end

    function text = exponential(node, x, d)
        % d is the digit count resolved by the caller: G without one keeps as many
        % digits as the kind carries, which the node itself does not say.
        code = node.code;
        probe = sprintf('%.17e', x);
        split = strfind(probe, 'e');
        power = str2double(probe(split + 1:end));
        if strcmp(code, 'G')
            probe = sprintf('%.*e', d - 1, x);
            split = strfind(probe, 'e');
            power = str2double(probe(split + 1:end));
        end
        if strcmp(code, 'G') && (x == 0 || (power >= -1 && power < d))
            decimals = max(0, d - power - 1);
            if x == 0
                decimals = max(0, d - 1);
            end
            text = fixedDecimal(x, decimals, 0);
            if ~isempty(node.width) && node.width > 0
                trailing = 4;
                if ~isempty(node.exponent)
                    trailing = node.exponent + 2;
                end
                text = [text repmat(' ', 1, trailing)];
            end
            return
        end
        k = scale;
        if strcmp(code, 'ES')
            k = 1;
        end
        if strcmp(code, 'EN')
            k = mod(power, 3) + 1;
        end
        % ES and EN place the point themselves, so their digits come from d and k
        % rather than from the scale factor's own rule.
        if strcmp(code, 'ES')
            significant = d + 1;
        elseif strcmp(code, 'EN')
            significant = d + k;
        elseif k > 0
            significant = d + 1;
        elseif k < 0
            significant = d + k;
        else
            significant = d;
        end
        if significant < 1
            error('fortran:format:scale', 'Scale factor incompatible with precision.');
        end
        rounded = sprintf('%.*e', significant - 1, x);
        split = strfind(rounded, 'e');
        newPower = str2double(rounded(split + 1:end));
        if strcmp(code, 'EN') && newPower ~= power
            k = mod(newPower, 3) + 1;
        end
        digitsText = strrep(rounded(1:split - 1), '.', '');
        if strcmp(code, 'EN')
            digitsText = [digitsText repmat('0', 1, max(0, d + k - length(digitsText)))];
            digitsText = digitsText(1:d + k);
        end
        exponent = newPower - k + 1;
        if x == 0
            exponent = 0;
        end
        if k <= 0
            mantissa = ['0.' repmat('0', 1, -k) digitsText];
        else
            digitsText = [digitsText repmat('0', 1, max(0, k - length(digitsText)))];
            mantissa = [digitsText(1:k) '.' digitsText(k + 1:end)];
        end
        ew = node.exponent;
        letter = 'E';
        if strcmp(code, 'D')
            letter = 'D';
        end
        if isempty(ew)
            ew = 2;
            if isequal(node.width, 0)
                ew = max(3, length(sprintf('%d', abs(exponent))));
            elseif abs(exponent) >= 100
                ew = 3;
                letter = '';
            end
        elseif ew == 0
            ew = length(sprintf('%d', abs(exponent)));
        elseif length(sprintf('%d', abs(exponent))) > ew
            text = repmat('*', 1, node.width);
            return
        end
        sign = '+';
        if exponent < 0
            sign = '-';
        end
        text = [mantissa letter sign sprintf('%0*d', ew, abs(exponent))];
    end

    function value = readField(node, field, kind, bytes, len)
        code = node.code;
        if strcmp(code, 'A')
            if length(field) > len
                field = field(end - len + 1:end);
            end
            value = string([field repmat(' ', 1, max(0, len - length(field)))]);
            return
        end
        if strcmp(code, 'L')
            token = upper(strtrim(field));
            token = regexprep(token, '^\.', '');
            if isempty(token) || ~any(token(1) == 'TF')
                error('fortran:format:logical', 'Invalid logical field.');
            end
            value = token(1) == 'T';
            return
        end
        field = regexprep(field, '^ +', '');
        if blankZero
            field = strrep(field, ' ', '0');
        else
            field = strrep(field, ' ', '');
        end
        if isempty(field)
            value = 0;
            return
        end
        if any(strcmp(code, {'B', 'O', 'Z'}))
            base = radixBase(code);
            bits = uint64(0);
            alphabet = '0123456789ABCDEF';
            for ch = upper(field)
                digit = find(alphabet == ch, 1) - 1;
                if isempty(digit) || digit >= base
                    error('fortran:format:number', 'Invalid radix field.');
                end
                bits = bitor(bitshift(bits, log2(base)), uint64(digit));
            end
            value = fromBits(bits, kind, bytes);
            return
        end
        field = upper(strrep(field, 'd', 'E'));
        field = strrep(field, 'D', 'E');
        if any(strcmp(field, {'INF', '+INF', 'INFINITY', '+INFINITY'}))
            value = Inf;
            return
        end
        if any(strcmp(field, {'-INF', '-INFINITY'}))
            value = -Inf;
            return
        end
        field = regexprep(field, '(?<=[0-9.])([+-][0-9]+)$', 'E$1');
        parts = strsplit(field, 'E');
        adjustment = 0;
        if ~strcmp(code, 'I') && ~contains(parts{1}, '.') && ~isempty(node.precision)
            adjustment = -node.precision;
        end
        if numel(parts) == 1 && ~strcmp(code, 'I')
            adjustment = adjustment - scale;
        end
        if numel(parts) == 2
            adjustment = adjustment + str2double(parts{2});
        end
        value = str2double([parts{1} 'e' num2str(adjustment)]);
        if isnan(value) && ~contains(field, 'NAN')
            error('fortran:format:number', 'Invalid numeric field %s.', field);
        end
        if any(strcmp(kind, {'real', 'complex'})) && bytes == 4
            value = double(single(value));
        end
    end

    function bytes = byteCount(entry)
        if isnumeric(entry)
            bytes = entry;
            return
        end
        names = {'int8', 'int16', 'int32', 'int64', 'single', 'double'};
        sizes = [1 2 4 8 4 8];
        index = find(strcmp(char(entry), names), 1);
        if isempty(index)
            error('fortran:format:kind', 'Unsupported kind %s.', char(entry));
        end
        bytes = sizes(index);
    end

    function base = radixBase(code)
        if strcmp(code, 'B')
            base = 2;
        elseif strcmp(code, 'O')
            base = 8;
        else
            base = 16;
        end
    end

    function bits = toBits(value, kind, bytes)
        if any(strcmp(kind, {'real', 'complex'}))
            if bytes == 4
                bits = uint64(typecast(single(value), 'uint32'));
            else
                bits = typecast(double(value), 'uint64');
            end
        else
            signed = feval(['int' num2str(8 * bytes)], value);
            bits = uint64(typecast(signed, ['uint' num2str(8 * bytes)]));
        end
    end

    function value = fromBits(bits, kind, bytes)
        unsigned = feval(['uint' num2str(8 * bytes)], bits);
        if any(strcmp(kind, {'real', 'complex'}))
            if bytes == 4
                value = double(typecast(unsigned, 'single'));
            else
                value = typecast(unsigned, 'double');
            end
        else
            value = double(typecast(unsigned, ['int' num2str(8 * bytes)]));
        end
    end

end
