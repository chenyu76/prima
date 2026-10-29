function result = assignCharacterSubstring(value, lower, upper, replacement)
    % Inputs use the translator's string scalar/array representation.
    % A char vector represents neither a string-array shape nor its elements.
    % Fortran assignment replaces a fixed-width field, including its blanks.
    result = value;
    if ~isscalar(replacement)
        replacement = reshape(replacement, size(value));
    end
    for i = 1:numel(value)
        text = char(value(i));
        last = upper;
        % [] denotes the end of this element, not an array index.
        if isempty(last)
            last = length(text);
        end
        width = max(last - lower + 1, 0);
        if width == 0
            continue;
        end
        if isscalar(replacement)
            newtext = char(replacement);
        else
            newtext = char(replacement(i));
        end
        field = repmat(' ', 1, width);
        copied = min(length(newtext), width);
        field(1:copied) = newtext(1:copied);
        text(lower:last) = field;
        result(i) = string(text);
    end
end
