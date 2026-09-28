function backtrace()
    % BACKTRACE Print the current MATLAB call stack without stopping execution.

    stack = dbstack(1, '-completenames');

    fprintf(2, 'Backtrace:\n');
    for k = 1:numel(stack)
        fprintf(2, '  at %s (%s:%d)\n', ...
                stack(k).name, stack(k).file, stack(k).line);
    end
end
