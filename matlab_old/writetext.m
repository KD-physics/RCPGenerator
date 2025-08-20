function writetext(data, filename)
%WRITETEXT  Save matrix to a text file with full double precision.
%   writetext(DATA, FILENAME) writes the contents of DATA to the file
%   specified by FILENAME. Each row of DATA becomes one line in the file,
%   and entries are printed with enough precision (%.17g) to reconstruct
%   IEEE‑754 double values exactly.

    % Open file for writing
    fid = fopen(filename, 'w');
    if fid < 0
        error('writetext:cannotOpenFile', ...
              'Could not open file "%s" for writing.', filename);
    end

    % Determine number of columns
    [~, M] = size(data);

    % Build format string: '%.17g ' for each column except last, then '%.17g\n'
    fmt = [repmat('%.17g ', 1, M-1), '%.17g\n'];

    % fprintf prints column‑major, so transpose data
    fprintf(fid, fmt, data.');

    % Close file
    fclose(fid);
end
