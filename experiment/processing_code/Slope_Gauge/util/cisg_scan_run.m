function S = cisg_scan_run(coreview_num, raw_root, frame_size)
%CISG_SCAN_RUN List bad / headerless CISG frames of a run from the file sizes only (no data read).
%
%   S = cisg_scan_run(coreview_num)
%   S = cisg_scan_run(coreview_num, raw_root, frame_size)
%
% A good file is header + h*w*3 bytes (the most common size).  Reported:
%   empty      - smaller than h*w*3 (e.g. 20-byte header-only files): unrecoverable
%   headerless - exactly h*w*3 bytes: pixels without header (cisg_load_coreview(file, [h w]) reads them)
%   other      - any other size
% S fields: normal_size, frames (numbers), empty, headerless, other (frame numbers).
%
% frame_size defaults to [3072 4100] (Flare 12M125).

    if nargin < 2 || isempty(raw_root), raw_root = 'D:\HLAB_2026\Core'; end
    if nargin < 3 || isempty(frame_size), frame_size = [3072 4100]; end
    folder = fullfile(raw_root, sprintf('CoreView_%d', coreview_num), 'Flare 12M125 CCL C1112A00');
    L = dir(fullfile(folder, sprintf('CoreView_%d_Flare 12M125 CCL C1112A00_*.raw', coreview_num)));
    tok = regexp({L.name}, '_(\d+)\.raw$', 'tokens', 'once');
    frames = cellfun(@(t) str2double(t{1}), tok);
    sizes = [L.bytes];
    pix = prod(frame_size) * 3;
    normal_size = mode(sizes);
    S.normal_size = normal_size;
    S.frames = frames;
    S.empty      = frames(sizes < pix);
    S.headerless = frames(sizes == pix);
    S.other      = frames(sizes > pix & sizes ~= normal_size);
    fprintf('CoreView %d: %d files, normal size %d bytes\n', coreview_num, numel(L), normal_size);
    fprintf('  empty (< %d bytes): %d  %s\n', pix, numel(S.empty), ranges(S.empty));
    fprintf('  headerless (== %d bytes): %d  %s\n', pix, numel(S.headerless), ranges(S.headerless));
    fprintf('  other sizes: %d  %s\n', numel(S.other), ranges(S.other));
end

function s = ranges(v)
    % 3 5 6 7 -> '3, 5-7'
    if isempty(v), s = ''; return; end
    v = sort(v(:)');
    brk = [0 find(diff(v) > 1) numel(v)];
    parts = cell(1, numel(brk) - 1);
    for k = 1:numel(brk) - 1
        a = v(brk(k) + 1);  b = v(brk(k + 1));
        if a == b, parts{k} = sprintf('%d', a); else, parts{k} = sprintf('%d-%d', a, b); end
    end
    s = strjoin(parts, ', ');
end
