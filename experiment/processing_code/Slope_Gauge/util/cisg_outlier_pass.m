function S = cisg_outlier_pass(cache_dir, frames, opts)
%CISG_OUTLIER_PASS Run cisg_outlier_mask over a slope cache and store the masks in the cache files.
%
%   S = cisg_outlier_pass(cache_dir)              % every frame_NNNN.mat in the folder
%   S = cisg_outlier_pass(cache_dir, frames, opts)
%
% Each frame is tested against its neighbours in space and in frames n-1, n+1 (when cached).
% The mask is saved in the frame file as 'out' (sparse logical, with 'out_opts'); Sx, Sy
% are NOT changed, and cisg_cache_frame applies 'out' when it loads the frame.  Re-running
% replaces the masks (tests always use the raw Sx, Sy).  opts: see cisg_outlier_mask.
% S: frames, pixel_fraction (per frame).

    if nargin < 3, opts = struct(); end
    if nargin < 2 || isempty(frames)
        L = dir(fullfile(cache_dir, 'frame_*.mat'));
        frames = sort(cellfun(@(s) sscanf(s, 'frame_%d.mat'), {L.name}));
    end
    fn = @(n) fullfile(cache_dir, sprintf('frame_%04d.mat', n));
    buf = containers.Map('KeyType', 'double', 'ValueType', 'any');   % raw slopes of n-1, n, n+1
    frac = nan(numel(frames), 1);
    t0 = tic;
    for i = 1:numel(frames)
        n = frames(i);
        for m = [n - 1, n, n + 1]                                     % fill the rolling buffer
            if ~isKey(buf, m) && isfile(fn(m)), buf(m) = load(fn(m), 'Sx', 'Sy'); end
        end
        for m = cell2mat(keys(buf)), if m < n - 1, remove(buf, m); end, end
        if ~isKey(buf, n), continue; end
        cur = buf(n);  pv = [];  nx = [];
        if isKey(buf, n - 1), pv = buf(n - 1); end
        if isKey(buf, n + 1), nx = buf(n + 1); end
        [o, oi] = cisg_outlier_mask(cur.Sx, cur.Sy, pv, nx, opts);
        out = sparse(o);  out_opts = oi.opts;  %#ok<NASGU>
        save(fn(n), 'out', 'out_opts', '-append');
        frac(i) = oi.pixel_fraction;
        if mod(i, 100) == 0 || i == numel(frames)
            fprintf('outlier pass: %d/%d frames, median %.2f%% of ROI flagged, %.0f s\n', i, numel(frames), ...
                    100 * median(frac(1:i), 'omitnan'), toc(t0));
        end
    end
    S = struct('frames', frames, 'pixel_fraction', frac);
end
