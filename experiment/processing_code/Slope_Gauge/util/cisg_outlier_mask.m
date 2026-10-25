function [out, info] = cisg_outlier_mask(Sx, Sy, prev, next, opts)
%CISG_OUTLIER_MASK Space-time normalised median test of one slope frame (PIV-style validation).
%
%   [out, info] = cisg_outlier_mask(Sx, Sy, prev, next, opts)
%
% Sx, Sy     - ROI slopes of the frame tested (NaN = masked)
% prev, next - structs with Sx, Sy of frames n-1 and n+1 ([] = not available)
% out        - logical, size(Sx): pixels flagged as outliers (whole blocks)
%
% The slopes are block-averaged on a q x q grid (the colour smoothing makes neighbouring
% pixels nearly equal, so a pixel-level test would see nothing).  For each block u0, the
% neighbours are the w x w blocks around it in the same frame plus the same block in frames
% n-1 and n+1, and (Westerweel & Scarano 2005)
%
%   r = |u0 - med(u_i)| / (med|u_i - med(u_i)| + eps0),   outlier if r > thr (for Sx or Sy).
%
% Blocks with fewer than min_nb valid neighbours are not tested.  Flagged pixels should be
% masked, not filled.  Large failed patches (wider than w blocks) are not seen by a median test.
%
% opts fields (defaults): q 8 (binned px), w 5, thr 3, eps0 0.01 (slope), min_nb 8

    if nargin < 5, opts = struct(); end
    d = struct('q', 8, 'w', 5, 'thr', 3, 'eps0', 0.01, 'min_nb', 8);
    for f = fieldnames(d).', if ~isfield(opts, f{1}), opts.(f{1}) = d.(f{1}); end, end
    q = opts.q;  h = floor(opts.w / 2);
    [nr0, nc0] = size(Sx);  nr = floor(nr0 / q);  nc = floor(nc0 / q);
    blk = @(A) squeeze(mean(reshape(double(A(1:nr*q, 1:nc*q)), q, nr, q, nc), [1 3], 'omitnan'));

    bad = false(nr, nc);  rmax = zeros(nr, nc);
    for comp = {'Sx', 'Sy'}
        if strcmp(comp{1}, 'Sx'), A = blk(Sx); else, A = blk(Sy); end
        P = nan(nr + 2*h, nc + 2*h);  P(h+1:h+nr, h+1:h+nc) = A;
        N = nan(nr, nc, (2*h + 1)^2 + 1);  k = 0;
        for di = -h:h
            for dj = -h:h
                if di == 0 && dj == 0, continue; end
                k = k + 1;  N(:, :, k) = P(h+1+di:h+nr+di, h+1+dj:h+nc+dj);
            end
        end
        for T = {prev, next}
            if ~isempty(T{1}), k = k + 1;  N(:, :, k) = blk(T{1}.(comp{1})); end
        end
        N = N(:, :, 1:k);
        med = median(N, 3, 'omitnan');
        rm  = median(abs(N - med), 3, 'omitnan');
        r = abs(A - med) ./ (rm + opts.eps0);
        tested = isfinite(A) & sum(isfinite(N), 3) >= opts.min_nb;
        r(~tested) = 0;
        bad = bad | r > opts.thr;
        rmax = max(rmax, r);
    end
    out = false(nr0, nc0);
    out(1:nr*q, 1:nc*q) = repelem(bad, q, q);
    out = out & (isfinite(Sx) | isfinite(Sy));
    info = struct('block_fraction', mean(bad(:)), 'pixel_fraction', mean(out(:)), 'r_max', rmax, 'opts', opts, ...
                  'temporal', [~isempty(prev) ~isempty(next)]);
end
