function zs = lif_surface_edge(Is, z_cm, lo, hi, dx_mm, opts)
%LIF_SURFACE_EDGE Water surface per column, as the side view shows it.
%
%   zs = LIF_SURFACE_EDGE(Is, z_cm, lo, hi, dx_mm, opts)
%
% Is     - [ny x nx] (filtered) counts, row 1 at the top (largest z)
% z_cm   - [1 x ny] row heights (cm)
% lo, hi - search window per column (cm), scalars or [1 x nx]
% dx_mm  - pixel size (mm)
%
% The side view looks THROUGH the tank width, so the surface is not a thin
% line but a dark band: the surface everywhere across the width, crests and
% troughs, seen together.  opts.method:
%   'edge'    - TOP edge of the band: going down from the window top, the
%               first level where the counts fall halfway from the backlit
%               air to the band minimum (sub-pixel).  = the highest water
%               anywhere across the width, i.e. the surface the view shows.
%               A sharp edge, stable column to column without a median.
%   'darkest' - darkest row in the window (old): can sit anywhere inside the
%               band and jumps between its top and bottom.
% Columns where the air above the window is not clearly brighter than the
% band (outside the lit panel) return NaN: opts.contrast = minimum
% (I_air - I_band) / I_air, air measured opts.above_mm above the window.
%
% opts (defaults): method 'edge', above_mm 8, contrast 0.15

    def = struct('method', 'edge', 'above_mm', 8, 'contrast', 0.15);
    if nargin < 6, opts = struct(); end
    for f = fieldnames(def).'
        if ~isfield(opts, f{1}), opts.(f{1}) = def.(f{1}); end
    end

    nx = size(Is, 2);
    z  = z_cm(:);
    if isscalar(lo), lo = repmat(lo, 1, nx); end
    if isscalar(hi), hi = repmat(hi, 1, nx); end
    na = max(1, round(opts.above_mm / dx_mm));
    zs = nan(1, nx);

    for j = 1:nx
        if isnan(lo(j)) || isnan(hi(j)), continue; end
        r = find(z >= lo(j) & z <= hi(j));             % top row first
        if numel(r) < 3 || r(1) < 2, continue; end
        [Imin, k] = min(Is(r, j));
        Ia = median(Is(max(1, r(1) - na) : r(1) - 1, j));
        if (Ia - Imin) / max(Ia, 1) < opts.contrast, continue; end
        if strcmp(opts.method, 'darkest')
            zs(j) = z(r(k));  continue
        end
        thr = (Ia + Imin) / 2;
        seg = Is(r(1) : r(k), j);
        q   = find(seg < thr, 1, 'first');
        rr  = r(1) + q - 1;
        if rr > r(1)                                   % sub-pixel crossing
            f = (Is(rr-1, j) - thr) / max(Is(rr-1, j) - Is(rr, j), eps);
            zs(j) = z(rr-1) + f * (z(rr) - z(rr-1));
        else
            zs(j) = z(rr);
        end
    end
end
