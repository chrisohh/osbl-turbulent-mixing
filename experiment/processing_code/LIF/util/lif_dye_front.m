function [z_front, z_surf, thr_col] = lif_dye_front(I, z_cm, dx_mm, opts)
%LIF_DYE_FRONT Lower edge of a dark (absorbing) dye layer under the surface.
%
%   [z_front, z_surf, thr_col] = LIF_DYE_FRONT(I, z_cm, dx_mm, opts)
%
% I     - [ny x nx] counts, row 1 at the top (largest z)
% z_cm  - [1 x ny] height of each row (cm)
% dx_mm - pixel size of I (mm), for the smoothing / window lengths
%
% 1. Gaussian smooth (opts.smooth_mm) to kill speckle.
% 2. Surface = darkest smoothed row in opts.surf_window_cm, per column,
%    opts.surf_median_cm running median along x.
% 3. Each column is scanned DOWN from opts.surf_gap_cm under the surface
%    to the first clear pixel.  Column-wise rather than a connected-region
%    fill, so dark patches lower down (uneven lighting, vignetting) can't
%    leak in through a thin path.  Clear means
%      opts.dye_thr empty: Is >= I_dye + dye_frac*(I_clear - I_dye), with
%        I_dye   = median of opts.top_band_cm just under the gap (median, since
%                  that band can clip the dark surface line),
%        I_clear = 98th percentile of the column below (clear plateau).
%        Needed because the clear-water level varies across the frame.
%      opts.dye_thr given: Is >= dye_thr (fixed counts).
%
% opts fields (defaults): dye_frac 0.9, dye_thr [], smooth_mm 1.5,
%   surf_window_cm [4 10], surf_gap_cm 0.4, top_band_cm 0.8, surf_median_cm 0.6
%
% Outputs (all [1 x nx], cm / counts):
%   z_front - front height (NaN where the column is clear under the gap)
%   z_surf  - surface height
%   thr_col - threshold used in each column

    def = struct('dye_frac', 0.9, 'dye_thr', [], 'smooth_mm', 1.5, ...
                 'surf_window_cm', [4 10], 'surf_gap_cm', 0.4, ...
                 'top_band_cm', 0.8, 'surf_median_cm', 0.6);
    if nargin < 4, opts = struct(); end
    for f = fieldnames(def).'
        if ~isfield(opts, f{1}), opts.(f{1}) = def.(f{1}); end
    end

    [ny, nx] = size(I);
    Is = imgaussfilt(I, opts.smooth_mm / dx_mm);

    rows   = find(z_cm >= opts.surf_window_cm(1) & z_cm <= opts.surf_window_cm(2));
    [~, k] = min(Is(rows, :), [], 1);
    z_surf = movmedian(z_cm(rows(k)), round(10*opts.surf_median_cm / dx_mm));

    below   = z_cm(:) < z_surf - opts.surf_gap_cm;   % [ny x nx]
    n_top   = max(1, round(10*opts.top_band_cm / dx_mm));   % band for I_dye
    z_front = nan(1, nx);
    thr_col = nan(1, nx);
    for j = 1:nx
        r0 = find(below(:, j), 1, 'first');
        if isempty(r0), continue; end
        col = Is(r0:end, j);
        if isempty(opts.dye_thr)
            I_dye   = median(col(1:min(n_top, end)));
            I_clear = prctile(col, 98);
            thr_col(j) = I_dye + opts.dye_frac * (I_clear - I_dye);
        else
            thr_col(j) = opts.dye_thr;
        end
        if col(1) >= thr_col(j), continue; end
        r1 = find(col >= thr_col(j), 1, 'first');
        if isempty(r1), r1 = ny - r0 + 2; end        % dyed to the frame bottom
        z_front(j) = z_cm(r0 + r1 - 2);
    end
end
