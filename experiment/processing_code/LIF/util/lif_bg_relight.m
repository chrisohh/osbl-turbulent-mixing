function [bg, info] = lif_bg_relight(bg, fname, avail, bg_frames, downsample, z_cm, opts)
%LIF_BG_RELIGHT Re-light the bg reference when the lighting steps after it.
%
%   [bg, info] = LIF_BG_RELIGHT(bg, fname, avail, bg_frames, downsample, z_cm, opts)
%
% bg (from lif_load_bg) is the theta reference.  lif_z99 removes drift
% between bg and a frame with ONE offset per column (or per frame), which
% is right for lamp / exposure drift (95: 1.5%, uniform).  A lighting STEP
% is not uniform: CoreView 102 halves its counts at t ~ 14.5 s (frames
% ~210-250, after its bg frames 1-91), 0.59 x in the air, 0.58 x at mid
% depth, 0.54 x at the bottom.  The 8% spread is as large as the faint dye
% signal (theta 0.05-0.1): the per-column offset taken from the bottom made
% the mid depths read theta ~ -0.07 and hid the dye there.
%
% Scan: median counts of the lit air (opts.band_cm(1) to opts.air_max_cm
% above the bg water level) in every opts.scan_step-th frame.
% If any frame differs from the bg frames by more than opts.tol, the first
% lighting step is located, opts.n_frames frames
% from opts.settle_s after it are averaged (P) and
%       bg <- bg .* G,   G = smooth(P ./ bg)
% G is smoothed by normalised convolution (opts.sigma_cm) over the rows
% away from the surface (opts.band_cm = [above below] the bg water level:
% the surface line, its drop after the bg frames -- 102: 2.6 cm -- and any
% surface dye) and the dark lid (bg < dark_frac x median), then filled
% across them.  The air above must stay in: lif_z99 takes its drift
% reference there.
% bg keeps its own dye-free image; only its lighting is changed.
%
% opts (defaults): scan_step 25, tol 0.03, settle_s 0.5, n_frames 10,
% frame_step 5, sigma_cm 0.25 (light: P is 10 frames, ~1% shot noise; the
% air strip under the lid is only ~4 cm), fill_sigma_cm 2 (across band / lid), band_cm [1 4], dark_frac 0.3, air_max_cm 6, fs 50, z_water 0 (bg water level, cm)
%
% info: .frames (scanned), .rel (air counts / bg), .step_frame, .gain (median G),
%       .gain_range ([5 95] percentiles of G), .applied

    def = struct('scan_step', 25, 'tol', 0.03, 'settle_s', 0.5, 'n_frames', 10, ...
                 'frame_step', 5, 'sigma_cm', 0.25, 'fill_sigma_cm', 2, 'band_cm', [1 4], 'dark_frac', 0.3, 'air_max_cm', 6, 'fs', 50, 'z_water', 0);
    if nargin < 7, opts = struct(); end
    f = fieldnames(def);
    for k = 1:numel(f), if ~isfield(opts, f{k}), opts.(f{k}) = def.(f{k}); end, end

    info = struct('applied', false, 'step_frame', NaN, 'gain', 1, 'gain_range', [1 1]);
    % judged on the lit AIR above the water only: it is never dyed, while the
    % whole-image median also falls 4-7% late in every run as the dye fills
    % the water (95, 118, 120) -- that is dye, not lighting
    air = z_cm(:) > opts.z_water + opts.band_cm(1) & z_cm(:) < opts.z_water + opts.air_max_cm;
    A   = double(bg(air, :));
    lit = A > opts.dark_frac * median(double(bg(:)));
    if nnz(lit) < 1000, warning('lif_bg_relight:air', 'No lit air above the water: lighting not checked.'); return; end
    b0  = median(A(lit));
    rows = find(air) * downsample - (downsample - 1);      % full-resolution rows
    sc = avail(1:opts.scan_step:end);
    rel = nan(size(sc));
    for k = 1:numel(sc)
        I = lif_load_raw(fname(sc(k)));
        Ia = double(I(rows, 1:downsample:end));
        rel(k) = median(Ia(lit)) / b0;
    end
    info.frames = sc;  info.rel = rel;
    if all(abs(rel - 1) <= opts.tol), return; end

    % the step: last scanned frame still matching bg -> first one that does not,
    % then take frames once the new level has settled (it may ramp over a second)
    k1 = find(abs(rel - 1) > opts.tol, 1, 'first');
    ns = round(opts.settle_s * opts.fs);
    lvl = rel(k1:min(k1 + 3, end));
    k2 = k1 - 1 + find(abs(rel(k1:end) - median(lvl)) < opts.tol/3, 1, 'first');
    if isempty(k2), k2 = k1; end
    f0 = sc(k2) + ns;
    pf = intersect(f0 : opts.frame_step : f0 + opts.frame_step*(opts.n_frames - 1), avail);
    pf = setdiff(pf, bg_frames);
    if isempty(pf), warning('lif_bg_relight:frames', 'No frames after the lighting step.'); return; end
    P = lif_load_bg(fname, pf, downsample);

    bg  = double(bg);  P = double(P);                       % lif_load_bg can return uint16
    G   = P ./ max(bg, 1);
    dz  = abs(median(diff(z_cm)));
    sig = opts.sigma_cm / dz;                               % rows (= columns: square pixels)
    zc  = z_cm(:);
    M   = double(zc > opts.z_water + opts.band_cm(1) | zc < opts.z_water - opts.band_cm(2)) .* ones(1, size(bg, 2));
    M(~isfinite(G) | G <= 0 | bg < opts.dark_frac * median(bg(:))) = 0;   % lid / unlit: not lighting
    G(M == 0) = 0;
    num = imgaussfilt(G, sig);  den = imgaussfilt(M, sig);
    Gs  = num ./ max(den, eps);
    % the band and the lid: a much wider kernel reaches across them
    sw   = opts.fill_sigma_cm / dz;
    numW = imgaussfilt(G, sw);  denW = imgaussfilt(M, sw);
    fill = den < 0.5;
    Gs(fill) = numW(fill) ./ max(denW(fill), eps);

    bg = bg .* Gs;
    info.applied = true;  info.step_frame = sc(k1);  info.frames_used = pf;
    info.gain = median(Gs(:));  info.gain_range = prctile(Gs(:), [5 95]);
end
