function [z99_x, jc, z99, zeta, theta_bar, theta_img, T_img, stats] = lif_z99(I, z_cm, z_surf, z_front, dx_mm, opts)
%LIF_Z99 Depth above which 99% of the dye resides (paper eq. 3.14), vs x.
%
%   [z99_x, jc, z99, zeta, theta_bar, theta_img, T_img, stats] = LIF_Z99(I, z_cm, z_surf, z_front, dx_mm, opts)
%
%   int_{z99}^{0} theta dz = frac * int_{-H}^{0} theta dz
%
% with z = 0 at the (local) surface and -H the bottom of the camera frame
% (the dye stays well inside it).  The paper averages theta horizontally
% first, but its LIF plane was TRANSVERSE (y, z), where the dye is roughly
% uniform in y.  This view is LONGITUDINAL (x, z): the dye evolves along
% the tank, so z99 is computed locally, on theta averaged over x bins of
% width opts.bin_cm  ->  z99_x.  The frame-wide value (theta averaged over
% all x, as in the paper) is still returned as z99.
%
% Concentration: the dye absorbs, so theta = ln(I_clear / I) (Beer-Lambert,
% arbitrary units -- z99 only needs the shape).  There is no dye-free
% background frame, and the clear-water brightness falls with depth from
% the lighting alone (~600 -> ~500 counts over the frame), which an
% integral to 99% would read as deep dye.  So, unless opts.bg is given,
% I_clear(z) is fitted PER COLUMN to ln(I) over the clear water below the
% dye front (z < z_front - fit_margin_cm, polynomial of order fit_order)
% and extrapolated up through the dye layer.  Negative theta (noise) is
% kept so it averages out rather than biasing the tail.
%
% Each column is referenced to its own surface (zeta = z - z_surf) before
% the horizontal average, so surface waves don't smear the profile.
%
% Inputs:
%   I        - [ny x nx] counts (row 1 = top)
%   z_cm     - [1 x ny] row heights (cm)
%   z_surf   - [1 x nx] surface height per column (cm), from lif_dye_front
%   z_front  - [1 x nx] dye-front height per column (cm), from lif_dye_front
%   dx_mm    - pixel size of I (mm)
% opts fields (defaults):
%   frac 0.99, smooth_mm 1.5, surf_gap_cm 0.4, fit_order 2,
%   fit_margin_cm 0.5, bg [] (dye-free image, same size as I: replaces the fit),
%   bg_normalize true (with bg: zero theta on the clear water below the
%   front in each column, removing bg-vs-frame brightness drift)
%
%   bin_cm 1 (x bin width for z99_x)
%   dark [] (dark frame, same size as I: subtracted from I and bg before the
%   ratio -- any camera offset biases I/I_bg, most in the dark lower image)
%   bg_shift true (with bg: move each bg column so its surface sits on this
%   frame's -- the flat bg otherwise puts air under wave crests),
%   bg_search_cm 1 (bg surface search, +/- around this frame's surface),
%   surf_median_cm 0.6 (running median of the bg surface along x)
%   lit_frac 0.5 (columns whose bg water is darker than this fraction of the
%   frame median are unlit: blank, out of z99)
%   trough_cm 20 (window for the trough envelope; the band above the local
%   trough level is blank and out of the integrals -- look-through mixing)
%   drift_ref 'auto' ('auto' = 'clear' while >= clear_min_cm (6) of clear
%   water remains below the front, else 'air'; 'air' = backlit air above the water, never dyed;
%   'clear' = clear water below the dye front), air_window_cm [2 6],
%   air_above_cm 1 (air rows must be this far above the surface)
%   bottom_cm 2 (dye still detectable over the bottom bottom_cm of the frame
%   -> z99 deeper than the frame: NaN, flagged in stats.deep_x)
%   snr_min 4 (NaN z99 unless theta somewhere above the front exceeds this
%   many robust sigma of the clear-water noise: no dye -> no z99, rather
%   than a z99 at the bottom of the frame)
%   tail_cm 1 (the integral stops at the first depth below which theta stays
%   inside the noise for this many cm, so the clear-water tail cannot own
%   the last 1% -- see depth_stats)
%
% Outputs:
%   z99_x     - [1 x nbin] local z99 (cm, NEGATIVE = below the surface;
%               NaN where there is no detectable dye)
%   jc        - [1 x nbin] centre column of each bin (index into x / z_surf)
%   z99       - frame-wide z99 (cm): opts.z99_all 'median' (default) = median
%               of z99_x; 'profile' = from the all-x averaged theta (paper)
%   stats     - struct: z90 (90% depth), zrms (RMS dye depth, -sqrt(<z^2>)),
%               z90_x, zrms_x (per bin), peak (max of theta_bar).  All depths
%               NaN while peak < opts.peak_min (0.03): too little dye to measure
%   zeta      - [nz x 1] depth grid (cm, 0 = surface, negative down)
%   theta_bar - [nz x 1] horizontally averaged theta on zeta
%   theta_img - [ny x nx] theta in image coordinates, for display (NaN above
%               the surface; not zeroed in the meniscus gap)
%   T_img     - [ny x nx] transmission I/I_bg, UNSMOOTHED, with the same
%               per-column drift normalisation (1 = clear water, < 1 = dye);
%               exp(-theta_img) when there is no bg

    def = struct('frac', 0.99, 'smooth_mm', 1.5, 'surf_gap_cm', 0.4, ...
                 'fit_order', 2, 'fit_margin_cm', 0.5, 'bg', [], 'bg_normalize', true, 'bin_cm', 1, 'step_cm', 0.25, 'dark', [], ...
                 'bg_shift', true, 'bg_search_cm', 1, 'surf_median_cm', 0.6, ...
                 'snr_min', 4, 'tail_cm', 1, 'tail_snr', 2, 'tail_rel', 0.02, 'z99_all', 'median', 'peak_min', 0.03, ...
                 'lit_frac', 0.5, 'trough_cm', 20, 'band_max_cm', 3, 'drift_ref', 'auto', 'clear_min_cm', 6, 'air_delta', 0, 'delta_tol', 0.01, ...
                 'air_window_cm', [2 6], 'air_above_cm', 1, 'bottom_cm', 2);
    if nargin < 6, opts = struct(); end
    for f = fieldnames(def).'
        if ~isfield(opts, f{1}), opts.(f{1}) = def.(f{1}); end
    end

    if ~isempty(opts.dark)
        I = I - opts.dark;
        if ~isempty(opts.bg), opts.bg = opts.bg - opts.dark; end
    end
    nx  = size(I, 2);
    sig = opts.smooth_mm / dx_mm;
    lnI = log(max(imgaussfilt(I, sig), 1));
    if ~isempty(opts.bg)
        lnB = log(max(imgaussfilt(opts.bg, sig), 1));
    end

    % Background surface alignment.  The bg frame's surface is flat; where
    % this frame's surface is higher (crests), the bg pixels just under it
    % are air / the dark surface line, I > I_bg and theta goes negative
    % ("undyed" white band at the top).  Find the bg surface per column the
    % same way as lif_dye_front (darkest smoothed row, here within
    % bg_search_cm of this frame's surface), and shift each bg column so
    % its surface sits on this frame's.  The shift is a few mm, so the
    % lighting gradient it also moves is negligible (and normalised below).
    z = z_cm(:);
    bg_shift = zeros(1, nx);
    if ~isempty(opts.bg) && opts.bg_shift
        zb = nan(1, nx);
        for j = 1:nx
            r = find(abs(z - z_surf(j)) <= opts.bg_search_cm);
            if isempty(r), continue; end
            [~, k] = min(lnB(r, j));
            zb(j) = z(r(k));
        end
        zb = movmedian(zb, max(1, round(10*opts.surf_median_cm / dx_mm)), 'omitnan');
        bg_shift = z_surf - zb;
        bg_shift(isnan(bg_shift)) = 0;
    end

    dz   = abs(z_cm(2) - z_cm(1));
    zeta = (0 : -dz : (min(z_cm) - max(z_surf))).';
    TH   = nan(numel(zeta), nx);
    theta_img = nan(size(I));
    T_img     = nan(size(I));

    % Dye front for choosing each column's clear water, filled and smoothed
    % along x: a single column whose front detection failed (NaN, or wildly
    % off) would otherwise take dye as "clear water", get over-corrected and
    % show as a full-height white streak in theta.
    wx = max(1, round(10*opts.surf_median_cm / dx_mm));
    zf = movmedian(z_front, 4*wx, 'omitnan');
    if any(~isnan(zf)), zf = fillmissing(zf, 'nearest'); end

    % Lit columns.  The far edge of the frame (CoreView 95 Water_SURF:
    % x > ~23 cm) is outside the backlight, nearly black in water and air;
    % theta there is the log-ratio of two dark numbers and drew z99(x) up to
    % the surface.  Judged on the bg frame: water (and, for the drift
    % reference, air) at least opts.lit_frac of the frame median.
    lit = true(1, nx);  air_lit = false(1, nx);
    if ~isempty(opts.bg)
        wr  = z < min(z_surf) - 1;
        Bw  = median(opts.bg(wr, :), 1);
        lit = Bw >= opts.lit_frac * median(Bw);
        arw = z >= opts.air_window_cm(1) & z <= opts.air_window_cm(2) & z > max(z_surf) + opts.air_above_cm;
        if nnz(arw) >= 20
            Ba = median(opts.bg(arw, :), 1);
            air_lit = Ba >= opts.lit_frac * median(Ba) & lit;
        end
    end

    % Pass 1: theta (and transmission) per column, plus the per-column
    % drift offset measured on the clear water below the front and on the
    % backlit air above the water
    off_th = nan(1, nx);  off_T = nan(1, nx);
    off_air_th = nan(1, nx);  off_air_T = nan(1, nx);
    for j = 1:nx
        if isnan(z_surf(j)), continue; end
        if isnan(zf(j)), clear_rows = z < z_surf(j) - opts.surf_gap_cm;
        else,            clear_rows = z < zf(j) - opts.fit_margin_cm; end
        if isempty(opts.bg)
            if nnz(clear_rows) < 10*(opts.fit_order + 1), continue; end
            zc = z(clear_rows);
            [p, ~, mu] = polyfit(zc, lnI(clear_rows, j), opts.fit_order);
            lnC = polyval(p, z, [], mu);
        else
            lnC = lnB(:, j);
            Bj  = opts.bg(:, j);
            if bg_shift(j) ~= 0                        % bg surface onto this frame's
                lnC = fillmissing(interp1(z, lnC, z - bg_shift(j), 'linear'), 'nearest');
                Bj  = fillmissing(interp1(z, Bj,  z - bg_shift(j), 'linear'), 'nearest');
            end
        end
        th = lnC - lnI(:, j);
        if isempty(opts.bg), Tc = exp(-th); else, Tc = I(:, j) ./ max(Bj, 1); end
        theta_img(:, j) = th;  T_img(:, j) = Tc;
        if ~isempty(opts.bg) && opts.bg_normalize && nnz(clear_rows) >= 20
            off_th(j) = median(th(clear_rows));
            off_T(j)  = median(Tc(clear_rows));
        end
        if ~isempty(opts.bg) && air_lit(j)
            ar = z >= opts.air_window_cm(1) & z <= opts.air_window_cm(2) ...
                 & z > z_surf(j) + opts.air_above_cm;
            if nnz(ar) >= 20
                off_air_th(j) = median(th(ar));
                off_air_T(j)  = median(Tc(ar));
            end
        end
    end

    % Lighting/exposure drift between bg and frame (CoreView 95: frame 1875
    % is ~1.5% darker than frame 1 even in deep clear water) adds a constant
    % theta at every depth, which the 99% integral turns into "dye" down to
    % the frame bottom.  Reference (opts.drift_ref):
    %   'air'   - the backlit air above the water: same light, same exposure,
    %             never dyed.  Default: the clear-water reference fails once
    %             the dye fills the frame (it then subtracts real dye).
    %   'clear' - the clear water below the dye front (old behaviour).
    % Either way the per-column offsets are median-filtered along x before
    % being removed, so one bad column can no longer whiten itself.
    % 'auto': clear water while there is plenty of it below the dye front
    % (median thickness >= clear_min_cm), the air once the dye fills the
    % frame.  The two agree within ~0.4 cm in z99 mid-run (CoreView 95,
    % 47-66 s), but early on the air reference leaves a small positive theta
    % in the water (z99 -3.7 vs -2.0 cm at 34 s), and late the clear-water
    % one subtracts real dye (z99 -13.7 vs -16.0 cm at 70 s).
    %
    % Water-minus-air offset difference: constant (delta0, a small genuine
    % water/air difference) while the reference water is truly clear, grows
    % once faint dye reaches it.  With opts.air_delta = delta0 (learned by
    % the caller from earlier frames), 'auto' switches to air + delta0 as
    % soon as delta exceeds delta0 by opts.delta_tol.  Without it, 'auto'
    % falls back to the clear-water thickness test.  CoreView 95: delta stays
    % within +/-0.005 from 18 to 66 s, then jumps to 0.022-0.036 at 70-74 s
    % as faint dye reaches the reference water -> air_delta 0, delta_tol 0.01.
    delta = median(off_th - off_air_th, 'omitnan');
    ref = opts.drift_ref;
    if strcmp(ref, 'auto')
        if ~isempty(opts.air_delta) && ~isnan(delta)
            if delta - opts.air_delta > opts.delta_tol, ref = 'air'; else, ref = 'clear'; end
        else
            clear_thk = median(zf - opts.fit_margin_cm - min(z), 'omitnan');
            if isnan(clear_thk) || clear_thk >= opts.clear_min_cm, ref = 'clear'; else, ref = 'air'; end
        end
    end
    d0 = 0;
    if strcmp(ref, 'air') && ~isempty(opts.air_delta), d0 = opts.air_delta; end
    drift_used = 'none';
    if ~isempty(opts.bg) && opts.bg_normalize
        if strcmp(ref, 'air') && nnz(~isnan(off_air_th)) >= 4*wx
            off_th = off_air_th + d0;  off_T = off_air_T * exp(-d0);  drift_used = 'air';
        elseif any(~isnan(off_th))
            drift_used = 'clear';
        end
    end
    if ~strcmp(drift_used, 'none')
        off_th = fillmissing(movmedian(off_th, 4*wx, 'omitnan'), 'nearest');
        off_T  = fillmissing(movmedian(off_T,  4*wx, 'omitnan'), 'nearest');
        theta_img = theta_img - off_th;
        T_img     = T_img ./ off_T;
    end

    % Trough envelope: lowest surface within +/- trough_cm/2.  The side view
    % looks THROUGH the tank width, so between the troughs and the crests
    % the line of sight mixes water near the front wall with air / the
    % surface line further back -- no background matches it ("white peaks"
    % at the crests).  That band is left out of the depth integrals; depths
    % stay referenced to the local surface.
    z_tr = movmin(z_surf, max(1, round(10*opts.trough_cm / dx_mm)));

    % For DISPLAY the band is drawn as a uniform layer under the local
    % surface, of thickness band_cm = crest-to-trough range of the surface
    % (5th-95th percentile over the lit x), clipped to [surf_gap_cm,
    % band_max_cm].  The trough envelope above changes shape every frame and
    % made a blanked patch that flickered in the video; one thickness per
    % frame grows smoothly with the waves.  The images themselves are only
    % blanked above the surface and in unlit columns; lif_surf_figure covers
    % the band with a grey layer.
    zl = z_surf(lit & ~isnan(z_surf));
    if isempty(zl), band_cm = opts.surf_gap_cm;
    else, band_cm = min(max(prctile(zl, 95) - prctile(zl, 5), opts.surf_gap_cm), opts.band_max_cm); end

    % Pass 2: mask unlit columns and above the surface, depth-referenced
    % profiles for z99
    for j = 1:nx
        if isnan(z_surf(j)) || ~lit(j)
            theta_img(:, j) = NaN;  T_img(:, j) = NaN;  continue
        end
        above = z > z_surf(j);
        theta_img(above, j) = NaN;  T_img(above, j) = NaN;
        band = z > min(z_tr(j), z_surf(j) - opts.surf_gap_cm);
        th = theta_img(:, j);
        th(band) = 0;                                  % surface band / meniscus
        TH(:, j) = interp1(z - z_surf(j), th, zeta, 'linear', NaN);
    end
    % Depth below which the water is clear, i.e. where theta is pure noise:
    % the dye front is the deepest dye the image shows.  NaN front = the
    % column never left the clear-water threshold, so there is no dye at all.
    depth_front = z_front - z_surf;                    % [1 x nx], <= 0, NaN = clear

    % Local depths in a window bin_cm wide that slides along x in step_cm
    % steps (overlapping), each from that window's own noise.  Adjacent
    % non-overlapping bins smoothed the line over the bin width; the
    % sliding window resolves features down to ~bin_cm.  A 3-point median
    % removes single-window outliers.
    nb    = max(1, round(10*opts.bin_cm / dx_mm));     % columns per window
    ns    = max(1, round(10*opts.step_cm / dx_mm));    % columns per step
    jc    = ceil(nb/2) : ns : nx - floor(nb/2);
    nbin  = numel(jc);
    z99_x = nan(1, nbin);  z90_x = nan(1, nbin);  zrms_x = nan(1, nbin);  deep_x = false(1, nbin);
    for b = 1:nbin
        cols  = max(1, jc(b) - floor(nb/2)) : min(nx, jc(b) + ceil(nb/2) - 1);
        s = depth_stats(zeta, mean(TH(:, cols), 2, 'omitnan'), ...
                        median(depth_front(cols), 'omitnan'), opts);
        z99_x(b) = s.z99;  z90_x(b) = s.z90;  zrms_x(b) = s.zrms;  deep_x(b) = s.deep;
    end
    if nbin >= 3
        z99_x  = movmedian(z99_x,  3, 'omitnan', 'Endpoints', 'shrink');
        z90_x  = movmedian(z90_x,  3, 'omitnan', 'Endpoints', 'shrink');
        zrms_x = movmedian(zrms_x, 3, 'omitnan', 'Endpoints', 'shrink');
    end
    z99_x(~lit(jc)) = NaN;  z90_x(~lit(jc)) = NaN;  zrms_x(~lit(jc)) = NaN;

    theta_bar = mean(TH, 2, 'omitnan');
    ok = ~isnan(theta_bar);
    peak = max(theta_bar(ok));

    % Not enough dye to measure a depth: the all-x profile peaks below
    % peak_min (CoreView 95: theta_bar ~0.01 at t ~ 17 s, i.e. ~1%
    % absorption -- the size of the background residuals).  No line beats a
    % noise-driven one.
    if ~(peak >= opts.peak_min)
        z99_x(:) = NaN;  z90_x(:) = NaN;  zrms_x(:) = NaN;
    end

    switch opts.z99_all
        case 'median'
            % Frame value = median of the local values.  The all-x profile
            % (paper's way) has tiny random noise, so while the layer is thin
            % systematic residuals at the 1% level own the tail: frame 321
            % gave -3.9 cm from the profile vs -0.9 cm median of the bins;
            % from t ~ 30 s on the two agree within 0.2 cm.
            z99 = median(z99_x, 'omitnan');  z90 = median(z90_x, 'omitnan');
            % Most bins have dye down to the frame bottom: the median of the
            % rest would be biased shallow, so the frame z99 is undefined too
            if mean(deep_x(~isnan(z90_x))) > 0.5, z99 = NaN; end
            zrms = median(zrms_x, 'omitnan');
        case 'profile'
            s = depth_stats(zeta(ok), theta_bar(ok), median(depth_front, 'omitnan'), opts);
            z99 = s.z99;  z90 = s.z90;  zrms = s.zrms;
            if ~(peak >= opts.peak_min), z99 = NaN; z90 = NaN; zrms = NaN; end
    end
    stats = struct('z90', z90, 'zrms', zrms, 'z90_x', z90_x, 'zrms_x', zrms_x, 'peak', peak, ...
                   'deep_x', deep_x, 'deep_frac', mean(deep_x), 'drift_ref', drift_used, ...
                   'lit_frac', mean(lit), 'drift_delta', delta, ...
                   'lit', lit, 'band_cm', band_cm);
    zeta = zeta(ok);  theta_bar = theta_bar(ok);
end
function s = depth_stats(zeta, th, z_front, opts)
    % Level above which opts.frac of the integral of th lies (zeta <= 0, top
    % first).  z_front = dye-front depth for this profile (<= 0, NaN = none).
    %
    % The 99% criterion is brutally sensitive to the deep tail: the last 1%
    % of the integral is the same size as the noise collected over a few cm
    % of clear water.  Taken literally it over-predicts in exactly the two
    % cases below, so both are handled before the integral is formed.
    %   no dye        - C(end) is a zero-mean random walk that only reaches
    %                   its final value at the last row, so the 99% crossing
    %                   lands at the bottom of the frame.  C(end) > 0 half
    %                   the time, so the old `C(end) <= 0` guard caught only
    %                   half of these.
    %   thin layer    - the layer contributes A near the surface, the noise
    %                   below contributes ~sigma*sqrt(H*dz) with either sign,
    %                   and once that is ~A/100 the crossing is set by the
    %                   noise rather than by the dye.
    s  = struct('z99', NaN, 'z90', NaN, 'zrms', NaN, 'deep', false);
    ok = ~isnan(th);
    if nnz(ok) < 2, return; end
    zeta = zeta(ok);  th = th(ok);
    dz   = abs(zeta(2) - zeta(1));

    % Noise floor of theta, from the clear water below the front
    if isnan(z_front), z_front = 0; end
    z_clear = min(z_front - opts.fit_margin_cm, 0.75 * min(zeta));
    nz = zeta <= z_clear;
    if nnz(nz) < 10, return; end
    sig = 1.4826 * median(abs(th(nz) - median(th(nz))));   % robust sigma
    if ~(sig > 0), return; end

    % No dye: theta never rises out of the noise -> z99 undefined, not deep
    if max(th(zeta > z_clear)) < opts.snr_min * sig, return; end

    % Integrate only over the depth where dye is still detectable: below the
    % dye peak, the first level where the MEAN of theta over the next
    % opts.tail_cm has dropped into the noise of such a mean.  This is what
    % stops the clear-water tail from owning the last 1%.  (A max-over-the-
    % window < sigma test never fires: of ~40 noisy samples per cm one almost
    % always exceeds sigma, so the integral ran to the frame bottom.)
    nrun  = max(2, round(opts.tail_cm / dz));
    ths   = movmean(th, [0 nrun-1]);                   % mean of th and the next tail_cm
    sigm  = 1.4826 * median(abs(ths(nz) - median(ths(nz))));   % noise of that mean
    if ~(sigm > 0), sigm = sig / sqrt(nrun); end
    [thp, kp] = max(th);                               % dye peak
    % Plus a floor relative to the peak: averaged over all x the random
    % noise is tiny, and small SYSTEMATIC residuals (bg mismatch, lighting)
    % then pass as dye and drag the all-x z99 deep while the layer is thin
    % (frame 321: all-x -3.9 cm vs local median -0.9 cm).
    quiet = ths < max(opts.tail_snr * sigm, opts.tail_rel * thp);
    quiet(1:kp) = false;
    k_end = find(quiet, 1, 'first');
    % Dye still detectable at the bottom of the frame: the integral is cut
    % by the field of view, so z99 is deeper than the camera sees -> NaN,
    % flagged.  z90 / z_rms are kept (lower bounds, far less tail-driven).
    nbot = max(2, round(opts.bottom_cm / dz));
    if isempty(k_end)
        k_end = numel(zeta);
        s.deep = mean(th(max(1, end-nbot+1):end)) > max(opts.tail_snr * sigm, opts.tail_rel * thp);
    end
    if k_end < 2, return; end
    zeta = zeta(1:k_end);  th = th(1:k_end);

    C = cumtrapz(-zeta, th);                           % integral from the surface down
    if C(end) <= 0, return; end
    s.z99 = zeta(find(C >= opts.frac * C(end), 1, 'first'));
    if s.deep, s.z99 = NaN; end
    s.z90 = zeta(find(C >= 0.90 * C(end), 1, 'first'));
    % RMS dye depth: every part of the profile weighted by its dye, so the
    % faint tail barely matters (negative = below the surface)
    s.zrms = -sqrt(max(trapz(-zeta, zeta.^2 .* th) / C(end), 0));
end
