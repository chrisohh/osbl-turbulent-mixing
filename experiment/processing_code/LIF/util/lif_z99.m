function [z99_x, jc, z99, zeta, theta_bar, theta_img, T_img] = lif_z99(I, z_cm, z_surf, z_front, dx_mm, opts)
%LIF_Z99 Depth above which 99% of the dye resides (paper eq. 3.14), vs x.
%
%   [z99_x, jc, z99, zeta, theta_bar, theta_img, T_img] = LIF_Z99(I, z_cm, z_surf, z_front, dx_mm, opts)
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
%
% Outputs:
%   z99_x     - [1 x nbin] local z99 (cm, NEGATIVE = below the surface)
%   jc        - [1 x nbin] centre column of each bin (index into x / z_surf)
%   z99       - frame-wide z99 (cm), all x averaged as in the paper
%   zeta      - [nz x 1] depth grid (cm, 0 = surface, negative down)
%   theta_bar - [nz x 1] horizontally averaged theta on zeta
%   theta_img - [ny x nx] theta in image coordinates, for display (NaN above
%               the surface; not zeroed in the meniscus gap)
%   T_img     - [ny x nx] transmission I/I_bg, UNSMOOTHED, with the same
%               per-column drift normalisation (1 = clear water, < 1 = dye);
%               exp(-theta_img) when there is no bg

    def = struct('frac', 0.99, 'smooth_mm', 1.5, 'surf_gap_cm', 0.4, ...
                 'fit_order', 2, 'fit_margin_cm', 0.5, 'bg', [], 'bg_normalize', true, 'bin_cm', 1, 'dark', []);
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

    dz   = abs(z_cm(2) - z_cm(1));
    zeta = (0 : -dz : (min(z_cm) - max(z_surf))).';
    TH   = nan(numel(zeta), nx);
    theta_img = nan(size(I));
    T_img     = nan(size(I));
    z    = z_cm(:);

    for j = 1:nx
        if isnan(z_surf(j)), continue; end
        if isempty(opts.bg)
            clear_rows = z < z_front(j) - opts.fit_margin_cm;
            if isnan(z_front(j)), clear_rows = z < z_surf(j) - opts.surf_gap_cm; end
            if nnz(clear_rows) < 10*(opts.fit_order + 1), continue; end
            zc = z(clear_rows);
            [p, ~, mu] = polyfit(zc, lnI(clear_rows, j), opts.fit_order);
            lnC = polyval(p, z, [], mu);
        else
            lnC = lnB(:, j);
        end
        th = lnC - lnI(:, j);
        if ~isempty(opts.bg), Tc = I(:, j) ./ max(opts.bg(:, j), 1); end
        if ~isempty(opts.bg) && opts.bg_normalize
            % Lighting/exposure drift between bg and frame (CoreView 95:
            % frame 1875 is ~1.5% darker than frame 1 even in deep clear
            % water) adds a constant theta at every depth, which the 99%
            % integral turns into "dye" down to the frame bottom.  Zero
            % theta on the clear water below the front.
            clear_rows = z < z_front(j) - opts.fit_margin_cm;
            if isnan(z_front(j)), clear_rows = z < z_surf(j) - opts.surf_gap_cm; end
            if any(clear_rows)
                th = th - median(th(clear_rows));
                Tc = Tc / median(Tc(clear_rows));
            end
        end
        if isempty(opts.bg), Tc = exp(-th); end
        T_img(z <= z_surf(j), j) = Tc(z <= z_surf(j));
        theta_img(z <= z_surf(j), j) = th(z <= z_surf(j));
        th(z > z_surf(j) - opts.surf_gap_cm) = 0;      % above the surface / meniscus
        TH(:, j) = interp1(z - z_surf(j), th, zeta, 'linear', NaN);
    end

    % Local z99(x): same integral on the theta averaged over each x bin
    nb    = max(1, round(10*opts.bin_cm / dx_mm));     % columns per bin
    edges = 1 : nb : nx;
    jc    = nan(1, numel(edges));
    z99_x = nan(1, numel(edges));
    for b = 1:numel(edges)
        cols  = edges(b) : min(edges(b) + nb - 1, nx);
        jc(b) = round(mean(cols));
        z99_x(b) = z99_of(zeta, mean(TH(:, cols), 2, 'omitnan'), opts.frac);
    end

    theta_bar = mean(TH, 2, 'omitnan');
    ok = ~isnan(theta_bar);
    z99  = z99_of(zeta(ok), theta_bar(ok), opts.frac);
    zeta = zeta(ok);  theta_bar = theta_bar(ok);
end

function z = z99_of(zeta, th, frac)
    % Level above which frac of the integral of th lies (zeta <= 0, top first)
    ok = ~isnan(th);
    z  = NaN;
    if nnz(ok) < 2, return; end
    zeta = zeta(ok);  th = th(ok);
    C = cumtrapz(-zeta, th);                           % integral from the surface down
    if C(end) <= 0, return; end
    k = find(C >= frac * C(end), 1, 'first');
    z = zeta(k);
end
