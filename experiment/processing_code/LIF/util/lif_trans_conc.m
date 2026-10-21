function R = lif_trans_conc(I, z_cm, z_surf, dx_mm, opts)
%LIF_TRANS_CONC Dye concentration, profile and z99 from a transverse LIF image.
%
%   R = LIF_TRANS_CONC(I, z_cm, z_surf, dx_mm, opts)
%
% This is the emission (true LIF) counterpart of lif_z99.  The side views
% are back-lit and the dye ABSORBS, so there the concentration has to be
% recovered through Beer-Lambert, theta = ln(I_bg/I), with a clear-water
% reference per column.  A laser sheet through rhodamine FLUORESCES, so in
% the dilute, optically thin limit
%
%       C(y, z)  propto  I(y, z) - I_dark
%
% directly -- no background frame, no logarithm.  Everything below is
% therefore about the two things that do bite:
%
%   Noise floor.  The frames are black (median 0 counts, clear water
%     ~0.3 counts of read noise) but the 99% integral runs over thousands
%     of rows, so an un-removed 0.3 counts/row swamps the dye tail and
%     drags z99 to the bottom of the frame.  opts.deep_cm sets a band of
%     clear water well below the dye whose median is subtracted from each
%     column.  Negatives are KEPT so the residual noise averages out
%     instead of biasing the tail upward (same reasoning as lif_z99).
%
%   Integration depth.  The paper takes H as the bottom of the camera
%     frame, because there the frame was matched to the mixed layer.  Here
%     the frame is ~19 cm tall and the dye reaches ~2 cm, so 90% of the
%     integration path is clear water: the floor subtraction above is
%     never perfect (the sheet scatters a faint glow), and the 1% the
%     99% level is looking for is the same size as what is left.  So the
%     integral stops opts.int_depth_cm below the surface -- far below the
%     dye, well above the noise-dominated remainder.  It is a real choice:
%     print z99 for a couple of values before trusting a number.
%
%   Sheet non-uniformity.  The laser enters from one side and is absorbed
%     on the way across, and the sheet waist varies, so counts are not
%     comparable between columns.  There is no uniform-dye calibration
%     frame, but the optics do not move during a run, so the sheet profile
%     S(y) is steady and the reference frame measures it: frame 1 is an
%     undisturbed dye layer of uniform concentration lying on the surface,
%     so its near-surface counts ARE S(y).  opts.normalize = 'ref' divides
%     every column by that, which both removes the sheet profile and
%     leaves C in units of the initial surface concentration.
%
%     The alternatives are worse here and are kept only for comparison.
%     'column' (divide each column by its own depth-integral) forces every
%     y to weigh the same -- but a column the wind has stripped of dye then
%     has its leftover noise rescaled to unit integral, and in frame 1843
%     that alone moved z99 from -1.9 cm to -4.5 cm.  'peak' divides by the
%     column's own maximum, with the same flaw in milder form; 'none'
%     leaves raw counts and lets the brightest part of the sheet dominate
%     the y-average.
%
% z99 (paper eq. 3.14) is then the usual
%       int_{z99}^{0} C dz = frac * int_{-H}^{0} C dz
% with z = 0 at the LOCAL surface, so the camera roll and any cross-tank
% wave are removed before the average.  This plane IS the paper's
% transverse one, so the frame-wide z99 (y averaged) is the directly
% comparable number; z99 vs y is returned too, in opts.bin_cm bins.
%
% Inputs:
%   I      - [ny x nx] counts (row 1 = top)
%   z_cm   - [1 x ny] row heights (cm)
%   z_surf - [1 x nx] surface height per column (cm), from lif_trans_surface
%   dx_mm  - pixel size of I (mm)
%
% opts fields (defaults):
%   frac 0.99, smooth_mm 1.0, surf_gap_cm 0 (skip this much below the
%   surface, e.g. to drop a saturated surface line), top_cm 0.3 (band used
%   for the 'peak' normalisation), deep_cm 6 (clear water starts this far
%   below the surface -> noise floor; it has to be below the deepest dye,
%   so run this on the WHOLE frame and crop only for display),
%   deep_min_rows 200 (a column with less clear water than that is skipped:
%   its floor would be guesswork),
%   ref [] (the reference frame, same size as I -- the undisturbed dye sheet
%   lying on the surface.  Which columns the laser lights is a property of
%   the optics, not of the frame, so the lit mask is taken from the
%   reference: a column the wind has stripped of surface dye still belongs
%   in the average, while a column the sheet never reached does not, and
%   only the reference can tell them apart),
%   min_peak 8 (counts: the threshold that mask uses; without ref it is
%   applied to the frame's own brightest pixel instead),
%   int_depth_cm 5 (the 99% integral stops here below the surface;
%   [] = the whole frame, as in the paper),
%   normalize 'ref' (needs opts.ref; 'column' | 'peak' | 'none'), bin_cm 1,
%   min_bin_frac 0.5 (a z99(y) bin needs this fraction of its columns lit,
%   else NaN), dark [] (scalar counts or an image),
%   pen_frac 0.02 (dye-penetration contour, fraction of the column peak),
%   clip_above true (blank everything above the surface: the reflection of
%   the dye layer in the free surface is not dye)
%
% Output struct R:
%   C_img    [ny x nx] concentration in image coordinates, in the units set
%            by opts.normalize (NaN above the surface)
%   Craw_img [ny x nx] dark- and floor-corrected counts, same blanking
%   zeta     [nz x 1] depth grid (cm, 0 = local surface, negative down)
%   C_bar    [nz x 1] y-averaged concentration on zeta
%   z99      frame-wide z99 (cm, negative = below the surface)
%   z99_y    [1 x nbin] local z99 in each y bin
%   jc       [1 x nbin] centre column of each bin
%   z_pen    [1 x nx] deepest level reaching pen_frac of the column peak (cm)
%   peak_col [1 x nx] brightest counts anywhere in the column (NaN where the
%            sheet does not light it)
%   surf_col [1 x nx] counts in this frame's top_cm band: how much surface
%            dye is left
%   ref_col  [1 x nx] the same band in the reference frame: the sheet
%            profile S(y), and the divisor behind normalize = 'ref'
%   floor_col[1 x nx] subtracted noise floor per column (counts)

    def = struct('frac', 0.99, 'smooth_mm', 1.0, 'surf_gap_cm', 0, ...
                 'top_cm', 0.3, 'deep_cm', 6, 'deep_min_rows', 200, ...
                 'min_peak', 8, 'ref', [], 'int_depth_cm', 5, 'normalize', 'ref', ...
                 'bin_cm', 1, 'min_bin_frac', 0.5, 'dark', [], ...
                 'pen_frac', 0.02, 'pen_min_counts', 2, 'pen_smooth_cm', 0.5, ...
                 'clip_above', true);
    if nargin < 5, opts = struct(); end
    for f = fieldnames(def).'
        if ~isfield(opts, f{1}), opts.(f{1}) = def.(f{1}); end
    end

    if strcmpi(opts.normalize, 'ref') && isempty(opts.ref)
        error('lif_trans_conc:ref', 'normalize = ''ref'' needs opts.ref (the reference frame).');
    end
    if ~isempty(opts.dark), I = I - opts.dark; end
    [ny, nx] = size(I);
    sig = max(opts.smooth_mm / dx_mm, 0.5);
    Is  = imgaussfilt(I, sig);
    z   = z_cm(:);

    Rs = [];
    if ~isempty(opts.ref)
        R0 = opts.ref;
        if ~isempty(opts.dark), R0 = R0 - opts.dark; end
        Rs = imgaussfilt(R0, sig);
    end

    dz    = abs(z_cm(2) - z_cm(1));
    z_bot = min(z_cm) - max(z_surf, [], 'omitnan');
    if ~isempty(opts.int_depth_cm), z_bot = max(z_bot, -abs(opts.int_depth_cm)); end
    zeta = (0 : -dz : z_bot).';
    CC   = nan(numel(zeta), nx);

    C_img     = nan(ny, nx);
    Craw_img  = nan(ny, nx);
    z_pen     = nan(1, nx);
    peak_col  = nan(1, nx);
    surf_col  = nan(1, nx);
    ref_col   = nan(1, nx);
    floor_col = nan(1, nx);

    for j = 1:nx
        if isnan(z_surf(j)), continue; end
        c = Is(:, j);

        deep = z < z_surf(j) - opts.deep_cm;        % clear water -> noise floor
        if nnz(deep) < opts.deep_min_rows
            continue                                % not enough clear water to trust
        end
        floor_col(j) = median(c(deep));
        c = c - floor_col(j);

        wet = z <= z_surf(j) - opts.surf_gap_cm;
        if ~any(wet), continue; end

        % Normalise by the brightest point of the WHOLE column, not of the
        % near-surface band: once the wind strips the surface layer the band
        % holds almost nothing, and dividing the plume below by it sends
        % C/C_peak far past 1 and whites the panel out.  The near-surface
        % value is still returned separately, as the sheet's y profile.
        top = wet & z >= z_surf(j) - opts.surf_gap_cm - opts.top_cm;
        surf_col(j) = max(c(top));
        peak_col(j) = max(c(wet));

        Craw_img(wet, j) = c(wet);

        if isempty(Rs)
            ref_col(j) = peak_col(j);
        else
            % Same band, same floor treatment, on the reference frame: in
            % frame 1 all the dye is still at the surface, so this measures
            % the sheet's y profile and nothing else.
            r = Rs(:, j) - median(Rs(deep, j));
            ref_col(j) = max(r(top));
        end
        if ~(ref_col(j) >= opts.min_peak)
            peak_col(j) = NaN;                      % column not lit by the sheet
            ref_col(j)  = NaN;
            continue
        end

        win = wet & z >= z_surf(j) + z_bot;         % the z99 integration window
        cw  = c;  cw(~win) = 0;
        switch lower(opts.normalize)
            case 'ref'
                cn = c / ref_col(j);                % units of the frame-1 surface dye
            case 'column'
                s = trapz(flipud(z), flipud(cw));   % z decreases down the rows
                if ~(s > 0), continue; end
                cn = c / s;                         % 1/cm: integrates to 1 over win
            case 'peak'
                cn = c / peak_col(j);
            case 'none'
                cn = c;
            otherwise
                error('lif_trans_conc:normalize', 'Unknown normalize ''%s''.', opts.normalize);
        end

        C_img(wet, j) = cn(wet);

        % Dye penetration: deepest row still at pen_frac of the column peak.
        % The absolute floor matters as much as the fraction -- without it a
        % single hot pixel 10 cm down passes in a dim column and the contour
        % drops to the bottom of the frame.
        thr = max(opts.pen_frac * peak_col(j), opts.pen_min_counts);
        deepest = find(wet & c >= thr, 1, 'last');
        if ~isempty(deepest), z_pen(j) = z(deepest); end

        cz = cn;  cz(~wet) = 0;
        CC(:, j) = interp1(z - z_surf(j), cz, zeta, 'linear', NaN);
    end

    if ~opts.clip_above
        C_img(isnan(C_img) & ~isnan(Craw_img)) = NaN;   % no-op, kept for clarity
    end

    if opts.pen_smooth_cm > 0
        z_pen = movmedian(z_pen, max(1, round(10*opts.pen_smooth_cm / dx_mm)), ...
                          'omitnan', 'Endpoints', 'shrink');
    end

    % z99 vs y
    nb    = max(1, round(10*opts.bin_cm / dx_mm));
    edges = 1 : nb : nx;
    jc    = nan(1, numel(edges));
    z99_y = nan(1, numel(edges));
    for b = 1:numel(edges)
        cols  = edges(b) : min(edges(b) + nb - 1, nx);
        jc(b) = round(mean(cols));
        % A bin made of two or three lit columns gives a z99 that is mostly
        % noise, which then dominates the 5-95% band in the profile panel.
        if nnz(~isnan(peak_col(cols))) < opts.min_bin_frac * numel(cols), continue; end
        z99_y(b) = z99_of(zeta, mean(CC(:, cols), 2, 'omitnan'), opts.frac);
    end

    C_bar = mean(CC, 2, 'omitnan');
    ok    = ~isnan(C_bar);
    z99   = z99_of(zeta(ok), C_bar(ok), opts.frac);

    R = struct('C_img', C_img, 'Craw_img', Craw_img, 'zeta', zeta(ok), ...
               'C_bar', C_bar(ok), 'z99', z99, 'z99_y', z99_y, 'jc', jc, ...
               'z_pen', z_pen, 'peak_col', peak_col, 'surf_col', surf_col, 'ref_col', ref_col, ...
               'floor_col', floor_col, ...
               'normalize', opts.normalize, 'int_depth_cm', -z_bot);
end

function z = z99_of(zeta, c, frac)
    % Level above which frac of the integral of c lies (zeta <= 0, top first)
    ok = ~isnan(c);
    z  = NaN;
    if nnz(ok) < 2, return; end
    zeta = zeta(ok);  c = c(ok);
    C = cumtrapz(-zeta, c);                % integral from the surface down
    if C(end) <= 0, return; end
    k = find(C >= frac * C(end), 1, 'first');
    z = zeta(k);
end
