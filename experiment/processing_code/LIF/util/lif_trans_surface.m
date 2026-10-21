function [z_surf, S] = lif_trans_surface(I, z_cm, dx_mm, opts)
%LIF_TRANS_SURFACE Water surface in a transverse (y, z) LIF emission image.
%
%   [z_surf, S] = LIF_TRANS_SURFACE(I, z_cm, dx_mm, opts)
%
% The side-view helper lif_dye_front looks for the DARKEST row, because
% there the dye absorbs a back-lit scene.  Here the sheet excites
% rhodamine and the dye EMITS: the surface dye layer is the brightest
% feature, and plumes of dye hang below it.
%
% Finding it robustly matters because the dye detaches: in CoreView 105
% frame 1843 a plain per-column argmax lands on a plume as often as on the
% surface, and the resulting "surface" wanders by hundreds of pixels.  So:
%
%   1. smooth (opts.smooth_mm), then average along y over
%      opts.surf_smooth_cm -- a plume is ~1 cm wide, the surface line runs
%      the whole frame, so the along-y average keeps the line and washes
%      the plumes out;
%   2. per column, the brightest row of that y-averaged image inside
%      opts.search_cm (a prior band; [] = whole frame).  Columns whose
%      peak is below opts.min_counts are dropped -- the sheet does not
%      light the whole width;
%   3. a polynomial of order opts.fit_order through the surviving picks
%      vs column index, refitted opts.n_iter times after rejecting
%      residuals beyond opts.reject_sig standard deviations.  The camera
%      is rolled ~6 deg about its axis, so the surface images as a
%      straight slanted line (fit_order 1); raise fit_order if a long
%      cross-tank wave bends it.
%
% Pass opts.prior (a previous S.p) to keep the datum fixed across frames,
% e.g. from the first, undisturbed frame of a run: the fit is then skipped
% and only the offset is re-estimated when opts.track is true.
%
% Inputs:
%   I     - [ny x nx] counts, row 1 at the top (largest z)
%   z_cm  - [1 x ny] height of each row (cm)
%   dx_mm - pixel size of I (mm)
%
% opts fields (defaults):
%   smooth_mm 1.0, surf_smooth_cm 3, search_cm [] , min_counts 15,
%   fit_order 1, n_iter 3, reject_sig 2.5, prior [], track true,
%   offset_cm 0 (shifts the datum down by this much, e.g. if the emitting
%   layer sits below the free surface)
%
% Outputs:
%   z_surf - [1 x nx] surface height (cm) in every column, from the fit
%   S      - struct: p (polyfit coefficients, col -> row), mu (polyfit
%            centring), cols / rows (the accepted picks), n_used,
%            rms_px, roll_deg (tilt of the fitted line), peak (per-column
%            peak counts of the y-averaged image)

    def = struct('smooth_mm', 1.0, 'surf_smooth_cm', 3, 'search_cm', [], ...
                 'min_counts', 15, 'fit_order', 1, 'n_iter', 3, ...
                 'reject_sig', 2.5, 'prior', [], 'track', true, 'offset_cm', 0);
    if nargin < 4, opts = struct(); end
    for f = fieldnames(def).'
        if ~isfield(opts, f{1}), opts.(f{1}) = def.(f{1}); end
    end

    [ny, nx] = size(I);
    Is = imgaussfilt(I, max(opts.smooth_mm / dx_mm, 0.5));
    ny_av = max(1, round(10*opts.surf_smooth_cm / dx_mm));
    Iy = movmean(Is, ny_av, 2);                      % average along y: keeps the line

    rows = 1:ny;
    if ~isempty(opts.search_cm)
        rows = find(z_cm >= min(opts.search_cm) & z_cm <= max(opts.search_cm));
    end
    [pk, k] = max(Iy(rows, :), [], 1);
    r_pick  = rows(k);
    use     = pk(:).' >= opts.min_counts;

    cols = find(use);
    if numel(cols) < 10*(opts.fit_order + 1)
        error('lif_trans_surface:signal', ...
              'Only %d columns reach %g counts -- the sheet may be off or min_counts too high.', ...
              numel(cols), opts.min_counts);
    end

    if ~isempty(opts.prior)
        p = opts.prior.p;  mu = opts.prior.mu;
        if opts.track      % keep the shape, let the whole line ride up and down
            r_fit = polyval(p, cols, [], mu);
            d = r_pick(cols) - r_fit;
            p(end) = p(end) + median(d(abs(d - median(d)) < opts.reject_sig*1.4826*mad(d, 1) + eps));
        end
        r_all = polyval(p, 1:nx, [], mu);
        res   = r_pick(cols) - polyval(p, cols, [], mu);
    else
        for it = 1:max(1, opts.n_iter)
            [p, ~, mu] = polyfit(cols, r_pick(cols), opts.fit_order);
            res = r_pick(cols) - polyval(p, cols, [], mu);
            s   = 1.4826 * mad(res, 1);
            if it < opts.n_iter && s > 0
                keep = abs(res) <= opts.reject_sig * s;
                if nnz(keep) >= 10*(opts.fit_order + 1), cols = cols(keep); end
            end
        end
        r_all = polyval(p, 1:nx, [], mu);
    end

    % Rows -> cm: z_cm is a regular, decreasing grid, so interpolate linearly
    % and extrapolate rather than rounding to a row.
    z_surf = interp1(1:ny, z_cm(:).', r_all, 'linear', 'extrap') - opts.offset_cm;

    S = struct('p', p, 'mu', mu, 'cols', cols, 'rows', r_pick(cols), ...
               'n_used', numel(cols), 'rms_px', sqrt(mean(res.^2)), ...
               'roll_deg', atand((r_all(end) - r_all(1)) / max(nx-1, 1)), ...
               'peak', pk(:).');
end
