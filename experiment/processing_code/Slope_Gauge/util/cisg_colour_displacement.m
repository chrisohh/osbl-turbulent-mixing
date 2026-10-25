function [dx, dy, info, D] = cisg_colour_displacement(ref_image, obs_image, setup)
%CISG_COLOUR_DISPLACEMENT Local colour -> position inversion of the CISG target.
%
%   [dx, dy, info] = cisg_colour_displacement(ref_image, obs_image, setup)
%   [dx, dy, info, D] = cisg_colour_displacement(...)        % D: per-pixel diagnostics (ROI-sized)
%
% For every pixel x of the observed image, find the shift delta (pixels) such
% that the reference colour at x + delta equals the observed colour at x,
%
%   C_ref(x + delta) = C_obs(x),   C = [R G] (smoothed, optionally R/(R+G+B), G/(R+G+B)), or [R G B].
%
% Newton / Lucas-Kanade iterations with the local Jacobian of the reference,
%   J = [dR/dx dR/dy; dG/dx dG/dy],   delta <- delta + J^-1 (C_obs(x) - C_ref(x + delta)),
% or with [R G B] (colour_channels 'RGB') the 3 x 2 Jacobian and the least-squares (Gauss-Newton) step
%   delta <- delta + (J'J)^-1 J' (C_obs(x) - C_ref(x + delta)).
% sampled at x + delta (linear interpolation), starting from delta = 0.  This
% replaces the global linear colour gradient of the old path: it follows the
% real (non-linear) colour map and handles R/G cross-talk.
%
% dx = shift along image columns (+ right), dy = shift along image rows (+ down), in
% pixels (the sign used by the linear path: delta = dC / gradient).  NaN where
% the inversion is not trustworthy: J ill-conditioned (|det J| small or
% cond(J) large, e.g. the non-monotonic red peak), not converged, or a large
% final colour residual.
%
% setup fields (defaults):
%   colour_channels  'RG'    'RG' | 'RGB'.  'RGB' adds blue where the R/G map folds (a channel peaks: |det J| ~ 0), e.g. the
%                            L along the left / top of the run 101 target, where blue lifts sigma_min(J) ~4x (cisg_pipeline 5f).
%                            Dye-free runs only: rhodamine changes blue most (signature d, cisg_dye_unmix); not with 'sum'
%                            or dye_model 'joint'.  The |det J| / cond tests use sigma1 * sigma2 and sigma1 / sigma2 of J.
%                            'hybrid': solve with R, G; re-solve with R, G, B (from 0) only where the R, G solve failed
%                            (any test, incl. sigma_S) or sigma_min(J_RG) at its solution is below hybrid_frac x the median
%                            over the target (the folds); the R, G, B result is taken where it passes or R, G failed.
%                            On steep wave faces blue disagrees with R, G (a radial, lens-like pattern: cisg_pipeline 5i),
%                            so R, G is kept wherever it is well conditioned.  info.hybrid_fraction = share taken from RGB.
%   hybrid_frac      0.3     'hybrid': sigma_min(J_RG) threshold, fraction of its median over the target
%   colour_norm      'none'  'none' | 'sum' (R/(R+G+B), G/(R+G+B): cancels brightness changes)
%   colour_sigma     12      Gaussian smoothing of both images (px)
%   max_disp_px      600     clamp on |delta|
%   colour_n_iter    15
%   det_min_frac     0.1     mask |det J| < frac * median |det J| of the reference
%   cond_max         10      mask cond(J) above this
%   res_max          0.05    mask final colour residual (|C_obs - C_ref|, 0..1 units)
%   conv_px          2        mask last Newton step above this (px)
%   drift_corr       'none'  'gain': per-channel robust fit obs = a ref + b over the ROI, obs <- (obs - b)/a (lighting drift);
%                            'offset': b only (a = 1), b = median(obs - ref) per channel (frame-to-frame brightness flicker)
%   dye_mask         false   mask cross-stream dye stripes: rows whose median blue / brightness residual (obs vs ref at x+delta,
%                            0..1), minus its running median over dye_win rows (101), exceeds dye_thresh = [tB tS] (default
%                            [0.010 0.004]) are masked, +- dye_dilate_px rows (6).  R6G dims B and brightness, which [R G] ignores
%   dye_model        'none'  'joint': solve (shift, dye amount c) per pixel from R,G,B (cisg_dye_unmix.m); c in info.dye_c.
%                            Noisier than 'none' (adds ~8% rms to Sx on clean frames), see calibration/cisg_dye_unmix_check.m
%   init_dx, init_dy []      ROI-sized starting shift (px, NaN = none), e.g. the previous frame's: those pixels are solved
%                            from it AND from 0, and per pixel the passing / lower-residual solution is kept (warm start;
%                            info.warm_fraction = share of kept pixels that came from it)
%   roi              []      struct row_min/row_max/col_min/col_max: only these pixels
%                            are solved (output is still full size, NaN elsewhere)
%   target_cols      []      [left right] columns of the target (median |det J| uses them)
%
% info: mask_fraction (all reasons), dye_fraction (pixels masked by the dye test), thresholds, runtime.
% D (4th output): ROI-sized bits (1 x+delta out of bounds, 2 displacement clamped, 3 ill-conditioned J,
%   4 not converged, 5 colour residual, 6 dye), dx, dy (unmasked), eB, eS (blue / brightness residuals, 0..1),
%   e (rows x cols x channels: residual C_obs - C_ref(x + delta) of each solved channel, 0..1).
% The smoothed reference and its gradients are cached between calls.

    t0 = tic;
    o = getopts(setup);
    [H, W, ~] = size(ref_image);
    if isempty(o.roi)
        r1 = 1; r2 = H; c1 = 1; c2 = W;
    else
        r1 = o.roi.row_min; r2 = o.roi.row_max; c1 = o.roi.col_min; c2 = o.roi.col_max;
    end

    %% Reference: smoothed colours and Jacobian (cached)
    persistent cache
    need_bs = o.dye_mask || nargout > 3;                    % blue and brightness checks
    nch = o.nch;
    key = [H W o.sigma double(strcmp(o.norm, 'sum')) sum(double(ref_image(1:997:end))) nch];
    if isempty(cache) || ~isequal(cache.key, key) || (need_bs && ~isfield(cache, 'CB'))
        C = prep(ref_image, o);
        cache = struct('key', key, 'C', {cell(1, nch)}, 'Gx', {cell(1, nch)}, 'Gy', {cell(1, nch)});
        for k = 1:nch
            cache.C{k} = C(:,:,k);  [cache.Gx{k}, cache.Gy{k}] = gradient(C(:,:,k));
        end
        if need_bs
            [cache.CB, cache.CS] = blue_bright(ref_image, o);
        end
    end
    cs = 1:W;  if ~isempty(o.target_cols), cs = max(1, round(o.target_cols(1))):min(W, round(o.target_cols(2))); end
    sub = @(A) A(1:16:end, cs(1:16:end));
    % |det J| of R, G, and sigma1 * sigma2 = sqrt(det J'J) of R, G, B: the thresholds are relative to their medians
    dd2 = abs(sub(cache.Gx{1}) .* sub(cache.Gy{2}) - sub(cache.Gy{1}) .* sub(cache.Gx{2}));
    dd3 = [];
    if nch == 3
        p = 0;  q = 0;  r = 0;
        for k = 1:3, gx = sub(cache.Gx{k});  gy = sub(cache.Gy{k});  p = p + gx.^2;  q = q + gx .* gy;  r = r + gy.^2; end
        dd3 = sqrt(max(p .* r - q.^2, 0));
    end
    smin_med = [];
    if o.hybrid                                             % typical sigma_min of J_RG over the target
        gx1 = sub(cache.Gx{1});  gy1 = sub(cache.Gy{1});  gx2 = sub(cache.Gx{2});  gy2 = sub(cache.Gy{2});
        smin_med = median(smin_rg(gx1, gy1, gx2, gy2), 'all', 'omitnan');
    end
    det_frac = o.det_min_frac;
    % physical precision mask: propagate the camera noise through J^-1 to the slope uncertainty (needs setup.noise,
    % from cisg_noise_estimate, and the slope scale fields); the relative |det J| / cond tests are then only loose guards
    use_sig = ~isempty(o.noise) && strcmp(o.norm, 'none') && all(isfield(setup, {'pattern_scale', 'cm_per_pixel_x', 'cm_per_pixel_y', 'L_eff_cm', 'n_water', 'n_air'}));
    cond_eff = o.cond_max;  kS = [];
    if use_sig
        det_frac = 0.003;  cond_eff = 100;
        kS = [setup.pattern_scale * setup.cm_per_pixel_x, setup.pattern_scale * setup.cm_per_pixel_y] / ((setup.n_water - setup.n_air) * setup.L_eff_cm);
    end
    det_rg = det_frac * median(dd2(:));
    det_all = det_rg;  if nch == 3, det_all = det_frac * median(dd3(:)); end
    det_min = det_rg;  if strcmp(o.channels, 'RGB'), det_min = det_all; end   % reported in info

    %% Observed: smoothed colours (ROI + margin)
    pad = ceil(4 * o.sigma);
    pr = max(1, r1 - pad):min(H, r2 + pad);  pc = max(1, c1 - pad):min(W, c2 + pad);
    Co = prep(obs_image(pr, pc, :), o);
    Co = Co(r1 - pr(1) + 1 : r2 - pr(1) + 1, c1 - pc(1) + 1 : c2 - pc(1) + 1, :);
    if need_bs
        [OB, OS] = blue_bright(obs_image(pr, pc, :), o);
        OB = OB(r1 - pr(1) + 1 : r2 - pr(1) + 1, c1 - pc(1) + 1 : c2 - pc(1) + 1);
        OS = OS(r1 - pr(1) + 1 : r2 - pr(1) + 1, c1 - pc(1) + 1 : c2 - pc(1) + 1);
        bits = zeros(r2 - r1 + 1, c2 - c1 + 1, 'uint8');
        eBm = nan(r2 - r1 + 1, c2 - c1 + 1, 'single');  eSm = eBm;
        dxm = eBm;  dym = eBm;
        em = nan(r2 - r1 + 1, c2 - c1 + 1, nch, 'single');  % residual per solved channel at the final shift
    end

    %% Lighting-drift correction: obs ~ a * ref + b per channel (trimmed fit, waves do not bias it)
    % 'offset': obs ~ ref + b, b = median(obs - ref) per channel (a uniform brightness step, e.g. the ~0.15-count
    % frame-to-frame flicker of run 101, reads as a shift wherever the colour map is flat or has a sharp feature)
    drift_ab = repmat([1 0], nch, 1);
    if any(strcmp(o.drift, {'gain', 'offset'}))
        for k = 1:nch
            Rk = cache.C{k}(r1:r2, c1:c2);
            rcs = true(1, c2 - c1 + 1);
            if ~isempty(o.target_cols), rcs = (c1:c2) >= o.target_cols(1) + 100 & (c1:c2) <= o.target_cols(2) - 100; end
            xv = Rk(1:6:end, rcs(1:1:end)); xv = xv(:, 1:6:end); yv = Co(1:6:end, :, k); yv = yv(:, rcs); yv = yv(:, 1:6:end);
            xv = double(xv(:)); yv = double(yv(:));
            if strcmp(o.drift, 'offset')
                p = [1, median(yv - xv, 'omitnan')];
            else
                keep = true(size(xv));
                for it = 1:6
                    p = polyfit(xv(keep), yv(keep), 1);
                    res = abs(yv - polyval(p, xv));
                    keep = res <= quantile(res, 0.7);
                end
            end
            drift_ab(k, :) = p;                                % [a b]
            Co(:, :, k) = (Co(:, :, k) - p(2)) / p(1);
        end
    end

    %% Newton iterations (row chunks)
    dx = nan(H, W, 'single');  dy = nan(H, W, 'single');
    sig_sx = [];  sig_sy = [];
    if use_sig, sig_sx = nan(r2 - r1 + 1, c2 - c1 + 1, 'single');  sig_sy = sig_sx; end
    chunk = 256;
    n_tot = 0;  n_bad = 0;  n_warm = 0;  n_hyb = 0;
    if ~isempty(o.init_dx) && ~isequal(size(o.init_dx), [r2 - r1 + 1, c2 - c1 + 1])
        warning('cisg_colour_displacement:init', 'init_dx is %dx%d, ROI is %dx%d: warm start ignored', ...
                size(o.init_dx, 1), size(o.init_dx, 2), r2 - r1 + 1, c2 - c1 + 1);
        o.init_dx = [];  o.init_dy = [];
    end
    for rs = r1:chunk:r2
        re = min(rs + chunk - 1, r2);
        [X, Y] = meshgrid(single(c1:c2), single(rs:re));
        Mall = struct('X', X, 'Y', Y, 'O', Co(rs - r1 + 1 : re - r1 + 1, :, :), 'C', {cache.C}, 'Gx', {cache.Gx}, 'Gy', {cache.Gy}, ...
                      'W', W, 'H', H, 'det_min', det_all, 'cond', cond_eff);
        M = Mall;                                            % the main solve: all channels for 'RGB', R, G otherwise
        if ~strcmp(o.channels, 'RGB')
            M.O = M.O(:, :, 1:2);  M.C = M.C(1:2);  M.Gx = M.Gx(1:2);  M.Gy = M.Gy(1:2);  M.det_min = det_rg;
        end
        z = zeros(size(X), 'single');
        s = newton(z, z, M, o);                              % cold start: delta = 0
        if ~isempty(o.init_dx)                               % warm start from the previous frame's shift, where it had one
            ii = rs - r1 + 1 : re - r1 + 1;
            ix = single(o.init_dx(ii, :));  iy = single(o.init_dy(ii, :));
            has = isfinite(ix) & isfinite(iy);
            if any(has(:))
                ix(~has) = 0;  iy(~has) = 0;
                w = newton(ix, iy, M, o);
                % per pixel: a passing solution beats a failing one, then the smaller colour residual wins;
                % keeping the cold solve means a wrong previous shift cannot carry over by itself
                use = has & ((w.pass & ~s.pass) | (w.pass == s.pass & w.res < s.res));
                for f = fieldnames(s).'                      % gx, gy, e are rows x cols x channels
                    v = s.(f{1});  wv = w.(f{1});  m = repmat(use, 1, 1, size(v, 3));  v(m) = wv(m);  s.(f{1}) = v;
                end
                n_warm = n_warm + nnz(use & s.pass);
            end
        end
        s = with_precision(s, M.C, o, kS, use_sig);
        if o.hybrid
            % re-solve with R, G, B where R, G failed or sits on a fold of its colour map
            s.e(:, :, 3) = NaN;                              % blue residual only where the RGB solve is taken
            weak = ~s.good | smin_rg(s.gx(:,:,1), s.gy(:,:,1), s.gx(:,:,2), s.gy(:,:,2)) < o.hybrid_frac * smin_med;
            idx = find(weak);
            if ~isempty(idx)
                Ms = Mall;  Ms.X = Mall.X(idx);  Ms.Y = Mall.Y(idx);  Ms.O = zeros(numel(idx), 1, 3, 'single');
                for k = 1:3, Ok = Mall.O(:, :, k);  Ms.O(:, 1, k) = Ok(idx); end
                zs = zeros(numel(idx), 1, 'single');
                w = with_precision(newton(zs, zs, Ms, o), Ms.C, o, kS, use_sig);
                take = w.good | ~s.good(idx);                % RGB where it passes, or where R, G failed anyway
                t = idx(take);
                for f = {'ux', 'uy', 'step', 'xs', 'ys', 'ok', 'oob', 'res', 'pass', 'sSx', 'sSy', 'good'}
                    v = w.(f{1});  s.(f{1})(t) = v(take);
                end
                for k = 1:3
                    ek = s.e(:, :, k);  wk = w.e(:, 1, k);  ek(t) = wk(take);  s.e(:, :, k) = ek;
                end
                n_hyb = n_hyb + nnz(take & w.good);
            end
        end
        ux = s.ux;  uy = s.uy;  step = s.step;  xs = s.xs;  ys = s.ys;
        ok = s.ok;  oob = s.oob;
        good = s.good;
        precok = true(size(good));
        if use_sig
            precok = max(s.sSx, s.sSy) <= o.max_slope_sigma;
            ii = rs - r1 + 1 : re - r1 + 1;
            sig_sx(ii, :) = single(s.sSx);  sig_sy(ii, :) = single(s.sSy);
        end
        if need_bs
            ii = rs - r1 + 1 : re - r1 + 1;
            eBm(ii, :) = OB(ii, :) - interp2(cache.CB, xs, ys, 'linear');
            eSm(ii, :) = OS(ii, :) - interp2(cache.CS, xs, ys, 'linear');
            dxm(ii, :) = ux;  dym(ii, :) = uy;  em(ii, :, :) = s.e;
            bits(ii, :) = uint8(oob) + 2 * uint8(abs(ux) >= o.max_disp_px | abs(uy) >= o.max_disp_px) + 4 * uint8(~ok) ...
                          + 8 * uint8(step > o.conv_px) + 16 * uint8(s.res > o.res_max) + 128 * uint8(~precok);
        end
        ux(~good) = NaN;  uy(~good) = NaN;
        dx(rs:re, c1:c2) = ux;  dy(rs:re, c1:c2) = uy;
        n_tot = n_tot + numel(good);  n_bad = n_bad + nnz(~good);
    end

    dye_fraction = 0;  dye_rows = false(r2 - r1 + 1, 1);
    if o.dye_mask
        % Dye (R6G) crosses the view as cross-stream stripes: B and brightness drop a little (~2%) along whole
        % rows, too little to flag single pixels.  Row-wise test: the row median of the blue / brightness residual
        % (valid pixels), minus its running median over dye_win rows, exceeds dye_thresh -> the whole row (+- dilation) is masked.
        good_px = bits == 0;                                 % pixels the inversion accepted (before the dye test)
        fr_ok = mean(good_px, 2) > 0.2;
        eBr = eBm;  eBr(~good_px) = NaN;  eSr = eSm;  eSr(~good_px) = NaN;
        rB = median(eBr, 2, 'omitnan');  rS = median(eSr, 2, 'omitnan');
        rB(~fr_ok) = NaN;  rS(~fr_ok) = NaN;
        devB = rB - movmedian(rB, o.dye_win, 'omitnan');
        devS = rS - movmedian(rS, o.dye_win, 'omitnan');
        dye_rows = abs(devB) > o.dye_thresh(1) | abs(devS) > o.dye_thresh(2);
        dye_rows = imdilate(dye_rows, ones(2 * round(o.dye_dilate) + 1, 1));
        dx_roi = dx(r1:r2, c1:c2);  dy_roi = dy(r1:r2, c1:c2);
        sel = repmat(dye_rows, 1, c2 - c1 + 1) & isfinite(dx_roi);
        dx_roi(dye_rows, :) = NaN;  dy_roi(dye_rows, :) = NaN;
        dx(r1:r2, c1:c2) = dx_roi;  dy(r1:r2, c1:c2) = dy_roi;
        dye_fraction = nnz(sel) / numel(sel);
        n_bad = n_bad + nnz(sel);
    end

    if ~isempty(o.dark_mask)                                 % dark / outside the target in the reference itself
        dm = o.dark_mask(r1:r2, c1:c2);
        dxr = dx(r1:r2, c1:c2);  dyr = dy(r1:r2, c1:c2);
        n_bad = n_bad + nnz(dm & isfinite(dxr));
        dxr(dm) = NaN;  dyr(dm) = NaN;
        dx(r1:r2, c1:c2) = dxr;  dy(r1:r2, c1:c2) = dyr;
        if need_bs, bits(dm) = bits(dm) + 64; end
    end
    dye_c = [];  dye_applied = [];
    if strcmp(o.dye_model, 'joint')                          % per-pixel (shift, dye amount) solve from all 3 channels
        s_u = setup;  s_u.roi = struct('row_min', r1, 'row_max', r2, 'col_min', c1, 'col_max', c2);
        [dx, dy, dye_c, uinfo] = cisg_dye_unmix(ref_image, obs_image, s_u, dx, dy);   % hybrid: keeps the 2-channel mask
        n_bad = nnz(isnan(dx(r1:r2, c1:c2)));
        dye_applied = uinfo.applied;
    end

    info = struct('mask_fraction', n_bad / n_tot, 'dye_fraction', dye_fraction, 'det_min', det_min, 'opts', o, ...
                  'drift_ab', drift_ab, 'dye_c', dye_c, 'dye_applied', dye_applied, 'sig_sx', sig_sx, 'sig_sy', sig_sy, ...
                  'sigma_mask', use_sig, 'warm_fraction', n_warm / max(n_tot, 1), 'channels', o.channels, ...
                  'hybrid_fraction', n_hyb / max(n_tot, 1), 'seconds', toc(t0));
    if nargout > 3
        D = struct('bits', bits, 'dx', dxm, 'dy', dym, 'eB', eBm, 'eS', eSm, 'e', em);
        D.dye_rows = dye_rows;  if o.dye_mask, D.bits(dye_rows, :) = D.bits(dye_rows, :) + 32; end
    end
end

function [CB, CS] = blue_bright(img, o)
    img = single(img);
    CB = imgaussfilt(img(:,:,3) / 255, o.sigma);
    CS = imgaussfilt(mean(img, 3) / 255, o.sigma);
end

function s = newton(ux, uy, M, o)
    % Newton iterations from (ux, uy), then the final residual / conditioning at x + delta
    step = zeros(size(ux), 'single');
    for it = 1:o.n_iter
        xs = min(max(M.X + ux, 1), M.W);  ys = min(max(M.Y + uy, 1), M.H);
        [e, gx, gy] = sample(M, xs, ys);
        [ok, sx, sy] = gn_step(e, gx, gy, M);
        ux = min(max(ux + sx, -o.max_disp_px), o.max_disp_px);
        uy = min(max(uy + sy, -o.max_disp_px), o.max_disp_px);
        step = hypot(sx, sy);
    end
    xs = min(max(M.X + ux, 1), M.W);  ys = min(max(M.Y + uy, 1), M.H);
    [e, gx, gy] = sample(M, xs, ys);
    ok = gn_step(e, gx, gy, M);
    oob = ~((M.X + ux >= 1) & (M.X + ux <= M.W) & (M.Y + uy >= 1) & (M.Y + uy <= M.H));
    if size(e, 3) == 2, res = hypot(e(:,:,1), e(:,:,2)); else, res = sqrt(sum(e.^2, 3)); end
    pass = ok & res <= o.res_max & step <= o.conv_px & isfinite(ux) & isfinite(uy) & ~oob;
    s = struct('ux', ux, 'uy', uy, 'step', step, 'xs', xs, 'ys', ys, 'e', e, ...
               'gx', gx, 'gy', gy, 'ok', ok, 'oob', oob, 'res', res, 'pass', pass);
end

function [e, gx, gy] = sample(M, xs, ys)
    % colour residual C_obs(x) - C_ref(x + delta) and the reference gradients at x + delta, per channel
    n = numel(M.C);
    e = zeros([size(xs) n], 'single');  gx = e;  gy = e;
    for k = 1:n
        e(:,:,k)  = M.O(:,:,k) - interp2(M.C{k}, xs, ys, 'linear');
        gx(:,:,k) = interp2(M.Gx{k}, xs, ys, 'linear');
        gy(:,:,k) = interp2(M.Gy{k}, xs, ys, 'linear');
    end
end

function [ok, sx, sy] = gn_step(e, gx, gy, M)
    % step delta: J^-1 e for [R G] (2 x 2, as originally written), (J'J)^-1 J' e for [R G B] (3 x 2)
    if size(e, 3) == 2
        a = gx(:,:,1);  b = gy(:,:,1);  c = gx(:,:,2);  d = gy(:,:,2);
        [ok, dtm] = conditioning(a, b, c, d, M.det_min, M.cond);
        sx = (d .* e(:,:,1) - b .* e(:,:,2)) ./ dtm;
        sy = (a .* e(:,:,2) - c .* e(:,:,1)) ./ dtm;
    else
        p = sum(gx.^2, 3);  q = sum(gx .* gy, 3);  r = sum(gy.^2, 3);    % J'J = [p q; q r]
        bx = sum(gx .* e, 3);  by = sum(gy .* e, 3);                     % J' e
        dA = p .* r - q.^2;                                               % det J'J = (sigma1 sigma2)^2
        s12 = sqrt(max(dA, 0));
        s1sq = (p + r + sqrt((p - r).^2 + 4 * q.^2)) / 2;                 % sigma1^2
        ok = s12 >= M.det_min & s1sq ./ max(s12, eps('single')) <= M.cond & isfinite(dA);
        dA(~ok) = 1;
        sx = (r .* bx - q .* by) ./ dA;
        sy = (p .* by - q .* bx) ./ dA;
    end
    sx(~ok) = 0;  sy(~ok) = 0;
end

function s = with_precision(s, C, o, kS, use_sig)
    % slope uncertainty sSx, sSy of a solve and s.good = passed every test incl. sigma_S (sSx = NaN without the noise model).
    % Cov(delta) = (J' W J)^-1, W = diag(1 / sigma_c^2), sigma_c = a_c + b_c C_c(x + delta) per channel;
    % with R, G only this is J^-1 diag(sigma_c^2) J^-T (the original formula, kept as written)
    if ~use_sig
        s.sSx = nan(size(s.ux), 'single');  s.sSy = s.sSx;  s.good = s.pass;
        return;
    end
    n = numel(C);  sc = cell(1, n);
    for k = 1:n, sc{k} = max(o.noise.a(k) + o.noise.b(k) * interp2(C{k}, s.xs, s.ys, 'linear'), 1e-5); end
    if n == 2
        ja = s.gx(:,:,1);  jb = s.gy(:,:,1);  jc = s.gx(:,:,2);  jd = s.gy(:,:,2);
        dt2 = max(abs(ja .* jd - jb .* jc), eps('single'));
        s.sSx = kS(1) * sqrt(jd.^2 .* sc{1}.^2 + jb.^2 .* sc{2}.^2) ./ dt2;
        s.sSy = kS(2) * sqrt(jc.^2 .* sc{1}.^2 + ja.^2 .* sc{2}.^2) ./ dt2;
    else
        P = 0;  Q = 0;  Rr = 0;
        for k = 1:n
            wk = 1 ./ sc{k}.^2;  gx = s.gx(:,:,k);  gy = s.gy(:,:,k);
            P = P + wk .* gx.^2;  Q = Q + wk .* gx .* gy;  Rr = Rr + wk .* gy.^2;
        end
        dW = max(P .* Rr - Q.^2, eps('single'));
        s.sSx = kS(1) * sqrt(Rr ./ dW);
        s.sSy = kS(2) * sqrt(P ./ dW);
    end
    s.good = s.pass & max(s.sSx, s.sSy) <= o.max_slope_sigma;
end

function v = smin_rg(gx1, gy1, gx2, gy2)
    % smallest singular value of J_RG = [gx1 gy1; gx2 gy2]: sqrt of the smaller eigenvalue of J'J
    p = gx1.^2 + gx2.^2;  q = gx1 .* gy1 + gx2 .* gy2;  r = gy1.^2 + gy2.^2;
    v = sqrt(max((p + r - sqrt((p - r).^2 + 4 * q.^2)) / 2, 0));
end

function [ok, dtm] = conditioning(a, b, c, d, det_min, cond_max)
    dtm = a .* d - b .* c;
    f2 = a.^2 + b.^2 + c.^2 + d.^2;
    s1sq = (f2 + sqrt(max(f2.^2 - 4 * dtm.^2, 0))) / 2;    % largest singular value squared
    cnd = s1sq ./ max(abs(dtm), eps('single'));             % = sigma1 / sigma2
    ok = abs(dtm) >= det_min & cnd <= cond_max & isfinite(dtm);
    dtm(~ok) = 1;                                           % avoid 0/0 in the step; step is zeroed by ok
end

function C = prep(img, o)
    img = single(img);
    R = img(:,:,1);  G = img(:,:,2);  B = img(:,:,3);
    if strcmp(o.norm, 'sum')
        S = max(R + G + B, 1);
        C = cat(3, R ./ S, G ./ S);
    elseif o.nch == 3
        C = cat(3, R, G, B) / 255;
    else
        C = cat(3, R, G) / 255;
    end
    C = imgaussfilt(C, o.sigma);
end

function o = getopts(s)
    d = struct('colour_norm', 'none', 'colour_sigma', 12, 'max_disp_px', 600, 'colour_n_iter', 15, ...
               'det_min_frac', 0.1, 'cond_max', 10, 'res_max', 0.05, 'conv_px', 2, ...
               'roi', [], 'target_cols', [], 'drift_corr', 'none', ...
               'dye_mask', false, 'dye_thresh', [0.010 0.004], 'dye_dilate_px', 6, 'dye_win', 101, 'dye_model', 'none', 'dark_mask', [], 'noise', [], 'max_slope_sigma', 0.02, ...
               'init_dx', [], 'init_dy', [], 'colour_channels', 'RG', 'hybrid_frac', 0.3);
    for f = fieldnames(d).'
        if isfield(s, f{1}), d.(f{1}) = s.(f{1}); end
    end
    o = struct('norm', d.colour_norm, 'sigma', d.colour_sigma, 'max_disp_px', d.max_disp_px, ...
               'n_iter', d.colour_n_iter, 'det_min_frac', d.det_min_frac, 'cond_max', d.cond_max, ...
               'res_max', d.res_max, 'conv_px', d.conv_px, 'roi', d.roi, 'target_cols', d.target_cols, 'drift', d.drift_corr, ...
               'dye_mask', d.dye_mask, 'dye_thresh', d.dye_thresh, 'dye_dilate', d.dye_dilate_px, 'dye_win', d.dye_win, 'dye_model', d.dye_model, 'dark_mask', d.dark_mask, 'noise', d.noise, 'max_slope_sigma', d.max_slope_sigma, ...
               'init_dx', d.init_dx, 'init_dy', d.init_dy, 'channels', upper(d.colour_channels), 'hybrid_frac', d.hybrid_frac);
    switch o.channels
        case 'RG',     o.nch = 2;
        case {'RGB', 'HYBRID'}, o.nch = 3;
        otherwise,     error('cisg_colour_displacement:channels', 'colour_channels must be ''RG'', ''RGB'' or ''hybrid'', not ''%s''', d.colour_channels);
    end
    o.hybrid = strcmp(o.channels, 'HYBRID');
    if o.nch == 3 && strcmp(o.norm, 'sum')
        error('cisg_colour_displacement:channels', 'colour_channels ''%s'' needs colour_norm ''none'' (R/S, G/S, B/S sum to 1)', d.colour_channels);
    end
    if o.nch == 3 && strcmp(o.dye_model, 'joint')
        error('cisg_colour_displacement:channels', 'colour_channels ''%s'' cannot be combined with dye_model ''joint'' (blue is used for the dye)', d.colour_channels);
    end
end
