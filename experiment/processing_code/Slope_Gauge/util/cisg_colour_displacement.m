function [dx, dy, info, D] = cisg_colour_displacement(ref_image, obs_image, setup)
%CISG_COLOUR_DISPLACEMENT Local colour -> position inversion of the CISG target.
%
%   [dx, dy, info] = cisg_colour_displacement(ref_image, obs_image, setup)
%   [dx, dy, info, D] = cisg_colour_displacement(...)        % D: per-pixel diagnostics (ROI-sized)
%
% For every pixel x of the observed image, find the shift delta (pixels) such
% that the reference colour at x + delta equals the observed colour at x,
%
%   C_ref(x + delta) = C_obs(x),   C = [R G] (smoothed, optionally R/(R+G+B), G/(R+G+B)).
%
% Newton / Lucas-Kanade iterations with the local Jacobian of the reference,
%   J = [dR/dx dR/dy; dG/dx dG/dy],   delta <- delta + J^-1 (C_obs(x) - C_ref(x + delta)),
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
%   colour_norm      'none'  'none' | 'sum' (R/(R+G+B), G/(R+G+B): cancels brightness changes)
%   colour_sigma     12      Gaussian smoothing of both images (px)
%   max_disp_px      600     clamp on |delta|
%   colour_n_iter    15
%   det_min_frac     0.1     mask |det J| < frac * median |det J| of the reference
%   cond_max         10      mask cond(J) above this
%   res_max          0.05    mask final colour residual (|C_obs - C_ref|, 0..1 units)
%   conv_px          2        mask last Newton step above this (px)
%   drift_corr       'none'  'gain': per-channel robust fit obs = a ref + b over the ROI, obs <- (obs - b)/a (lighting drift)
%   dye_mask         false   mask cross-stream dye stripes: rows whose median blue / brightness residual (obs vs ref at x+delta,
%                            0..1), minus its running median over dye_win rows (101), exceeds dye_thresh = [tB tS] (default
%                            [0.010 0.004]) are masked, +- dye_dilate_px rows (6).  R6G dims B and brightness, which [R G] ignores
%   dye_model        'none'  'joint': solve (shift, dye amount c) per pixel from R,G,B (cisg_dye_unmix.m); c in info.dye_c.
%                            Noisier than 'none' (adds ~8% rms to Sx on clean frames), see calibration/cisg_dye_unmix_check.m
%   roi              []      struct row_min/row_max/col_min/col_max: only these pixels
%                            are solved (output is still full size, NaN elsewhere)
%   target_cols      []      [left right] columns of the target (median |det J| uses them)
%
% info: mask_fraction (all reasons), dye_fraction (pixels masked by the dye test), thresholds, runtime.
% D (4th output): ROI-sized bits (1 x+delta out of bounds, 2 displacement clamped, 3 ill-conditioned J,
%   4 not converged, 5 colour residual, 6 dye), dx, dy (unmasked), eB, eS (blue / brightness residuals, 0..1).
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
    key = [H W o.sigma double(strcmp(o.norm, 'sum')) sum(double(ref_image(1:997:end)))];
    if isempty(cache) || ~isequal(cache.key, key) || (need_bs && ~isfield(cache, 'CB'))
        C = prep(ref_image, o);
        [g1x, g1y] = gradient(C(:,:,1));
        [g2x, g2y] = gradient(C(:,:,2));
        cache = struct('key', key, 'C1', C(:,:,1), 'C2', C(:,:,2), ...
                       'a', g1x, 'b', g1y, 'c', g2x, 'd', g2y);
        if need_bs
            [cache.CB, cache.CS] = blue_bright(ref_image, o);
        end
    end
    C1 = cache.C1; C2 = cache.C2; ga = cache.a; gb = cache.b; gc = cache.c; gd = cache.d;
    cs = 1:W;  if ~isempty(o.target_cols), cs = max(1, round(o.target_cols(1))):min(W, round(o.target_cols(2))); end
    dd = abs(ga(1:16:end, cs(1:16:end)) .* gd(1:16:end, cs(1:16:end)) - gb(1:16:end, cs(1:16:end)) .* gc(1:16:end, cs(1:16:end)));
    det_min = o.det_min_frac * median(dd(:));

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
    end

    %% Lighting-drift correction: obs ~ a * ref + b per channel (trimmed fit, waves do not bias it)
    drift_ab = [1 0; 1 0];
    if strcmp(o.drift, 'gain')
        for k = 1:2
            Rk = cache.(sprintf('C%d', k))(r1:r2, c1:c2);
            rcs = true(1, c2 - c1 + 1);
            if ~isempty(o.target_cols), rcs = (c1:c2) >= o.target_cols(1) + 100 & (c1:c2) <= o.target_cols(2) - 100; end
            xv = Rk(1:6:end, rcs(1:1:end)); xv = xv(:, 1:6:end); yv = Co(1:6:end, :, k); yv = yv(:, rcs); yv = yv(:, 1:6:end);
            xv = double(xv(:)); yv = double(yv(:));
            keep = true(size(xv));
            for it = 1:6
                p = polyfit(xv(keep), yv(keep), 1);
                res = abs(yv - polyval(p, xv));
                keep = res <= quantile(res, 0.7);
            end
            drift_ab(k, :) = p;                                % [a b]
            Co(:, :, k) = (Co(:, :, k) - p(2)) / p(1);
        end
    end

    %% Newton iterations (row chunks)
    dx = nan(H, W, 'single');  dy = nan(H, W, 'single');
    chunk = 256;
    n_tot = 0;  n_bad = 0;
    for rs = r1:chunk:r2
        re = min(rs + chunk - 1, r2);
        [X, Y] = meshgrid(single(c1:c2), single(rs:re));
        O1 = Co(rs - r1 + 1 : re - r1 + 1, :, 1);
        O2 = Co(rs - r1 + 1 : re - r1 + 1, :, 2);
        ux = zeros(size(X), 'single');  uy = ux;  step = ux;
        for it = 1:o.n_iter
            xs = min(max(X + ux, 1), W);  ys = min(max(Y + uy, 1), H);
            e1 = O1 - interp2(C1, xs, ys, 'linear');
            e2 = O2 - interp2(C2, xs, ys, 'linear');
            a = interp2(ga, xs, ys, 'linear');  b = interp2(gb, xs, ys, 'linear');
            c = interp2(gc, xs, ys, 'linear');  d = interp2(gd, xs, ys, 'linear');
            [ok, dtm] = conditioning(a, b, c, d, det_min, o.cond_max);
            sx = (d .* e1 - b .* e2) ./ dtm;
            sy = (a .* e2 - c .* e1) ./ dtm;
            sx(~ok) = 0;  sy(~ok) = 0;
            ux = min(max(ux + sx, -o.max_disp_px), o.max_disp_px);
            uy = min(max(uy + sy, -o.max_disp_px), o.max_disp_px);
            step = hypot(sx, sy);
        end
        % final state: residual and conditioning at x + delta
        xs = min(max(X + ux, 1), W);  ys = min(max(Y + uy, 1), H);
        e1 = O1 - interp2(C1, xs, ys, 'linear');
        e2 = O2 - interp2(C2, xs, ys, 'linear');
        [ok, ~] = conditioning(interp2(ga, xs, ys, 'linear'), interp2(gb, xs, ys, 'linear'), ...
                               interp2(gc, xs, ys, 'linear'), interp2(gd, xs, ys, 'linear'), det_min, o.cond_max);
        oob = ~((X + ux >= 1) & (X + ux <= W) & (Y + uy >= 1) & (Y + uy <= H));
        good = ok & hypot(e1, e2) <= o.res_max & step <= o.conv_px & isfinite(ux) & isfinite(uy) & ~oob;
        if need_bs
            ii = rs - r1 + 1 : re - r1 + 1;
            eBm(ii, :) = OB(ii, :) - interp2(cache.CB, xs, ys, 'linear');
            eSm(ii, :) = OS(ii, :) - interp2(cache.CS, xs, ys, 'linear');
            dxm(ii, :) = ux;  dym(ii, :) = uy;
            bits(ii, :) = uint8(oob) + 2 * uint8(abs(ux) >= o.max_disp_px | abs(uy) >= o.max_disp_px) + 4 * uint8(~ok) ...
                          + 8 * uint8(step > o.conv_px) + 16 * uint8(hypot(e1, e2) > o.res_max);
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

    dye_c = [];
    if strcmp(o.dye_model, 'joint')                          % per-pixel (shift, dye amount) solve from all 3 channels
        s_u = setup;  s_u.roi = struct('row_min', r1, 'row_max', r2, 'col_min', c1, 'col_max', c2);
        [dx, dy, dye_c] = cisg_dye_unmix(ref_image, obs_image, s_u, dx, dy);
        n_bad = nnz(isnan(dx(r1:r2, c1:c2)));
    end

    info = struct('mask_fraction', n_bad / n_tot, 'dye_fraction', dye_fraction, 'det_min', det_min, 'opts', o, ...
                  'drift_ab', drift_ab, 'dye_c', dye_c, 'seconds', toc(t0));
    if nargout > 3
        D = struct('bits', bits, 'dx', dxm, 'dy', dym, 'eB', eBm, 'eS', eSm);
        D.dye_rows = dye_rows;  if o.dye_mask, D.bits(dye_rows, :) = D.bits(dye_rows, :) + 32; end
    end
end

function [CB, CS] = blue_bright(img, o)
    img = single(img);
    CB = imgaussfilt(img(:,:,3) / 255, o.sigma);
    CS = imgaussfilt(mean(img, 3) / 255, o.sigma);
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
    else
        C = cat(3, R, G) / 255;
    end
    C = imgaussfilt(C, o.sigma);
end

function o = getopts(s)
    d = struct('colour_norm', 'none', 'colour_sigma', 12, 'max_disp_px', 600, 'colour_n_iter', 15, ...
               'det_min_frac', 0.1, 'cond_max', 10, 'res_max', 0.05, 'conv_px', 2, ...
               'roi', [], 'target_cols', [], 'drift_corr', 'none', ...
               'dye_mask', false, 'dye_thresh', [0.010 0.004], 'dye_dilate_px', 6, 'dye_win', 101, 'dye_model', 'none');
    for f = fieldnames(d).'
        if isfield(s, f{1}), d.(f{1}) = s.(f{1}); end
    end
    o = struct('norm', d.colour_norm, 'sigma', d.colour_sigma, 'max_disp_px', d.max_disp_px, ...
               'n_iter', d.colour_n_iter, 'det_min_frac', d.det_min_frac, 'cond_max', d.cond_max, ...
               'res_max', d.res_max, 'conv_px', d.conv_px, 'roi', d.roi, 'target_cols', d.target_cols, 'drift', d.drift_corr, ...
               'dye_mask', d.dye_mask, 'dye_thresh', d.dye_thresh, 'dye_dilate', d.dye_dilate_px, 'dye_win', d.dye_win, 'dye_model', d.dye_model);
end
