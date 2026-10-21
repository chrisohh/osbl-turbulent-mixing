function [dx, dy, info] = cisg_colour_displacement(ref_image, obs_image, setup)
%CISG_COLOUR_DISPLACEMENT Local colour -> position inversion of the CISG target.
%
%   [dx, dy, info] = cisg_colour_displacement(ref_image, obs_image, setup)
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
%   roi              []      struct row_min/row_max/col_min/col_max: only these pixels
%                            are solved (output is still full size, NaN elsewhere)
%   target_cols      []      [left right] columns of the target (median |det J| uses them)
%
% info: mask fraction (inside the solved region), thresholds, runtime.
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
    key = [H W o.sigma double(strcmp(o.norm, 'sum')) sum(double(ref_image(1:997:end)))];
    if isempty(cache) || ~isequal(cache.key, key)
        C = prep(ref_image, o);
        [g1x, g1y] = gradient(C(:,:,1));
        [g2x, g2y] = gradient(C(:,:,2));
        cache = struct('key', key, 'C1', C(:,:,1), 'C2', C(:,:,2), ...
                       'a', g1x, 'b', g1y, 'c', g2x, 'd', g2y);
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
        good = ok & hypot(e1, e2) <= o.res_max & step <= o.conv_px & isfinite(ux) & isfinite(uy) ...
               & (X + ux >= 1) & (X + ux <= W) & (Y + uy >= 1) & (Y + uy <= H);
        ux(~good) = NaN;  uy(~good) = NaN;
        dx(rs:re, c1:c2) = ux;  dy(rs:re, c1:c2) = uy;
        n_tot = n_tot + numel(good);  n_bad = n_bad + nnz(~good);
    end

    info = struct('mask_fraction', n_bad / n_tot, 'det_min', det_min, 'opts', o, 'drift_ab', drift_ab, 'seconds', toc(t0));
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
               'roi', [], 'target_cols', [], 'drift_corr', 'none');
    for f = fieldnames(d).'
        if isfield(s, f{1}), d.(f{1}) = s.(f{1}); end
    end
    o = struct('norm', d.colour_norm, 'sigma', d.colour_sigma, 'max_disp_px', d.max_disp_px, ...
               'n_iter', d.colour_n_iter, 'det_min_frac', d.det_min_frac, 'cond_max', d.cond_max, ...
               'res_max', d.res_max, 'conv_px', d.conv_px, 'roi', d.roi, 'target_cols', d.target_cols, 'drift', d.drift_corr);
end
