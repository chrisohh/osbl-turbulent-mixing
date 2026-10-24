function [roi, info] = cisg_auto_roi(ref_image, setup, opts)
%CISG_AUTO_ROI ROI that avoids the regions that are already unusable in the REFERENCE image.
%
%   [roi, info] = cisg_auto_roi(ref_image, setup)
%   [roi, info] = cisg_auto_roi(ref_image, setup, opts)
%
% A pixel of the reference is usable when the colour map there can be inverted: bright enough
% (not in a wall shadow / outside the target), inside the target columns (setup.target_cols),
% and a well-conditioned Jacobian of [R G] (|det J| >= det_min_frac * median, cond(J) <= cond_max,
% the same tests as cisg_colour_displacement).  The ROI is the bounding box of the usable pixels,
% shrunk edge by edge (the edge row / column with the lowest usable fraction first) until every
% edge has at least opts.min_frac usable pixels, then inset by opts.margin_px (pattern shifts of up
% to max_disp_px must still land on usable colours).
%
% opts fields (defaults): min_frac 0.95, bright_min 0.15 (mean RGB / 255), margin_px 0,
%   scale 4 (analysis at 1/scale resolution)
% roi in full-resolution pixels (row_min row_max col_min col_max); info.usable = the logical map
% (analysis resolution), info.scale.

    def = struct('min_frac', 0.95, 'bright_min', 0.15, 'margin_px', 0, 'scale', 4);
    if nargin < 3, opts = struct(); end
    for f = fieldnames(def).'
        if ~isfield(opts, f{1}), opts.(f{1}) = def.(f{1}); end
    end
    s = opts.scale;
    Ir = single(imresize(ref_image, 1 / s, 'box')) / 255;
    sg = max(1, setup.colour_sigma / s);
    C1 = imgaussfilt(Ir(:,:,1), sg);  C2 = imgaussfilt(Ir(:,:,2), sg);
    S  = imgaussfilt(mean(Ir, 3), sg);
    [ax, ay] = gradient(C1);  [bx, by] = gradient(C2);
    dt = ax .* by - ay .* bx;
    f2 = ax.^2 + ay.^2 + bx.^2 + by.^2;
    s1 = (f2 + sqrt(max(f2.^2 - 4 * dt.^2, 0))) / 2;
    cnd = s1 ./ max(abs(dt), eps('single'));
    [h, w] = size(S);
    cols = true(1, w);
    if isfield(setup, 'target_cols') && ~isempty(setup.target_cols)
        tc = (setup.target_cols - 0.5) / s + 0.5;
        cols = (1:w) >= tc(1) & (1:w) <= tc(2);
    end
    inside = repmat(cols, h, 1);
    frac = 0.1;  if isfield(setup, 'det_min_frac'), frac = setup.det_min_frac; end
    det_min = frac * median(abs(dt(inside)));
    usable = inside & S >= opts.bright_min & abs(dt) >= det_min & cnd <= 10;

    r1 = find(any(usable, 2), 1);  r2 = find(any(usable, 2), 1, 'last');
    c1 = find(any(usable, 1), 1);  c2 = find(any(usable, 1), 1, 'last');
    for it = 1:h + w
        U = usable(r1:r2, c1:c2);
        fr = [mean(U(1, :)) mean(U(end, :)) mean(U(:, 1)) mean(U(:, end))];   % top bottom left right
        [m, k] = min(fr);
        if m >= opts.min_frac || r2 - r1 < 10 || c2 - c1 < 10, break; end
        switch k
            case 1, r1 = r1 + 1;
            case 2, r2 = r2 - 1;
            case 3, c1 = c1 + 1;
            case 4, c2 = c2 - 1;
        end
    end
    mg = ceil(opts.margin_px / s);
    r1 = r1 + mg;  r2 = r2 - mg;  c1 = c1 + mg;  c2 = c2 - mg;
    [H0, W0] = deal(size(ref_image, 1), size(ref_image, 2));
    roi = struct('row_min', max(1, round((r1 - 1) * s + 1)), 'row_max', min(H0, round(r2 * s)), ...
                 'col_min', max(1, round((c1 - 1) * s + 1)), 'col_max', min(W0, round(c2 * s)));
    info = struct('usable', usable, 'dark', ~(inside & S >= opts.bright_min), 'scale', s, 'opts', opts);   % dark: outside the target or in shadow
end
