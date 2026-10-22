function [dx, dy, cmap, info] = cisg_dye_unmix(ref_image, obs_image, setup, dx0, dy0)
%CISG_DYE_UNMIX Joint per-pixel solve for the pattern shift and a dye amount.
%
%   [dx, dy, cmap, info] = cisg_dye_unmix(ref_image, obs_image, setup, dx0, dy0)
%
% Wave motion only moves the pattern, C_obs = C_ref(x + delta).  Rhodamine 6G pulls the colour
% off that manifold; in run 95 the perturbation in the dye stripes is, to within the scatter,
% a constant direction d in (R, G, B) (R up, B down, G ~ 0: emission minus absorption), so
%
%   C_obs(x) = C_ref(x + delta) + c d,    c >= 0 (dye amount),   |d| = 1,
%
% solved for (delta_x, delta_y, c) from all three channels by Newton iterations with the 3x3
% Jacobian [dC_ref/dx, dC_ref/dy, d], starting from the 2-channel solution (dx0, dy0; NaN -> 0)
% and c = 0.  A frame-wide blue offset (lighting / exposure, ~2-4 counts even on clean frames)
% is removed first (robust median of the blue residual of the 2-channel solution), otherwise it
% would be read as dye.  Pixels where d lies (nearly) in the plane of the shift columns
% (|det J| / prod(column norms) < setup.unmix_min_sine) cannot be separated: NaN.
%
% dx0, dy0: [H x W] maps (full image size, as returned by cisg_colour_displacement), solved inside setup.roi.
% Hybrid output: the 2-channel shift (and its mask) is kept everywhere; the joint shift replaces it only where
% the 3x3 is well conditioned and the smoothed dye amount (imgaussfilt sigma unmix_cs_sigma = [4 40] px: stripes are
% streamwise) is >= unmix_c_min (0.015; clean frame 3000: applied on ~5% of pixels, Sx change 5% of rms).  info.applied marks them.
% setup fields (defaults): colour_sigma, roi, dye_dir ([0.50 0.07 -0.86]), unmix_iter 8,
%   unmix_min_sine 0.15, max_disp_px, res_max 0.05 (3-channel residual limit)
% cmap: [H x W] dye amount (NaN where not solved), info: angle statistics, mask fraction, offset.

    o = struct('dye_dir', [0.50 0.07 -0.86], 'unmix_iter', 8, 'unmix_min_sine', 0.15, 'res_max', 0.05, ...
               'unmix_cs_sigma', [4 40], 'unmix_c_min', 0.015);
    for f = fieldnames(o).'
        if isfield(setup, f{1}), o.(f{1}) = setup.(f{1}); end
    end
    d = single(o.dye_dir(:)' / norm(o.dye_dir));
    sig = setup.colour_sigma;
    [H, W, ~] = size(ref_image);
    r = setup.roi;  r1 = r.row_min; r2 = r.row_max; c1 = r.col_min; c2 = r.col_max;

    Cr = imgaussfilt(single(ref_image) / 255, sig);
    gx = zeros(size(Cr), 'single');  gy = gx;
    for k = 1:3, [gx(:,:,k), gy(:,:,k)] = gradient(Cr(:,:,k)); end
    pad = ceil(4 * sig);
    pr = max(1, r1 - pad):min(H, r2 + pad);  pc = max(1, c1 - pad):min(W, c2 + pad);
    Co = imgaussfilt(single(obs_image(pr, pc, :)) / 255, sig);
    Co = Co(r1 - pr(1) + 1 : r2 - pr(1) + 1, c1 - pc(1) + 1 : c2 - pc(1) + 1, :);

    ux = dx0(r1:r2, c1:c2);  uy = dy0(r1:r2, c1:c2);
    solved0 = isfinite(ux) & isfinite(uy);
    ux(~solved0) = 0;  uy(~solved0) = 0;
    [X, Y] = meshgrid(single(c1:c2), single(r1:r2));
    xs = min(max(X + ux, 1), W);  ys = min(max(Y + uy, 1), H);
    % frame-wide blue offset from the pixels the 2-channel solve accepted
    eB0 = Co(:,:,3) - interp2(Cr(:,:,3), xs, ys);
    boff = median(eB0(solved0), 'omitnan');
    Co(:,:,3) = Co(:,:,3) - boff;

    cc = zeros(size(ux), 'single');
    for it = 1:o.unmix_iter
        xs = min(max(X + ux, 1), W);  ys = min(max(Y + uy, 1), H);
        for k = 1:3
            M = interp2(Cr(:,:,k), xs, ys) + cc * d(k);
            E(:,:,k) = Co(:,:,k) - M;                                       %#ok<AGROW>
            Jx(:,:,k) = interp2(gx(:,:,k), xs, ys);                         %#ok<AGROW>
            Jy(:,:,k) = interp2(gy(:,:,k), xs, ys);                         %#ok<AGROW>
        end
        [du, dv, dc, ~] = solve3(Jx, Jy, d, E);
        ux = min(max(ux + du, -setup.max_disp_px), setup.max_disp_px);
        uy = min(max(uy + dv, -setup.max_disp_px), setup.max_disp_px);
        cc = max(cc + dc, 0);                                               % dye amount >= 0
    end
    xs = min(max(X + ux, 1), W);  ys = min(max(Y + uy, 1), H);
    for k = 1:3
        E(:,:,k) = Co(:,:,k) - interp2(Cr(:,:,k), xs, ys) - cc * d(k);
        Jx(:,:,k) = interp2(gx(:,:,k), xs, ys);  Jy(:,:,k) = interp2(gy(:,:,k), xs, ys);
    end
    [~, ~, ~, sine] = solve3(Jx, Jy, d, E);
    res = sqrt(sum(E.^2, 3));
    inb = isfinite(ux) & (X + ux >= 1) & (X + ux <= W) & (Y + uy >= 1) & (Y + uy <= H);
    good3 = sine >= o.unmix_min_sine & res <= o.res_max & inb;           % the 3x3 solve is trustworthy here

    % Hybrid: keep the 2-channel shift (and its own mask) everywhere; take the joint shift only where the
    % solve is well conditioned AND the (smoothed) dye amount is significant.  Nothing is NaN-ed by the joint step.
    cs = imgaussfilt(cc, o.unmix_cs_sigma);
    apply = good3 & solved0 & cs >= o.unmix_c_min;
    ux2 = dx0(r1:r2, c1:c2);  uy2 = dy0(r1:r2, c1:c2);
    ux2(apply) = ux(apply);  uy2(apply) = uy(apply);
    cc(~good3) = NaN;

    dx = nan(H, W, 'single');  dy = dx;  cmap = dx;
    dx(r1:r2, c1:c2) = ux2;  dy(r1:r2, c1:c2) = uy2;  cmap(r1:r2, c1:c2) = cc;
    applied = false(H, W);  applied(r1:r2, c1:c2) = apply;
    % angle between d and the plane of the shift columns (0 = d inside the shift subspace)
    ang = asind(min(sine, 1));
    info = struct('mask_fraction', mean(isnan(ux2(:))), 'applied', applied, 'applied_fraction', mean(apply(:)), ...
                  'joint_only_mask_fraction', mean(~good3(:)), 'frac_lowsine', mean(sine(:) < o.unmix_min_sine), ...
                  'frac_res', mean(res(:) > o.res_max), 'frac_oob', mean(~inb(:)), 'c_smooth', cs, ...
                  'blue_offset', boff, 'angle_deg', ang(isfinite(ang)), 'min_sine', o.unmix_min_sine, 'dye_dir', d);
end

function [du, dv, dc, sine] = solve3(Jx, Jy, d, E)
    % Solve [Jx Jy d] * [du dv dc]' = E per pixel (Cramer); sine = |det| / (|Jx||Jy||d|) = sin of the volume angle
    ax = Jx(:,:,1); ay = Jx(:,:,2); az = Jx(:,:,3);
    bx = Jy(:,:,1); by = Jy(:,:,2); bz = Jy(:,:,3);
    dxx = d(1); dyy = d(2); dzz = d(3);
    D = ax .* (by * dzz - bz * dyy) - bx .* (ay * dzz - az * dyy) + dxx * (ay .* bz - az .* by);
    nrm = sqrt(ax.^2 + ay.^2 + az.^2) .* sqrt(bx.^2 + by.^2 + bz.^2);
    sine = abs(D) ./ max(nrm, eps('single'));
    Dm = D;  Dm(abs(Dm) < eps('single')) = eps('single');
    e1 = E(:,:,1); e2 = E(:,:,2); e3 = E(:,:,3);
    du = (e1 .* (by * dzz - bz * dyy) - bx .* (e2 * dzz - e3 * dyy) + dxx * (e2 .* bz - e3 .* by)) ./ Dm;
    dv = (ax .* (e2 * dzz - e3 * dyy) - e1 .* (ay * dzz - az * dyy) + dxx * (ay .* e3 - az .* e2)) ./ Dm;
    dc = (ax .* (by .* e3 - bz .* e2) - bx .* (ay .* e3 - az .* e2) + e1 .* (ay .* bz - az .* by)) ./ Dm;
    bad = sine < 0.02;  du(bad) = 0;  dv(bad) = 0;  dc(bad) = 0;
end
