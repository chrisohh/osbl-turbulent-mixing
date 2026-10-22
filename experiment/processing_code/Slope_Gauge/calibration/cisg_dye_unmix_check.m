%% cisg_dye_unmix_check.m
% Feasibility of per-pixel dye unmixing (util/cisg_dye_unmix.m): signature of the dye in run 95,
% clean-frame noise cost, synthetic known-shift + known-dye test, and frame 1510.
% Writes calibration/cisg_dye_unmix_check.png (invisible figure).  Single MATLAB, no pool, no writes outside calibration/.
clear; clc;
here = fileparts(mfilename('fullpath'));
sg = fileparts(here);
addpath(here, fullfile(sg, 'util'), sg, fullfile(fileparts(sg), 'IR Camera', 'calibration'));

roi = struct('row_min', 793, 'row_max', 2536, 'col_min', 1033, 'col_max', 2940);     % 'auto' ROI, margin 100 px
R = cisg_setup_run(struct('roi_mode', roi));
[ref_b, st] = cisg_bin_setup(R.ref_image, R.setup, 2);
r = st.roi;  rg = @(A) A(r.row_min:r.row_max, r.col_min:r.col_max, :);
getobs = @(n) uint8(imresize(cisg_load_coreview(R.frame_file(n), R.img_size), 0.5, 'box'));
px2S = st.pattern_scale * st.cm_per_pixel_x / ((st.n_water - st.n_air) * st.L_eff_cm);
sig = st.colour_sigma;
[H, W, ~] = size(ref_b);

%% ---- 1. signature of the dye (frames with dye stripes) ----
Cr = imgaussfilt(single(ref_b) / 255, sig);
gx = zeros(size(Cr), 'single');  gy = gx;
for k = 1:3, [gx(:,:,k), gy(:,:,k)] = gradient(Cr(:,:,k)); end
[X, Y] = meshgrid(single(r.col_min:r.col_max), single(r.row_min:r.row_max));
Pall = [];  Call = [];  Ang = [];
for n = [1505 1510 1515]
    [~, ~, ~, D] = cisg_colour_displacement(ref_b, getobs(n), st);
    Co = rg(imgaussfilt(single(getobs(n)) / 255, sig));
    good = D.bits == 0;
    xs = min(max(X + D.dx, 1), W);  ys = min(max(Y + D.dy, 1), H);
    Cx = zeros([size(X) 3], 'single');  Jx = Cx;  Jy = Cx;
    for k = 1:3
        Cx(:,:,k) = interp2(Cr(:,:,k), xs, ys);  Jx(:,:,k) = interp2(gx(:,:,k), xs, ys);  Jy(:,:,k) = interp2(gy(:,:,k), xs, ys);
    end
    dC = Co - Cx;
    dxg = D.dx;  dxg(~good) = NaN;  dyg = D.dy;  dyg(~good) = NaN;
    rdx = median(dxg, 2, 'omitnan');  rdy = median(dyg, 2, 'omitnan');
    ddx = rdx - movmedian(rdx, 101, 'omitnan');  ddy = rdy - movmedian(rdy, 101, 'omitnan');   % row-coherent shift bias
    pert = dC;                                     % dye perturbation = residual + colour change from the biased shift
    for k = 1:2, pert(:,:,k) = dC(:,:,k) + Jx(:,:,k) .* ddx + Jy(:,:,k) .* ddy; end
    [dB, dS, valid] = rowdev(D);
    flag = imdilate(abs(dB) > 0.010 | abs(dS) > 0.004, ones(5, 1));
    sel = repmat(flag, 1, size(good, 2)) & good;  idx = find(sel);  idx = idx(1:max(1, round(numel(idx) / 15000)):end);
    P = reshape(pert, [], 3);  C = reshape(Cx, [], 3);  Pall = [Pall; P(idx, :)];  Call = [Call; C(idx, :)];   %#ok<AGROW>
    nx = reshape(cross(reshape(Jx, [], 3), reshape(Jy, [], 3)), [], 3);
    nx = nx ./ max(vecnorm(nx, 2, 2), eps);
end
mu = mean(Pall);  [~, S, V] = svd(Pall - mu, 'econ');  dir1 = V(:, 1)' * sign(-V(3, 1));
fprintf('1. dye signature (colour units 0..1; 1 count = 0.0039), %d pixels in dye rows of frames 1505/1510/1515\n', size(Pall, 1));
fprintf('   mean perturbation (R G B) = %s, std = %s\n', mat2str(mu, 3), mat2str(std(Pall), 3));
fprintf('   dominant direction d = %s (singular values %s)\n', mat2str(dir1, 3), mat2str(diag(S)' / sqrt(size(Pall, 1)), 3));
for k = 1:3
    A = [Call(:, k) ones(size(Call, 1), 1)];  th = A \ Pall(:, k);  res = Pall(:, k) - A * th;
    fprintf('   channel %d: perturbation = %.3f * C_ref %+.4f, R^2 %.2f\n', k, th(1), th(2), 1 - var(res) / var(Pall(:, k)));
end
fprintf('   -> perturbation is ~constant in C_ref (absorption term not resolved); model used: C_obs = C_ref(x+delta) + c d, d = %s\n', mat2str([0.50 0.07 -0.86]));

%% ---- 2./3a. clean frame 3000 and 3c. frame 1510 ----
for n = [3000 1510]
    obs = getobs(n);
    t1 = tic;  [dx2, dy2, i2] = cisg_colour_displacement(ref_b, obs, st);  t2 = toc(t1);
    t1 = tic;  [dxj, dyj, cm, ij] = cisg_dye_unmix(ref_b, obs, st, dx2, dy2);  tj = toc(t1);
    a2x = rg(dx2) * px2S;  a2y = rg(dy2) * px2S;  ajx = rg(dxj) * px2S;  ajy = rg(dyj) * px2S;  c = rg(cm);
    ok = isfinite(a2x) & isfinite(ajx);
    fprintf('\nframe %d: 2-channel %.1f s (masked %.1f%%), joint +%.1f s (masked %.1f%%), blue offset removed %.4f\n', ...
            n, t2, 100*i2.mask_fraction, tj, 100*ij.mask_fraction, ij.blue_offset);
    fprintf('   angle between d and the shift plane: median %.0f deg [p10 %.0f, p90 %.0f]\n', median(ij.angle_deg), prctile(ij.angle_deg, 10), prctile(ij.angle_deg, 90));
    fprintf('   c: median %.4f, p90 %.4f, p99 %.4f\n', median(c(ok)), prctile(c(ok), 90), prctile(c(ok), 99));
    fprintf('   rms Sx 2ch %.4f joint %.4f | Sy 2ch %.4f joint %.4f | rms difference Sx %.4f, Sy %.4f\n', ...
            rms(a2x(ok)), rms(ajx(ok)), rms(a2y(ok)), rms(ajy(ok)), rms(a2x(ok) - ajx(ok)), rms(a2y(ok) - ajy(ok)));
    fprintf('   var(row-mean): Sx 2ch %.5f joint %.5f | Sy 2ch %.5f joint %.5f\n', var(mean(a2x, 2, 'omitnan')), var(mean(ajx, 2, 'omitnan')), ...
            var(mean(a2y, 2, 'omitnan')), var(mean(ajy, 2, 'omitnan')));
    fprintf('   corr with c: Sx 2ch %.2f joint %.2f | Sy 2ch %.2f joint %.2f\n', corr(a2x(ok), c(ok)), corr(ajx(ok), c(ok)), corr(a2y(ok), c(ok)), corr(ajy(ok), c(ok)));
    if n == 1510, F = struct('a2x', a2x, 'a2y', a2y, 'ajx', ajx, 'ajy', ajy, 'c', c); end
end

%% ---- 3b. synthetic: known shift + known diffuse dye ----
rng(5);
mg = 60;  rows = r.row_min - mg : r.row_max + mg;  cols = r.col_min - mg : r.col_max + mg;
[Xb, Yb] = meshgrid(single(cols), single(rows));
fld = @(s) imgaussfilt(randn(numel(rows), numel(cols), 'single'), s);
nz = @(f) f / std(f(:));
tx = 90 * nz(fld(60));  ty = 55 * nz(fld(60));                    % smooth waves: rms 90 / 55 px
[gxx, gxy] = gradient(tx);  [gyx, gyy] = gradient(ty);
fprintf('\n3b. synthetic: rms delta %.0f/%.0f px, min det(I+grad delta) %.2f (fold if <= 0)\n', rms(tx(:)), rms(ty(:)), min(min((1 + gxx) .* (1 + gyy) - gxy .* gyx)));
g1 = nz(fld(80));  g2 = nz(fld(25));
c_true = 0.05 * (max(g1 - 0.6, 0) / max(max(g1(:)) - 0.6, eps) + 0.7 * exp(-(g2 / 0.12).^2) .* (nz(fld(150)) > 0.3));   % patches + filaments
c_true = min(c_true, 0.06);
d_true = [0.50 0.07 -0.86] / norm([0.50 0.07 -0.86]);
for variant = 1:2
    if variant == 2, d_true = [0.35 -0.15 -0.92]; d_true = d_true / norm(d_true); end         % wrong signature: robustness
    obs = ref_b;
    for k = 1:3
        v = interp2(double(ref_b(:,:,k)), double(Xb + tx), double(Yb + ty), 'linear') / 255 + double(c_true) * d_true(k);
        obs(rows, cols, k) = uint8(round(255 * min(max(v, 0), 1)));
    end
    [dx2, dy2] = cisg_colour_displacement(ref_b, obs, st);
    [dxj, dyj, cm] = cisg_dye_unmix(ref_b, obs, st, dx2, dy2);
    cut = @(A) A(mg + 1 : end - mg, mg + 1 : end - mg);
    TX = cut(tx) * px2S;  TY = cut(ty) * px2S;  CT = cut(c_true);
    a2x = rg(dx2) * px2S;  a2y = rg(dy2) * px2S;  ajx = rg(dxj) * px2S;  ajy = rg(dyj) * px2S;  cc = rg(cm);
    ok = isfinite(a2x) & isfinite(ajx);
    gain = @(t, a) [t(:) * 0 + 1, t(:)] \ a(:);
    sx2 = polyfit(TX(ok), a2x(ok), 1);  sxj = polyfit(TX(ok), ajx(ok), 1);  sy2 = polyfit(TY(ok), a2y(ok), 1);  syj = polyfit(TY(ok), ajy(ok), 1);
    fprintf('   dye direction %s (%s): c_true max %.3f, area c>0.01: %.0f%%\n', mat2str(d_true, 2), iff(variant == 1, 'as assumed', 'WRONG'), max(c_true(:)), 100*mean(c_true(:) > 0.01));
    fprintf('     Sx 2ch: gain %.2f rms err %.4f bias %+.4f | joint: gain %.2f rms err %.4f bias %+.4f\n', sx2(1), rms(a2x(ok) - TX(ok)), mean(a2x(ok) - TX(ok)), sxj(1), rms(ajx(ok) - TX(ok)), mean(ajx(ok) - TX(ok)));
    fprintf('     Sy 2ch: gain %.2f rms err %.4f bias %+.4f | joint: gain %.2f rms err %.4f bias %+.4f\n', sy2(1), rms(a2y(ok) - TY(ok)), mean(a2y(ok) - TY(ok)), syj(1), rms(ajy(ok) - TY(ok)), mean(ajy(ok) - TY(ok)));
    pm = ok & CT > 0.01;                                              % inside the dye patches only
    fprintf('     inside dye patches (c>0.01, %.0f%% of ok px): rms err Sx 2ch %.4f joint %.4f | Sy 2ch %.4f joint %.4f ; bias Sx %+.4f/%+.4f Sy %+.4f/%+.4f\n', ...
            100*mean(pm(ok)), rms(a2x(pm) - TX(pm)), rms(ajx(pm) - TX(pm)), rms(a2y(pm) - TY(pm)), rms(ajy(pm) - TY(pm)), ...
            mean(a2x(pm) - TX(pm)), mean(ajx(pm) - TX(pm)), mean(a2y(pm) - TY(pm)), mean(ajy(pm) - TY(pm)));
    okc = ok & isfinite(cc);  pc = polyfit(CT(okc), cc(okc), 1);
    fprintf('     recovered c vs true c: corr %.2f, gain %.2f, rms err %.4f (colour units)\n', corr(CT(okc), cc(okc)), pc(1), rms(cc(okc) - CT(okc)));
    if variant == 1, SY = struct('TX', TX, 'a2x', a2x, 'ajx', ajx, 'CT', CT, 'cc', cc); end
end

%% ---- figure ----
q = 8;  d8 = @(A) A(1:q:end, 1:q:end);
fig = figure('Visible', 'off', 'Position', [50 50 1600 900]);
subplot(2,4,1); imagesc(d8(F.a2x)); axis image; caxis([-0.3 0.3]); colorbar; title('S_x 2-channel, frame 1510');
subplot(2,4,2); imagesc(d8(F.ajx)); axis image; caxis([-0.3 0.3]); colorbar; title('S_x joint');
subplot(2,4,3); imagesc(d8(F.a2y)); axis image; caxis([-0.3 0.3]); colorbar; title('S_y 2-channel');
subplot(2,4,4); imagesc(d8(F.ajy)); axis image; caxis([-0.3 0.3]); colorbar; title('S_y joint');
subplot(2,4,5); imagesc(d8(F.c)); axis image; colorbar; title('dye amount c (frame 1510)');
subplot(2,4,6); imagesc(d8(SY.CT)); axis image; colorbar; title('synthetic: true c');
subplot(2,4,7); imagesc(d8(SY.cc)); axis image; colorbar; title('synthetic: recovered c');
subplot(2,4,8); s = 1:97:numel(SY.cc); plot(SY.CT(s), SY.cc(s), '.'); xlabel('true c'); ylabel('recovered c'); grid on; title('c recovery');
exportgraphics(fig, fullfile(here, 'cisg_dye_unmix_check.png'), 'Resolution', 90); close(fig);

function s = iff(c, a, b), if c, s = a; else, s = b; end, end
function [dB, dS, valid] = rowdev(D)
    g = D.bits == 0;  eB = D.eB;  eB(~g) = NaN;  eS = D.eS;  eS(~g) = NaN;
    rB = median(eB, 2, 'omitnan');  rS = median(eS, 2, 'omitnan');
    valid = mean(g, 2) > 0.2;  rB(~valid) = NaN;  rS(~valid) = NaN;
    dB = rB - movmedian(rB, 101, 'omitnan');  dS = rS - movmedian(rS, 101, 'omitnan');
    dB(~valid) = 0;  dS(~valid) = 0;
end
