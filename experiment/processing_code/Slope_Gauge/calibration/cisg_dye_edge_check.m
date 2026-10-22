%% cisg_dye_edge_check.m
% (1) R6G dye stripes: signature in the blue / brightness residuals and the row-wise dye mask
% (2) why the edges of the ROI are masked: breakdown by reason, and an ROI from the reference
% Writes calibration/dye_mask_check.png and edge_mask_reasons.png (invisible figures).
clear; clc;
here = fileparts(mfilename('fullpath'));
sg = fileparts(here);
addpath(here, fullfile(sg, 'util'), sg, fullfile(fileparts(sg), 'IR Camera', 'calibration'));

roi = struct('row_min', 18, 'row_max', 2983, 'col_min', 410, 'col_max', 3714);   % the ROI of the real run (full res)
bin = 2;
R = cisg_setup_run(struct('roi_mode', roi));
[ref_b, st] = cisg_bin_setup(R.ref_image, R.setup, bin);
r = st.roi;  rg = @(A) A(r.row_min:r.row_max, r.col_min:r.col_max, :);
getobs = @(n) uint8(imresize(cisg_load_coreview(R.frame_file(n), R.img_size), 1/bin, 'box'));

%% ---- 1. dye ----
dye_frames = [1500 1505 1510 1515 1520 1530];  clean_frames = [1000 2000 2500 3000 3100];
st_d = st;  st_d.dye_mask = true;
fprintf('dye test: thresholds [tB tS] = [0.010 0.004], dilate 6 rows, window 101 rows\n');
for n = [dye_frames clean_frames]
    [~, ~, info, D] = cisg_colour_displacement(ref_b, getobs(n), st_d);
    rows = find(D.dye_rows);
    if ismember(n, dye_frames), tag = 'dye?'; else, tag = 'clean'; end
    fprintf('frame %4d (%s): dye-masked %.1f%% of ROI, total masked %.1f%%, %d flagged rows %s\n', n, tag, ...
            100*info.dye_fraction, 100*info.mask_fraction, numel(rows), rng(rows));
end
% signature of frame 1510 vs 3000 (dye off): row medians, per-pixel distributions, separability
[~, ~, ~, D1] = cisg_colour_displacement(ref_b, getobs(1510), st);
[~, ~, ~, D3] = cisg_colour_displacement(ref_b, getobs(3000), st);
[dB1, dS1, v1] = rowdev(D1, 101);  [dB3, dS3, v3] = rowdev(D3, 101);
flag = abs(dB1) > 0.010 | abs(dS1) > 0.004;                  % rows flagged on 1510
fprintf('\nframe 1510: %d rows flagged (%s)\n', nnz(flag), rng(find(flag)));
pf = valid_px(D1, flag);  pc = valid_px(D1, ~flag & v1);  p3 = valid_px(D3, v3);
fprintf('per-pixel residuals (0..1) median [p10 p90]:\n');
fprintf('  eB: dye rows %.4f [%.4f %.4f] | 1510 other rows %.4f [%.4f %.4f] | frame 3000 %.4f [%.4f %.4f]\n', stats(pf.eB), stats(pc.eB), stats(p3.eB));
fprintf('  eS: dye rows %.4f [%.4f %.4f] | 1510 other rows %.4f [%.4f %.4f] | frame 3000 %.4f [%.4f %.4f]\n', stats(pf.eS), stats(pc.eS), stats(p3.eS));
fprintf('separability (AUC, dye-row pixels vs clean pixels of 3000, 0.5 = none): per-pixel eB %.2f, eS %.2f\n', auc(-pf.eB, -p3.eB), auc(-pf.eS, -p3.eS));
fprintf('row-median deviation: dye rows |dB| min %.4f vs frame 3000 max |dB| %.4f ; |dS| min %.4f vs %.4f\n', ...
        min(abs(dB1(flag))), max(abs(dB3(v3))), min(abs(dS1(flag))), max(abs(dS3(v3))));
fprintf('false-positive rows on frame 3000: %d of %d valid rows\n', nnz(v3 & (abs(dB3) > 0.010 | abs(dS3) > 0.004)), nnz(v3));

% slopes with / without the dye mask on 1510, and eta from the masked slopes
obs1510 = getobs(1510);
s_off = st;  s_on = st;  s_on.dye_mask = true;
[Sx0, Sy0] = cisg_calculate_slopes(ref_b, obs1510, s_off);
[Sx1, Sy1] = cisg_calculate_slopes(ref_b, obs1510, s_on);
Sx0 = rg(Sx0);  Sx1 = rg(Sx1);  Sy1 = rg(Sy1);  Sy0 = rg(Sy0);
[eta1, ei] = cisg_slope_to_eta(Sx1, Sy1, 10*st.cm_per_pixel_x, 10*st.cm_per_pixel_y, 8);
[eta0, ~] = cisg_slope_to_eta(Sx0, Sy0, 10*st.cm_per_pixel_x, 10*st.cm_per_pixel_y, 8);
fprintf('eta (frame 1510): rms %.2f mm without dye mask, %.2f mm with (filled %.0f%% of coarse cells); rms of difference %.2f mm\n', ...
        rms(eta0(:)), rms(eta1(:)), 100*ei.nan_fraction, rms(eta1(:) - eta0(:)));

fig = figure('Visible', 'off', 'Position', [50 50 1500 900]);
subplot(2,2,1); plot(1:numel(dB1), 100*dB1, 'r', 1:numel(dS1), 100*dS1, 'b', 1:numel(dB3), 100*dB3, 'r:', 1:numel(dS3), 100*dS3, 'b:');
yline([1 -1 0.4 -0.4], 'k--'); xlabel('ROI row (binned)'); ylabel('row-median residual deviation (%)');
legend('B, 1510', 'brightness, 1510', 'B, 3000', 'brightness, 3000'); title('dye stripes: row-median B / brightness residual'); grid on;
q = 8;
subplot(2,2,2); imagesc(Sx0(1:q:end, 1:q:end)); axis image; caxis([-0.4 0.4]); colorbar; title('S_x frame 1510, no dye mask');
subplot(2,2,3); imagesc(Sx1(1:q:end, 1:q:end)); axis image; caxis([-0.4 0.4]); colorbar; title('S_x frame 1510, dye mask on');
subplot(2,2,4); imagesc(eta1); axis image; colorbar; title('\eta (mm) from the masked slopes (gaps filled)');
exportgraphics(fig, fullfile(here, 'dye_mask_check.png'), 'Resolution', 90); close(fig);

%% ---- 2. edge masks ----
fprintf('\nmasked pixels by reason (ROI %d-%d x %d-%d, bin %d; percent of ROI):\n', roi.row_min, roi.row_max, roi.col_min, roi.col_max, bin);
Sr = imgaussfilt(single(mean(ref_b, 3)) / 255, st.colour_sigma);
Cr = imgaussfilt(single(ref_b(:,:,1:2)) / 255, st.colour_sigma);
[ax, ay] = gradient(Cr(:,:,1));  [bx, by] = gradient(Cr(:,:,2));
dtm = ax .* by - ay .* bx;  f2 = ax.^2 + ay.^2 + bx.^2 + by.^2;
s1 = (f2 + sqrt(max(f2.^2 - 4 * dtm.^2, 0))) / 2;  cnd = s1 ./ max(abs(dtm), eps('single'));
tcb = (st.target_cols - 0.5);  [Hb, Wb] = size(Sr);
colsin = repmat((1:Wb) >= tcb(1) & (1:Wb) <= tcb(2), Hb, 1);
det_min = 0.1 * median(abs(dtm(colsin)));
ref_dark = rg(Sr < 0.15 | ~colsin);                         % dark or off-target in the reference itself
ref_cond = rg(abs(dtm) < det_min | cnd > 10) & ~ref_dark;   % reference colour map not invertible at the pixel
names = {'reference dark / outside target', 'reference ill-conditioned (at pixel)', 'x+delta out of bounds', 'displacement clamped', ...
         'ill-conditioned at x+delta', 'not converged', 'colour residual'};
figs = figure('Visible', 'off', 'Position', [50 50 1500 800]);
k = 0;
for n = [1510 3000]
    [~, ~, ~, D] = cisg_colour_displacement(ref_b, getobs(n), st);
    b = D.bits;  masked = b > 0 | isnan(D.dx);
    cat = zeros(size(b), 'uint8');                          % exclusive, in priority order
    cat(masked) = 8;
    flags = {ref_dark, ref_cond, bitand(b, 1) > 0, bitand(b, 2) > 0, bitand(b, 4) > 0, bitand(b, 8) > 0, bitand(b, 16) > 0};
    for j = numel(flags):-1:1, cat(masked & flags{j}) = j; end
    nr = size(b, 1);  nc = size(b, 2);
    bands = {'top 10% rows', 'bottom 10%', 'left 10% cols', 'right 10%', 'interior'};
    bm = false(nr, nc, 5);
    bm(1:round(.1*nr), :, 1) = true;  bm(round(.9*nr):end, :, 2) = true;
    bm(:, 1:round(.1*nc), 3) = true;  bm(:, round(.9*nc):end, 4) = true;
    bm(:, :, 5) = ~any(bm(:, :, 1:4), 3);
    fprintf('frame %d: total masked %.1f%%\n', n, 100*mean(masked(:)));
    for j = 1:7
        sel = cat == j;
        fprintf('  %-38s %5.1f%%   by band (%% of band): ', names{j}, 100*mean(sel(:)));
        for kb = 1:5, fprintf('%s %.0f  ', bands{kb}, 100*mean(sel(bm(:,:,kb)))); end, fprintf('\n');
    end
    other = cat == 8;  fprintf('  %-38s %5.1f%%\n', 'other (not explained by flags)', 100*mean(other(:)));
    k = k + 1;  subplot(1, 2, k);  imagesc(cat(1:8:end, 1:8:end)); axis image; colormap(gca, [1 1 1; lines(7)]); caxis([0 8]);
    title(sprintf('frame %d: 0 ok, 1 ref dark, 2 ref cond, 3 oob, 4 clamp, 5 ill-cond, 6 no conv, 7 residual', n));
end
exportgraphics(figs, fullfile(here, 'edge_mask_reasons.png'), 'Resolution', 90); close(figs);

% ROI from the reference
[ar, ~] = cisg_auto_roi(R.ref_image, R.setup);
[ar2, ~] = cisg_auto_roi(R.ref_image, R.setup, struct('margin_px', 100));
fprintf('\nreference-based ROI: rows %d-%d, cols %d-%d (margin 0); rows %d-%d, cols %d-%d (margin 100 px)\n', ...
        ar.row_min, ar.row_max, ar.col_min, ar.col_max, ar2.row_min, ar2.row_max, ar2.col_min, ar2.col_max);
for rr = {ar, ar2}
    [~, stR] = cisg_bin_setup(R.ref_image, setfield(R.setup, 'roi', rr{1}), bin);      %#ok<SFLD>
    for n = [1510 3000]
        [~, ~, iR] = cisg_colour_displacement(ref_b, getobs(n), stR);
        fprintf('ROI rows %d-%d cols %d-%d, frame %d: masked %.1f%%\n', rr{1}.row_min, rr{1}.row_max, rr{1}.col_min, rr{1}.col_max, n, 100*iR.mask_fraction);
    end
end

function s = rng(v)
    if isempty(v), s = ''; return; end
    brk = [0 find(diff(v(:)') > 1) numel(v)];  p = cell(1, numel(brk) - 1);
    for k = 1:numel(brk) - 1, p{k} = sprintf('%d-%d', v(brk(k) + 1), v(brk(k + 1))); end
    s = strjoin(p, ',');
end
function [dB, dS, valid] = rowdev(D, win)
    g = isfinite(D.dx) & D.bits == 0;
    eB = D.eB;  eB(~g) = NaN;  eS = D.eS;  eS(~g) = NaN;
    rB = median(eB, 2, 'omitnan');  rS = median(eS, 2, 'omitnan');
    valid = mean(g, 2) > 0.2;  rB(~valid) = NaN;  rS(~valid) = NaN;
    dB = rB - movmedian(rB, win, 'omitnan');  dS = rS - movmedian(rS, win, 'omitnan');
    dB(~valid) = 0;  dS(~valid) = 0;
end
function p = valid_px(D, rows)
    g = isfinite(D.dx) & D.bits == 0;  g(~rows, :) = false;
    idx = find(g);  idx = idx(1:max(1, round(numel(idx) / 20000)):end);
    p.eB = D.eB(idx);  p.eS = D.eS(idx);
end
function v = stats(x), v = [median(x) prctile(x, 10) prctile(x, 90)]; end
function a = auc(x, y)
    % P(x > y) by ranks
    n1 = numel(x);  n2 = numel(y);
    r = tiedrank([x(:); y(:)]);  a = (sum(r(1:n1)) - n1 * (n1 + 1) / 2) / (n1 * n2);
end
