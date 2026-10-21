%% cisg_frame3000_check.m
% Linear vs local colour mapping (with / without lighting-drift correction and
% mean-slope removal) on a wavy frame.  Writes calibration/cisg_frame3000_check.png.
clear; clc;
here = fileparts(mfilename('fullpath'));
sg = fileparts(here);
addpath(here, fullfile(sg, 'util'), sg, fullfile(fileparts(sg), 'IR Camera', 'calibration'));

coreview_num = 95;
frame = 3000;
folder = fullfile('D:\HLAB_2026\Core', sprintf('CoreView_%d', coreview_num), 'Flare 12M125 CCL C1112A00');
fname = @(n) fullfile(folder, sprintf('CoreView_%d_Flare 12M125 CCL C1112A00_%04d.raw', coreview_num, n));

S = load(fullfile(here, 'cisg_cal_CoreView_97.mat'));  cal = S.cal;
ref = cisg_load_coreview(fname(1));
obs = cisg_load_coreview(fname(frame));
[H, W, ~] = size(ref);
depth = cisg_apparent_depth(ref, cal);
setup = struct('n_water', 1.333, 'n_air', 1, 'use_apparent_depth', true, 'water_depth', 48, ...
               'colour_map', 'local', 'colour_norm', 'none', 'colour_sigma', 12, 'max_disp_px', 600, ...
               'L_eff_cm', depth.L_eff_mm / 10, 'pattern_scale', depth.r, 'target_cols', depth.edge_cols);
setup.roi = struct('row_min', round(0.2*H), 'row_max', round(0.8*H), 'col_min', round(0.2*W), 'col_max', round(0.8*W));
r = setup.roi;
rc = [mean([r.col_min r.col_max]) mean([r.row_min r.row_max])];
[setup.cm_per_pixel_x, setup.cm_per_pixel_y] = cisg_local_scale(cal, rc);
rg = @(A) A(r.row_min:r.row_max, r.col_min:r.col_max, :);

d = double(rg(obs)) - double(rg(ref));
fprintf('frame %d - frame 1 colour difference (ROI): rms R %.2f G %.2f B %.2f counts; mean %+.2f %+.2f %+.2f\n', frame, ...
        rms(d(:,:,1), 'all'), rms(d(:,:,2), 'all'), rms(d(:,:,3), 'all'), ...
        mean(d(:,:,1), 'all'), mean(d(:,:,2), 'all'), mean(d(:,:,3), 'all'));

%% ---- four ways ----
names = {'linear', 'local, drift none', 'local, drift gain', 'local, gain, mean removed'};
Sx = cell(1, 4);  Sy = cell(1, 4);
s = setup;  s.colour_map = 'linear';                      [Sx{1}, Sy{1}] = cisg_calculate_slopes(ref, obs, s);
s = setup;  s.drift_corr = 'none';                        [Sx{2}, Sy{2}] = cisg_calculate_slopes(ref, obs, s);
s = setup;  s.drift_corr = 'gain';                        [Sx{3}, Sy{3}] = cisg_calculate_slopes(ref, obs, s);
[~, ~, info] = cisg_colour_displacement(ref, obs, setfield(s, 'roi', r));
Sx{4} = Sx{3} - mean(rg(Sx{3}), 'all', 'omitnan');  Sy{4} = Sy{3} - mean(rg(Sy{3}), 'all', 'omitnan');
fprintf('drift fit [a b] R: %.4f %+.4f ; G: %.4f %+.4f  (colour units 0..1; b*255 = %+.2f, %+.2f counts)\n', ...
        info.drift_ab(1,:), info.drift_ab(2,:), 255*info.drift_ab(1,2), 255*info.drift_ab(2,2));
for i = 1:4
    a = rg(Sx{i});  b = rg(Sy{i});  ok = isfinite(a) & isfinite(b);
    la = rg(Sx{1});  lb = rg(Sy{1});  ok1 = ok & isfinite(la);
    fprintf('%-26s rms Sx %.4f Sy %.4f | mean Sx %+.4f Sy %+.4f | masked %.1f%% | corr with linear Sx %.3f Sy %.3f\n', ...
            names{i}, rms(a(ok)), rms(b(ok)), mean(a(ok)), mean(b(ok)), 100*mean(~ok(:)), ...
            corr(a(ok1), la(ok1)), corr(b(ok1), lb(ok1)));
end

%% ---- wavelength along x (tank x = image columns) from the Sx spectrum ----
mmpx = 10 * setup.cm_per_pixel_x;
for i = [1 4]
    A = rg(Sx{i});  A = A - mean(A, 'all', 'omitnan');  A(isnan(A)) = 0;
    rows = 1:16:size(A, 1);
    n = size(A, 2);  win = hann(n)';
    P = zeros(1, floor(n/2) + 1);
    for k = rows
        F = abs(fft((A(k,:) - mean(A(k,:))) .* win)).^2;
        P = P + F(1:floor(n/2) + 1);
    end
    f = (0:floor(n/2)) / n;                                   % cycles/px
    ok = f > 1.5 / n;                                         % >= 1.5 cycles across the ROI
    [~, ip] = max(P(ok));  fo = f(ok);
    pk = fo(ip);
    Pw = P(ok);  fm = sum(fo .* Pw) / sum(Pw);                % spectral centroid
    fprintf('%-26s Sx peak wavelength %.0f mm (centroid %.0f mm); ROI %.0f mm wide\n', names{i}, mmpx / pk, mmpx / fm, n * mmpx);
end

%% ---- figure ----
q = 8;  d8 = @(A) A(1:q:end, 1:q:end);
cols = r.col_min:q:r.col_max;  rows = r.row_min:q:r.row_max;
xv = cisg_pixel_grid(cal, [round(rc(2)) round(rc(2)) cols(1) cols(end)]);           % x (cm) along the mid row
[~, yv] = cisg_pixel_grid(cal, [rows(1) rows(end) round(rc(1)) round(rc(1))]);      % y (cm) along the mid column
xv = double(xv(:)');  yv = double(yv(:));
cl = 0.15;
fig = figure('Visible', 'off', 'Position', [100 100 1500 1100]);
tiledlayout(3, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
panels = {Sx{1}, 'S_x linear'; Sx{3}, 'S_x local + gain'; Sy{1}, 'S_y linear'; Sy{3}, 'S_y local + gain'};
for p = 1:4
    nexttile;
    imagesc(xv, yv, d8(rg(panels{p, 1})));  set(gca, 'YDir', 'normal');  axis image;
    caxis([-cl cl]);  colorbar;  title(panels{p, 2});  xlabel('x (cm)');  ylabel('y (cm)');
end
Sxl = rg(Sx{1});  Sxn = rg(Sx{3});  Syl = rg(Sy{1});  Syn = rg(Sy{3});
mr = round(size(Sxl, 1) / 2);
xf = linspace(xv(1), xv(end), size(Sxl, 2));
nexttile;  plot(xf, Sxl(mr, :), 'k', xf, Sxn(mr, :), 'r');  legend('linear', 'local + gain');
title('S_x, middle row');  xlabel('x (cm)');  grid on;
nexttile;  plot(xf, Syl(mr, :), 'k', xf, Syn(mr, :), 'r');  legend('linear', 'local + gain');
title('S_y, middle row');  xlabel('x (cm)');  grid on;
exportgraphics(fig, fullfile(here, 'cisg_frame3000_check.png'), 'Resolution', 100);
close(fig);
fprintf('Saved calibration/cisg_frame3000_check.png\n');
