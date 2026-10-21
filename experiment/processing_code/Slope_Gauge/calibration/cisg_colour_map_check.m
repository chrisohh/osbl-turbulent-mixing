%% cisg_colour_map_check.m
% Validation of the local colour -> position inversion (util/cisg_colour_displacement.m):
%   a. reference vs itself        -> delta ~ 0
%   b. synthetic known warp        -> recovered vs true displacement (with and without noise)
%   c. real frame vs the linear path (old global gradient fit)
% Writes calibration/cisg_colour_map_check.png.
clear; clc;
here = fileparts(mfilename('fullpath'));
sg = fileparts(here);
addpath(here, fullfile(sg, 'util'), sg, fullfile(fileparts(sg), 'IR Camera', 'calibration'));

%% ---- config ----
coreview_num = 95;
frame = 1101;
folder = fullfile('D:\HLAB_2026\Core', sprintf('CoreView_%d', coreview_num), 'Flare 12M125 CCL C1112A00');
fname = @(n) fullfile(folder, sprintf('CoreView_%d_Flare 12M125 CCL C1112A00_%04d.raw', coreview_num, n));

S = load(fullfile(here, 'cisg_cal_CoreView_97.mat'));  cal = S.cal;
ref = cisg_load_coreview(fname(1));
obs = cisg_load_coreview(fname(frame));
[H, W, ~] = size(ref);

depth = cisg_apparent_depth(ref, cal);
setup = struct('n_water', 1.333, 'n_air', 1, 'use_apparent_depth', true, 'water_depth', 48, ...
               'colour_map', 'local', 'colour_norm', 'sum', 'colour_sigma', 6, 'max_disp_px', 200, ...
               'L_eff_cm', depth.L_eff_mm / 10, 'pattern_scale', depth.r, 'target_cols', depth.edge_cols);
setup.roi = struct('row_min', round(0.2*H), 'row_max', round(0.8*H), 'col_min', round(0.2*W), 'col_max', round(0.8*W));
[setup.cm_per_pixel_x, setup.cm_per_pixel_y] = cisg_local_scale(cal, [mean([setup.roi.col_min setup.roi.col_max]) mean([setup.roi.row_min setup.roi.row_max])]);
rg = @(A) A(setup.roi.row_min:setup.roi.row_max, setup.roi.col_min:setup.roi.col_max);

%% ---- a. reference vs itself ----
[dx, dy, info] = cisg_colour_displacement(ref, ref, setup);
fprintf('a. ref vs itself: rms dx %.4f px, dy %.4f px; mask %.1f%%; %.1f s\n', ...
        rms(rg(dx), 'all', 'omitnan'), rms(rg(dy), 'all', 'omitnan'), 100*info.mask_fraction, info.seconds);

%% ---- b. synthetic warp ----
r = setup.roi;
[X, Y] = meshgrid(single(r.col_min - 40 : r.col_max + 40), single(r.row_min - 40 : r.row_max + 40));
tx = 8 * sin(2*pi*X / 600);  ty = 5 * cos(2*pi*Y / 500);          % true delta: obs(x) = ref(x + delta)
syn = ref;
for k = 1:3
    syn(r.row_min-40:r.row_max+40, r.col_min-40:r.col_max+40, k) = ...
        uint8(round(interp2(double(ref(:,:,k)), X + tx, Y + ty, 'linear')));
end
tx = tx(41:end-40, 41:end-40);  ty = ty(41:end-40, 41:end-40);
[sx, sy, info_s] = cisg_colour_displacement(ref, syn, setup);
ex = rg(sx) - tx;  ey = rg(sy) - ty;
fprintf('b. synthetic (clean): rms err x %.3f y %.3f px; bias x %+.3f y %+.3f px; mask %.1f%%; %.1f s\n', ...
        rms(ex, 'all', 'omitnan'), rms(ey, 'all', 'omitnan'), mean(ex, 'all', 'omitnan'), mean(ey, 'all', 'omitnan'), ...
        100*mean(isnan(ex), 'all'), info_s.seconds);
rng(1);
syn_n = uint8(double(syn) + 1.5 * randn(size(syn)));
[nx_, ny_] = cisg_colour_displacement(ref, syn_n, setup);
fprintf('   synthetic (+1.5 count noise): rms err x %.3f y %.3f px; bias x %+.3f y %+.3f; mask %.1f%%\n', ...
        rms(rg(nx_) - tx, 'all', 'omitnan'), rms(rg(ny_) - ty, 'all', 'omitnan'), ...
        mean(rg(nx_) - tx, 'all', 'omitnan'), mean(rg(ny_) - ty, 'all', 'omitnan'), 100*mean(isnan(rg(nx_) - tx), 'all'));

%% ---- c. real frame: local vs linear ----
setup_lin = setup;  setup_lin.colour_map = 'linear';
[Sx_l, Sy_l] = cisg_calculate_slopes(ref, obs, setup_lin);
tic; [Sx_n, Sy_n] = cisg_calculate_slopes(ref, obs, setup);  t_loc = toc;
[dxo, dyo, info_o] = cisg_colour_displacement(ref, obs, setup);
fprintf('c. frame %d: local rms Sx %.4f Sy %.4f (mask %.1f%%, %.1f s/frame) | linear rms Sx %.4f Sy %.4f\n', ...
        frame, rms(rg(Sx_n), 'all', 'omitnan'), rms(rg(Sy_n), 'all', 'omitnan'), 100*info_o.mask_fraction, t_loc, ...
        rms(rg(Sx_l), 'all', 'omitnan'), rms(rg(Sy_l), 'all', 'omitnan'));
ok = isfinite(rg(Sx_n));
a1 = rg(Sx_n);  a2 = rg(Sx_l);  b1 = rg(Sy_n);  b2 = rg(Sy_l);
fprintf('   same pixels (valid local): rms Sx local %.4f linear %.4f; Sy local %.4f linear %.4f; corr Sx %.3f Sy %.3f\n', ...
        rms(a1(ok)), rms(a2(ok)), rms(b1(ok)), rms(b2(ok)), corr(a1(ok), a2(ok)), corr(b1(ok), b2(ok)));
setup_n = setup;  setup_n.colour_norm = 'none';
[Sx_0, Sy_0] = cisg_calculate_slopes(ref, obs, setup_n);
fprintf('   colour_norm none: rms Sx %.4f Sy %.4f; mask %.1f%%\n', rms(rg(Sx_0), 'all', 'omitnan'), ...
        rms(rg(Sy_0), 'all', 'omitnan'), 100*mean(isnan(rg(Sx_0)), 'all'));

%% ---- figure ----
q = 4;  d4 = @(A) A(1:q:end, 1:q:end, :);
Cs = imgaussfilt(single(d4(ref)), 6/q);
Cn = cat(3, Cs(:,:,1) ./ max(sum(Cs, 3), 1), Cs(:,:,2) ./ max(sum(Cs, 3), 1));
[a, b] = gradient(Cn(:,:,1));  [c, d] = gradient(Cn(:,:,2));
dj = abs(a .* d - b .* c);          % per coarse px^2
fig = figure('Visible', 'off', 'Position', [100 100 1700 900]);
subplot(2,3,1); imagesc(dj); axis image; colorbar; title('|det J| of the reference (coarse px^{-2}), ROI in red');
hold on; rectangle('Position', [setup.roi.col_min setup.roi.row_min diff([setup.roi.col_min setup.roi.col_max]) diff([setup.roi.row_min setup.roi.row_max])] / q, 'EdgeColor', 'r');
subplot(2,3,2); imagesc(d4(isnan(dxo))); axis image; colormap(gca, gray); title(sprintf('mask (white = NaN), frame %d, %.1f%% in ROI', frame, 100*info_o.mask_fraction));
subplot(2,3,3); imagesc(ex); axis image; colorbar; title('synthetic: \deltax error (px)');
subplot(2,3,4); imagesc(d4(dxo)); axis image; colorbar; title(sprintf('recovered \deltax (px), frame %d', frame));
subplot(2,3,5); imagesc(d4(dyo)); axis image; colorbar; title(sprintf('recovered \deltay (px), frame %d', frame));
subplot(2,3,6); imagesc(ey); axis image; colorbar; title('synthetic: \deltay error (px)');
exportgraphics(fig, fullfile(here, 'cisg_colour_map_check.png'), 'Resolution', 100);
close(fig);
fprintf('Saved calibration/cisg_colour_map_check.png\n');
