%% MAKE_ALL_PLOTS_CSIG  Veron-style wave-field plotting suite from CSIG eta(x,y,t).
%
% Reads the per-frame slope cache built by cisg_main.m
% (D:\HLAB_2026\SlopeGauge\CISG_slopes_CoreView_NN\frame_NNNN.mat),
% reconstructs eta from (Sx, Sy), and renders the CSIG plot suite.

clear; clc;

addpath(fullfile(pwd, 'util'));
addpath(fullfile(pwd, 'plotting'));

coreview_num = 22;
cache_root   = 'D:\HLAB_2026\SlopeGauge';
cache_dir    = fullfile(cache_root, sprintf('CISG_slopes_CoreView_%d', coreview_num));

%% Load cache metadata
meta              = load(fullfile(cache_dir, 'metadata.mat'));
setup             = meta.setup;
time              = meta.time;

dx = setup.cm_per_pixel_x / 100;     % m
dy = setup.cm_per_pixel_y / 100;     % m
dt = setup.dt;

fprintf('Cache: %s\n', cache_dir);
fprintf('Outer ROI baked in cache: rows [%d %d] cols [%d %d]\n', ...
        setup.roi.row_min, setup.roi.row_max, setup.roi.col_min, setup.roi.col_max);

%% Tighten the ROI for plotting (coords are within the cached frame)
% The saved CISG_ROI file holds an inner ROI in ORIGINAL full-image coords.
% Cached frames are already cropped to the outer ROI, so translate by the
% outer-ROI offset and clamp to the cached-frame bounds.
outer_roi = setup.roi;   % outer ROI from cache metadata (full-image coords)
roi_file  = '\\Airseaserver28\D\HLAB_2026\SlopeGauge\Processed\CISG_ROI_CoreView_22.mat';
saved     = load(roi_file, 'roi');

inner_roi = struct( ...
    'row_min', max(1,         saved.roi.row_min - outer_roi.row_min + 1), ...
    'row_max', min(cached_ny, saved.roi.row_max - outer_roi.row_min + 1), ...
    'col_min', max(1,         saved.roi.col_min - outer_roi.col_min + 1), ...
    'col_max', min(cached_nx, saved.roi.col_max - outer_roi.col_min + 1));

if inner_roi.row_min >= inner_roi.row_max || inner_roi.col_min >= inner_roi.col_max
    error('make_all_plots_csig:roi_outside_cache', ...
          'Saved inner ROI does not overlap the cached outer ROI.');
end

fprintf('Saved inner ROI (full coords): rows [%d %d] cols [%d %d]\n', ...
        saved.roi.row_min, saved.roi.row_max, saved.roi.col_min, saved.roi.col_max);
fprintf('Translated to cached coords:   rows [%d %d] cols [%d %d]\n', ...
        inner_roi.row_min, inner_roi.row_max, inner_roi.col_min, inner_roi.col_max);

% inner_roi = pick_roi_interactive(Sx_one, [], ...
%     setup.cm_per_pixel_x, setup.cm_per_pixel_y, ...
%     sprintf('Inner ROI within cached frame %d', frame_to_show));

setup.roi        = inner_roi;
setup.valid_mask = false(cached_ny, cached_nx);
setup.valid_mask(inner_roi.row_min:inner_roi.row_max, ...
                 inner_roi.col_min:inner_roi.col_max) = true;

%% Load one cached frame for inner-ROI preview
frame_to_show = 2281;     % raw frame number; must exist in the cache
[Sx_one, Sy_one] = load_slope_frame(cache_dir, frame_to_show);
[cached_ny, cached_nx]  = size(Sx_one);

t_one = time(frame_to_show);
x_cm = ((1:cached_nx) - cached_nx/2) * setup.cm_per_pixel_x;
y_cm = ((1:cached_ny) - cached_ny/2) * setup.cm_per_pixel_y;

% Single-frame Sx, Sy preview (from cache)
figure('Position', [100 100 1000 600], 'Color', 'white');

subplot(1,2,1);
imagesc(x_cm(inner_roi.col_min:inner_roi.col_max), ...
        y_cm(inner_roi.row_min:inner_roi.row_max), ...
        Sx_one(inner_roi.row_min:inner_roi.row_max, inner_roi.col_min:inner_roi.col_max));
axis equal tight; set(gca, 'YDir', 'normal');
c = colorbar; colormap(gca, brewermap([], 'Spectral')); caxis([-0.4, 0.4]);
title(sprintf('$S_x$ at $t = %.2f$ s', t_one), 'Interpreter', 'latex');
xlabel('$x$ (cm)', 'interpreter', 'latex');
ylabel('$y$ (cm)', 'interpreter', 'latex');
set(gca, 'fontsize', 14, 'fontname', 'times');
c.Label.String = '$ak$'; c.Label.Interpreter = 'latex'; c.Label.FontSize = 16;

subplot(1,2,2);
imagesc(x_cm(inner_roi.col_min:inner_roi.col_max), ...
        y_cm(inner_roi.row_min:inner_roi.row_max), ...
        Sy_one(inner_roi.row_min:inner_roi.row_max, inner_roi.col_min:inner_roi.col_max));
axis equal tight; set(gca, 'YDir', 'normal');
c = colorbar; colormap(gca, brewermap([], 'Spectral')); caxis([-0.2, 0.2]);
title(sprintf('$S_y$ at $t = %.2f$ s', t_one), 'Interpreter', 'latex');
xlabel('$x$ (cm)', 'interpreter', 'latex');
ylabel('$y$ (cm)', 'interpreter', 'latex');
set(gca, 'fontsize', 14, 'fontname', 'times');
c.Label.String = '$ak$'; c.Label.Interpreter = 'latex'; c.Label.FontSize = 16;

%% Single-frame Sx, Sy preview from RAW (sanity check vs cached)
coreview_ref_num = 19;     % must match cisg_main
raw_root         = '\\Airseaserver28\D\HLAB_2026\SlopeGauge';
ref_fileName     = fullfile(raw_root, sprintf('CoreView_%d\\Flare 12M125 CCL C1112A00\\CoreView_%d_Flare 12M125 CCL C1112A00_01.raw', ...
                                              coreview_ref_num, coreview_ref_num));
obs_folderName   = fullfile(raw_root, sprintf('CoreView_%d\\Flare 12M125 CCL C1112A00\\', coreview_num));
obs_fileName     = sprintf('CoreView_%d_Flare 12M125 CCL C1112A00_%04d.raw', coreview_num, frame_to_show);

fprintf('Loading raw ref + obs frame %d for sanity check ...\n', frame_to_show);
[ref_image, ~] = cisg_load_coreview(ref_fileName);
[obs_image, ~] = cisg_load_coreview(strcat(obs_folderName, obs_fileName));

[sx_raw, sy_raw] = cisg_calculate_slopes(ref_image, obs_image, setup);

% Crop to the outer ROI (matches what the cache stores), display in cm.
r = outer_roi;
[ny_full, nx_full, ~] = size(ref_image);
x_cm_full = ((1:nx_full) - nx_full/2) * setup.cm_per_pixel_x;
y_cm_full = ((1:ny_full) - ny_full/2) * setup.cm_per_pixel_y;

figure('Position', [100 100 1000 600], 'Color', 'white', 'Name', 'Raw-derived S_x, S_y');

subplot(1,2,1);
imagesc(x_cm_full(r.col_min:r.col_max), y_cm_full(r.row_min:r.row_max), ...
        sx_raw(r.row_min:r.row_max, r.col_min:r.col_max));
axis equal tight; set(gca, 'YDir', 'normal');
c = colorbar; colormap(gca, brewermap([], 'Spectral')); caxis([-0.4, 0.4]);
title(sprintf('$S_x$ (raw) at $t = %.2f$ s', (frame_to_show-1)*dt), 'Interpreter', 'latex');
xlabel('$x$ (cm)', 'interpreter', 'latex');
ylabel('$y$ (cm)', 'interpreter', 'latex');
set(gca, 'fontsize', 14, 'fontname', 'times');
c.Label.String = '$ak$'; c.Label.Interpreter = 'latex'; c.Label.FontSize = 16;

subplot(1,2,2);
imagesc(x_cm_full(r.col_min:r.col_max), y_cm_full(r.row_min:r.row_max), ...
        sy_raw(r.row_min:r.row_max, r.col_min:r.col_max));
axis equal tight; set(gca, 'YDir', 'normal');
c = colorbar; colormap(gca, brewermap([], 'Spectral')); caxis([-0.2, 0.2]);
title(sprintf('$S_y$ (raw) at $t = %.2f$ s', (frame_to_show-1)*dt), 'Interpreter', 'latex');
xlabel('$x$ (cm)', 'interpreter', 'latex');
ylabel('$y$ (cm)', 'interpreter', 'latex');
set(gca, 'fontsize', 14, 'fontname', 'times');
c.Label.String = '$ak$'; c.Label.Interpreter = 'latex'; c.Label.FontSize = 16;

%% Single-frame eta preview
fprintf('Reconstructing eta for preview frame %d ...\n', frame_to_show);
eta_one = cisg_reconstruct_eta( ...
    Sx_one(inner_roi.row_min:inner_roi.row_max, inner_roi.col_min:inner_roi.col_max), ...
    Sy_one(inner_roi.row_min:inner_roi.row_max, inner_roi.col_min:inner_roi.col_max), ...
    dx, dy);

figure;
imagesc(x_cm(inner_roi.col_min:inner_roi.col_max), ...
        y_cm(inner_roi.row_min:inner_roi.row_max), eta_one * 100);
axis equal tight; set(gca, 'YDir', 'normal');
caxis([-0.4,0.4])
c = colorbar; colormap(gca, brewermap([], 'Spectral'));
title(sprintf('$\\eta$ at $t = %.2f$ s', t_one), 'Interpreter', 'latex');
xlabel('$x$ (cm)', 'interpreter', 'latex');
ylabel('$y$ (cm)', 'interpreter', 'latex');
set(gca, 'fontsize', 14, 'fontname', 'times');
c.Label.String = '$\eta$ (cm)'; c.Label.Interpreter = 'latex'; c.Label.FontSize = 16;
set(gcf, 'Color', 'white');

%% Materialize the inner-ROI cube and run the plot suite
% load_slope_cube crops to inner_roi as it loads, so allocations are
% (inner_ny * inner_nx * Nt * 4 bytes) * 3 arrays (Sx, Sy, eta).
frame_subset = 1111:10:length(time);     % e.g. frames_to_process(1:500) for a quick test
fprintf('Loading slope cube (%d frames, inner-ROI cropped) ...\n', numel(frame_subset));
[Sx, Sy] = load_slope_cube(cache_dir, frame_subset, inner_roi);
time_subset = (frame_subset - 1) * dt;

%% Rest-period per-pixel offset subtraction
% Removes time-invariant artifacts (ref-cal error, stuck spots present at rest).
% rest_frames must be cached and use the same outer ROI (same cached layout).
rest_frames = 1:50;
fprintf('Loading rest cube (%d frames) for per-pixel offset ...\n', numel(rest_frames));
[Sx_rest, Sy_rest] = load_slope_cube(cache_dir, rest_frames, inner_roi);
Sx_offset = mean(Sx_rest, 3);
Sy_offset = mean(Sy_rest, 3);
clear Sx_rest Sy_rest;
fprintf('Subtracting per-pixel rest offsets.\n');
Sx = Sx - Sx_offset;
Sy = Sy - Sy_offset;

fprintf('Reconstructing eta cube ...\n');
eta = cisg_reconstruct_eta(Sx, Sy, dx, dy);

if ~exist('figures', 'dir'), mkdir figures; end

% Cube is already inner-ROI sized — plot functions read setup.roi to slice,
% so reset to the trivial full-extent ROI on the cropped layout.
[inner_ny, inner_nx, ~] = size(Sx);
setup.roi        = struct('row_min', 1, 'row_max', inner_ny, ...
                          'col_min', 1, 'col_max', inner_nx);
setup.valid_mask = true(inner_ny, inner_nx);

data = struct('Sx', Sx, 'Sy', Sy, 'eta', eta, 'time', time_subset, ...
              'setup', setup, 'dx', dx, 'dy', dy);

%% Compute (ak)^2 time series
ak2_Sx = squeeze(mean(Sx.^2, [1 2]));
ak2_Sy = squeeze(mean(Sy.^2, [1 2]));

figure('Name', 'Mean Square Slope Evolution');
plot(time_subset, ak2_Sx, 'LineWidth', 1.5); hold on;
plot(time_subset, ak2_Sy, 'LineWidth', 1.5);
xlabel('Time (s)');
ylabel('$(ak)^2$', 'Interpreter', 'latex');
legend('$S_x$', '$S_y$', 'Interpreter', 'latex','location','northwest');
grid off;
set(gca, 'fontsize', 14, 'fontname', 'times');
set(gcf, 'Color', 'white');
%%
% Optional 2nd arg: time indices into data.eta (1..numel(time)). Omit for 5 evenly spaced.
target_times = [25,30,40,49];     % seconds
[~, idx_set] = min(abs(time_subset(:) - target_times(:).'), [], 1);
plot_A1_snapshots_csig(data, idx_set);
% plot_A1_snapshots_csig(data, [1 25 60 100]);             % raw index version

%%
% k_x-vs-t intensity of S_x, S_y (per-row FFT in x, averaged over y).
% Omit 2nd arg to use all snapshots, or pass idx_set for the chosen times.
plot_kx_t_slopes_csig(data);
%%
plot_f_t_slopes_csig(data);
%%
plot_kx_omega_Sx_csig(data);
%%
% plot_A2_spectrograms_csig(data);
% %%
% plot_A3_kspectra_csig(data);
% %%
% plot_A4_slope_variance(data);
% %%
% plot_A6_directional_B(data);
% plot_B1_kx_omega_csig(data);
% 
% plot_B3_phase_speed_dev(data);
% plot_B4_angular_dist(data);
% plot_C1_local_k_map(data);
% plot_C4_cross_stream_coh(data);
% plot_C5_heterogeneity_idx(data);
% plot_C6_local_stokes_map(data);
% plot_E3_bandwidth_correction(data);
% 
% fprintf('\nAll CSIG plots written to figures/\n');
