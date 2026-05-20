%% MAKE_ALL_PLOTS_CSIG  Veron-style wave-field plotting suite from CSIG eta(x,y,t).
%
% Loads CISG_slopes_CoreView_22.mat (produced by cisg_main.m), reconstructs
% eta from (Sx, Sy), and renders every CSIG-derivable figure from
% D:\Scripps\OSBL-notes\plotting_plan.tex categories A, B, C, E.3.
%
% Stage windows are deferred — first cut runs on the full record. Hardcode
% stage onset times later and add overlays.

clear; clc;

addpath(fullfile(pwd, 'util'));
addpath(fullfile(pwd, 'plotting'));

mat_file = 'Y:\HLAB_2026\SlopeGauge\Processed\CISG_slopes_CoreView_22.mat';
fprintf('Loading %s ...\n', mat_file);
% S = load(mat_file);

% Sx    = S.Sx;
% Sy    = S.Sy;
% time  = S.time;
% setup = S.setup;
% clear S;
load(mat_file,'setup');
%%
dx = setup.cm_per_pixel_x / 100;     % m
dy = setup.cm_per_pixel_y / 100;     % m
dt = setup.dt;
% Zero outside ROI for FFT cleanliness
mask = setup.valid_mask;
Sx(~mask) = 0;
Sy(~mask) = 0;

%% Pre-eta: preview ROI on raw obs frame and let user tighten it
frame_num10x = 228;       % index into processed array
frame_num = frame_num10x*10+1;      % corresponding raw CoreView frame number
coreview_num  = 22;        % matches mat_file name
obs_folderName = sprintf('Y:\\HLAB_2026\\SlopeGauge\\CoreView_%d\\Flare 12M125 CCL C1112A00\\', coreview_num);
obs_fileName   = sprintf("CoreView_%d_Flare 12M125 CCL C1112A00_%04d.raw", ...
                          coreview_num, frame_num);
[obs_image, ~] = cisg_load_coreview(strcat(obs_folderName, obs_fileName));
%%
coreview_ref_num = 19;%coreview_num
ref_fileName= sprintf('Y:\\HLAB_2026\\SlopeGauge\\CoreView_%d\\Flare 12M125 CCL C1112A00\\CoreView_%d_Flare 12M125 CCL C1112A00_01.raw', coreview_ref_num,coreview_ref_num);
[ref_image, meta_ref] = cisg_load_coreview(ref_fileName);
%%
[ny, nx, ~] = size(obs_image);
x_cm = ((1:nx) - nx/2) * setup.cm_per_pixel_x;
y_cm = ((1:ny) - ny/2) * setup.cm_per_pixel_y;

figure('Name', sprintf('Frame raw %d - redraw ROI',frame_num));
imagesc(x_cm, y_cm, obs_image); axis image; set(gca,'YDir','normal');
hold on;

r = setup.roi;
x_box_cm = ([r.col_min r.col_max] - nx/2) * setup.cm_per_pixel_x;
y_box_cm = ([r.row_min r.row_max] - ny/2) * setup.cm_per_pixel_y;
rectangle('Position', [x_box_cm(1) y_box_cm(1) diff(x_box_cm) diff(y_box_cm)], ...
          'EdgeColor','r','LineWidth',1.5);
title('Red = current ROI. Click top-left then bottom-right of new ROI.');

fprintf('Click top-left then bottom-right of new ROI...\n');
[cx_cm, cy_cm] = ginput(2);
cx = cx_cm / setup.cm_per_pixel_x + nx/2;
cy = cy_cm / setup.cm_per_pixel_y + ny/2;

col_min = ceil(  min(cx) );
col_max = floor( max(cx) );
row_min = ceil(  min(cy) );
row_max = floor( max(cy) );

setup.roi = struct('row_min',row_min,'row_max',row_max, ...
                   'col_min',col_min,'col_max',col_max);
setup.valid_mask = false(ny, nx);
setup.valid_mask(row_min:row_max, col_min:col_max) = true;

mask = setup.valid_mask;
Sx(~mask) = 0;
Sy(~mask) = 0;

fprintf('New ROI: rows [%d %d], cols [%d %d]\n', row_min, row_max, col_min, col_max);

roi_out_file = 'Y:\HLAB_2026\SlopeGauge\Processed\CISG_ROI_CoreView_22.mat';
roi        = setup.roi;          %#ok<NASGU>
valid_mask = setup.valid_mask;   %#ok<NASGU>
%%
figure; 
subplot(1,2,1)
imagesc(x_cm(col_min:col_max), y_cm(row_min:row_max), ref_image(row_min:row_max,col_min:col_max,:)); axis equal tight; set(gca,'YDir','normal');
xlabel('$x$ (cm)', 'interpreter','latex');
ylabel('$y$ (cm)', 'interpreter','latex');
set(gca, 'fontsize',14,'fontname','times')
title('Reference Image','Interpreter','latex');

subplot(1,2,2)
imagesc(x_cm(col_min:col_max), y_cm(row_min:row_max), obs_image(row_min:row_max,col_min:col_max,:)); axis equal tight; set(gca,'YDir','normal');
xlabel('$x$ (cm)', 'interpreter','latex');
ylabel('$y$ (cm)', 'interpreter','latex');
set(gca, 'fontsize',14,'fontname','times')
title(sprintf('Frame at $t = %.2f$ s', (frame_num-1)*dt*10),'Interpreter','latex');
%%
save(roi_out_file, 'roi', 'valid_mask', 'frame_to_show', 'raw_frame_num');
fprintf('Saved tightened ROI to %s\n', roi_out_file);

%%
[sx, sy] = cisg_calculate_slopes(ref_image, obs_image, setup);

%%
figure('Position', [100 100 1000 600], 'Color', 'white');

% Sx plot
subplot(1,2,1)
imagesc(x_cm(col_min:col_max), y_cm(row_min:row_max), sx(row_min:row_max,col_min:col_max));
axis equal tight;
c=colorbar;
colormap(gca, brewermap([],'Spectral'));
caxis([-0.4, 0.4]);
title(sprintf('$S_x$ at $t = %.2f$ s', (frame_num-1)*dt*10),'Interpreter','latex');
xlabel('$x$ (cm)', 'interpreter','latex');
ylabel('$y$ (cm)', 'interpreter','latex');
set(gca, 'fontsize',14,'fontname','times')
c.Label.String='$ak$';
c.Label.Interpreter='latex';
c.Label.FontSize=16;
% Sy plot
subplot(1,2,2)
imagesc(x_cm(col_min:col_max), y_cm(row_min:row_max), sy(row_min:row_max,col_min:col_max));
axis equal tight;
c=colorbar;
colormap(gca, brewermap([],'Spectral'));
caxis([-0.2, 0.2]);
title(sprintf('$S_y$ at $t = %.2f$ s', (frame_num-1)*dt*10),'Interpreter','latex');
xlabel('$x$ (cm)', 'interpreter','latex');
ylabel('$y$ (cm)', 'interpreter','latex');
set(gca, 'fontsize',14,'fontname','times')
c.Label.String='$ak$';
c.Label.Interpreter='latex';
c.Label.FontSize=16;
%%
fprintf('Reconstructing eta from slopes ...\n');
eta = cisg_reconstruct_eta(sx, sy, dx, dy);

if ~exist('figures', 'dir'), mkdir figures; end

data = struct('Sx', Sx, 'Sy', Sy, 'eta', eta, 'time', time, ...
              'setup', setup, 'dx', dx, 'dy', dy);

fprintf('\n=== Generating Veron-style CSIG plots ===\n');
%% plot_A1_snapshots_csig(data);

%%

plot_A2_spectrograms_csig(data);
plot_A3_kspectra_csig(data);
plot_A4_slope_variance(data);
plot_A6_directional_B(data);
plot_B1_kx_omega_csig(data);
plot_B3_phase_speed_dev(data);
plot_B4_angular_dist(data);
plot_C1_local_k_map(data);
plot_C4_cross_stream_coh(data);
plot_C5_heterogeneity_idx(data);
plot_C6_local_stokes_map(data);
plot_E3_bandwidth_correction(data);

fprintf('\nAll CSIG plots written to figures/\n');
