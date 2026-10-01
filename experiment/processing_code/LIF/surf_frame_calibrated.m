%% Water_SURF frame with physical scale + dye front, timed against wind start
% Longitudinal (side) view, CoreView 95.  The image is shown AS RECORDED
% (no rectification); the checkerboard plate in CoreView 94 (same camera,
% same mount, 2.5 cm toward the camera) only sets the mm/px scale -- see
% util/lif_plate_scale.m.  Dye front from util/lif_dye_front.m.
%
% Time: frame n is trigger n of the 50 Hz camera counter, so
%     t = (camStartElapsed - fanStartElapsed) + (n - 1) / fs
% seconds since wind start, from the run_experiment.m hot-wire file.
clear; close all;
addpath(fullfile(fileparts(mfilename('fullpath')), 'util'));

%% ---------------- User settings ----------------
raw_root  = 'D:\HLAB_2026\Core';
cam       = 'Water_Surf';
run_num   = 95;     frame     = 1875;
plate_num = 94;     plate_frm = 1;          % plate files are numbered _01

% Plate scale (see lif_plate_scale): printed SCALED DOWN from the design.
sc.plate_scale = 1.0;      % <-- printed / designed size
sc.square_mm   = 11;       % small plate: coarse 11x10 layout with 11 mm squares
                           % (big plate = same layout at the design 40 mm); overrides plate_scale
sc.offset_mm   = 25;       % plate was this far TOWARD the camera
sc.n_medium    = 1.33;     % plate in water, viewed through the wall
sc.lens_f_mm   = 35;       % <-- lens focal length (CHECK); only sets D

hw_file = 'D:\HLAB_2026\hotwire\hotwire_20260923_113939.mat';
fs      = 50;

downsample = 2;            % display every Nth pixel
clims      = [350 600];    % counts; [] = percentiles below
pct_lims   = [1 99.5];
cmap       = 'bone';

% Dye front (see lif_dye_front)
dy.dye_frac       = 0.90;  % 0.9 ~ the hand-drawn line; 0.99 = outermost edge
dy.dye_thr        = [];    % counts; set e.g. 530 for a fixed threshold instead
dy.smooth_mm      = 1.5;
dy.surf_window_cm = [4 10];
dy.surf_gap_cm    = 0.4;

%% ---------------- Scale ----------------
fplate = fullfile(raw_root, sprintf('CoreView_%d', plate_num), cam, ...
                  sprintf('CoreView_%d_%s_%02d.raw', plate_num, cam, plate_frm));
[mmpp, c0] = lif_plate_scale(fplate, sc);

%% ---------------- Load frame ----------------
fname = fullfile(raw_root, sprintf('CoreView_%d', run_num), cam, ...
                 sprintf('CoreView_%d_%s_%04d.raw', run_num, cam, frame));
I = lif_load_raw(fname);
I = I(1:downsample:end, 1:downsample:end);
[ny, nx] = size(I);
x_cm =  ((1:nx)*downsample - c0(1)) * mmpp / 10;
z_cm = -((1:ny)*downsample - c0(2)) * mmpp / 10;   % row 1 (top) -> largest z

%% ---------------- Time since wind start ----------------
H = load(hw_file);
if isfield(H, 'camStartElapsed') && ~isnan(H.camStartElapsed)
    camDelay = H.camStartElapsed - H.fanStartElapsed;
else
    camDelay = H.runConfig.DELAY_BEFORE_TRIG;
    warning('No camStartElapsed in %s -- using DELAY_BEFORE_TRIG = %g s.', hw_file, camDelay);
end
t_frame = camDelay + (frame - 1) / fs;
fprintf('Frame %d: t = %.3f s since wind start (cameras started at %.3f s)\n', ...
        frame, t_frame, camDelay);

%% ---------------- Dye front ----------------
[z_front, z_surf, thr_col] = lif_dye_front(I, z_cm, mmpp*downsample, dy);
depth = z_surf - z_front;
if isempty(dy.dye_thr)
    lab = sprintf('dye front (%.0f\\%% recovery)', 100*dy.dye_frac);
else
    lab = sprintf('dye front ($<%g$ counts)', dy.dye_thr);
end
fprintf('Dye front: depth below surface %.1f cm mean, %.1f..%.1f cm range; threshold %.0f..%.0f counts\n', ...
        mean(depth, 'omitnan'), min(depth), max(depth), min(thr_col), max(thr_col));

%% ---------------- Plot ----------------
figure('Name', sprintf('CoreView_%d %s frame %d', run_num, cam, frame), ...
       'Position', [60 60 1300 700], 'Color', 'w');
imagesc(x_cm, z_cm, I);                        % unwarped: image as recorded
axis image; set(gca, 'YDir', 'normal');
colormap(gca, cmap);
if isempty(clims), clims = prctile(I(:), pct_lims); end
caxis(clims);
c = colorbar; c.Label.String = 'counts';
hold on;
plot(x_cm, z_front, 'r-', 'LineWidth', 1.5);
plot(x_cm, z_surf, 'c--', 'LineWidth', 1);
hold off;
legend({lab, 'surface'}, 'Interpreter', 'latex', 'Location', 'southeast');
xlabel('$x$ (cm)', 'Interpreter', 'latex');
ylabel('$z$ (cm, from board centre)', 'Interpreter', 'latex');
title(sprintf('CoreView %d, Water\\_SURF, frame %d, $t = %.2f$ s since wind start', ...
      run_num, frame, t_frame), 'Interpreter', 'latex');
set(gca, 'FontSize', 13, 'FontName', 'times');
