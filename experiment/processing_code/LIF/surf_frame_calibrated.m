%% Water_SURF frame: raw | theta | profile, timed against wind start
% Longitudinal (side) view, CoreView 95.  Images are shown AS RECORDED (no
% rectification); the checkerboard plate in CoreView 94 (same camera, same
% mount, 2.5 cm toward the camera) only sets the mm/px scale -- see
% util/lif_plate_scale.m.  Dye front: util/lif_dye_front.m.  Concentration
% theta = ln(I_bg/I) and z99 (paper eq. 3.14, local in x): util/lif_z99.m.
% Figure: util/lif_surf_figure.m + lif_surf_update.m, shared with
% surf_theta_video.m so the still and the video match.
%
% Time: frame n is trigger n of the 50 Hz camera counter, so
%     t = (camStartElapsed - fanStartElapsed) + (n - 1) / fs
% seconds since wind start, from the run_experiment.m hot-wire file.
clear; close all;
addpath(fullfile(fileparts(mfilename('fullpath')), 'util'));

%% ---------------- User settings ----------------
raw_root  = 'C:\Users\Administrator\Downloads\';
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

hw_file = 'C:\Users\Administrator\Downloads\hotwire_20260923_113939.mat';
fs      = 50;

% Hot-wire panel under the images, vertical line at the frame time
show_hw  = true;
hw_avg_s = 0.02;           % block average (s); 0.02 = one camera frame
hw_xlim  = [];             % s since wind start; [] = whole record

downsample = 2;            % every Nth pixel

% Figure (see lif_surf_figure).  At z99 the dye darkens the water by only
% ~1% (~6 counts), below the ~100-count lighting gradient, so it shows in
% the theta panel, not the raw one.
fo.cmap            = 'bone';
fo.left            = 'raw';%'transmission';   % left panel: 'transmission' I/I_bg (lighting removed) | 'raw'
fo.trans_clims     = [0.75 1.05];
fo.theta_clims     = [0 0.25];      % theta panel, linear: saturates the core, shows the edge
fo.theta_log       = false;         % true: log colour scale (core AND edge)
fo.theta_clims_log = [0.005 0.3];
fo.profile_xlim    = [-0.01 0.3];
fo.profile_zlim    = [-15 0];       % cm below the surface

% Dye front (see lif_dye_front)
dy.dye_frac       = 0.90;  % 0.9 ~ the hand-drawn line; 0.99 = outermost edge
dy.dye_thr        = [];    % counts; set e.g. 530 for a fixed threshold instead
dy.smooth_mm      = 1.5;
dy.surf_window_cm = [4 10];
dy.surf_gap_cm    = 0.4;

% z99 (see lif_z99): level above which 99% of the dye resides (eq. 3.14)
zo.frac          = 0.99;
zo.smooth_mm     = dy.smooth_mm;
zo.surf_gap_cm   = dy.surf_gap_cm;
zo.fit_order     = 2;      % per-column clear-water fit, used only without bg
zo.fit_margin_cm = 0.5;    % clear water starts this far below the dye front
zo.bin_cm        = 1;      % x bin for the local z99(x) (longitudinal view)
bg_frame         = 1;      % dye-free Water_SURF frame -> theta = ln(I_bg/I); [] = per-column fit

%% ---------------- Scale ----------------
fplate = fullfile(raw_root, sprintf('CoreView_%d', plate_num), cam, ...
                  sprintf('CoreView_%d_%s_%02d.raw', plate_num, cam, plate_frm));
[mmpp, c0] = lif_plate_scale(fplate, sc);

%% ---------------- Load frame + background ----------------
fname = @(n) fullfile(raw_root, sprintf('CoreView_%d', run_num), cam, ...
                      sprintf('CoreView_%d_%s_%04d.raw', run_num, cam, n));
I = lif_load_raw(fname(frame));
I = I(1:downsample:end, 1:downsample:end);
[ny, nx] = size(I);
x_cm =  ((1:nx)*downsample - c0(1)) * mmpp / 10;
z_cm = -((1:ny)*downsample - c0(2)) * mmpp / 10;   % row 1 (top) -> largest z

zo.bg = [];
if ~isempty(bg_frame)
    zo.bg = lif_load_raw(fname(bg_frame));
    zo.bg = zo.bg(1:downsample:end, 1:downsample:end);
    fprintf('Background: frame %d (dye-free)\n', bg_frame);
    [z0, z_cm, dy] = lif_surface_datum(zo.bg, z_cm, mmpp*downsample, dy);   % z = 0 at the bg-frame water level
else
    warning('No bg_frame: z stays measured from the plate (board centre).');
end

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

%% ---------------- Dye front + z99 ----------------
[z_front, z_surf, thr_col] = lif_dye_front(I, z_cm, mmpp*downsample, dy);
depth = z_surf - z_front;
fprintf('Dye front: depth below surface %.1f cm mean, %.1f..%.1f cm range; threshold %.0f..%.0f counts\n', ...
        mean(depth, 'omitnan'), min(depth), max(depth), min(thr_col), max(thr_col));

[z99_x, jc, z99, zeta, theta_bar, theta_img, T_img] = lif_z99(I, z_cm, z_surf, z_front, mmpp*downsample, zo);
fprintf('z99(x) = %.2f..%.2f cm (mean %.2f) below the surface, %g cm bins; all-x average z99 = %.2f cm\n', ...
        min(z99_x), max(z99_x), mean(z99_x, 'omitnan'), zo.bin_cm, z99);

%% ---------------- Plot ----------------
if show_hw
    fo.hw = lif_hw_series(H, hw_avg_s);
    fo.hw.xlim = hw_xlim;
end
h = lif_surf_figure(x_cm, z_cm, fo);
set(h.fig, 'Name', sprintf('CoreView_%d %s frame %d', run_num, cam, frame));
if strcmp(fo.left, 'transmission'), L = T_img; else, L = I; end
lif_surf_update(h, L, theta_img, x_cm, z_front, z_surf, z99_x, jc, zeta, theta_bar, z99, ...
    sprintf('CoreView %d, Water\\_SURF, frame %d, $t = %.2f$ s since wind start', run_num, frame, t_frame), t_frame);
