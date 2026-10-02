%% Water_SURF video: raw | theta | profile, with z99 time series
% Same three-panel figure as surf_frame_calibrated.m (util/lif_surf_figure
% + lif_surf_update), one video frame per camera frame:
%   1. raw counts as recorded, dye front / surface / local z99(x)
%   2. theta = ln(I_bg/I) (dye concentration, lighting removed)
%   3. all-x averaged theta vs depth, z99 and the 5-95% spread of z99(x)
% Title: time since wind start, t = camDelay + (n-1)/fs, from the hot-wire
% file.  The z99 / dye-front time series is saved next to the video
% (<out_file>_z99.mat) for the z99(t) plot.
clear; close all;
addpath(fullfile(fileparts(mfilename('fullpath')), 'util'));

%% ---------------- User settings ----------------
raw_root  = 'D:\HLAB_2026\Core';
cam       = 'Water_Surf';
run_num   = 95;

frame_start = 1;
frame_step  = 5;           % every Nth acquired frame
frame_end   = [];          % [] = last available

plate_num = 94;  plate_frm = 1;
sc.plate_scale = 1.0;
sc.square_mm   = 11;       % small plate: coarse 11x10 layout with 11 mm squares
sc.offset_mm   = 25;       % plate was this far TOWARD the camera
sc.n_medium    = 1.33;
sc.lens_f_mm   = 35;       % <-- CHECK

hw_file = 'D:\HLAB_2026\hotwire\hotwire_20260923_113939.mat';
fs      = 50;

downsample = 2;

fo.cmap            = 'bone';    % see lif_surf_figure
fo.left            = 'transmission';   % left panel: 'transmission' I/I_bg (lighting removed) | 'raw'
fo.trans_clims     = [0.75 1.05];
fo.theta_clims     = [0 0.05];
fo.theta_log       = false;
fo.theta_clims_log = [0.005 0.3];
fo.profile_xlim    = [-0.01 0.3];
fo.profile_zlim    = [-15 0];
fo.fig_pos         = [40 40 2000 650];
fo.visible         = 'off';

dy.dye_frac       = 0.90;  % see lif_dye_front
dy.dye_thr        = [];
dy.smooth_mm      = 1.5;
dy.surf_window_cm = [4 10];
dy.surf_gap_cm    = 0.4;

zo.frac          = 0.99;   % see lif_z99
zo.smooth_mm     = dy.smooth_mm;
zo.surf_gap_cm   = dy.surf_gap_cm;
zo.fit_order     = 2;
zo.fit_margin_cm = 0.5;
zo.bin_cm        = 1;
bg_frame         = 1;      % dye-free frame -> theta = ln(I_bg/I)

out_file      = fullfile(raw_root, sprintf('CoreView_%d', run_num), ...
                         sprintf('CoreView_%d_surf_theta.mp4', run_num));
video_fps     = 10;
video_quality = 95;

%% ---------------- Frames ----------------
cam_dir = fullfile(raw_root, sprintf('CoreView_%d', run_num), cam);
fname   = @(n) fullfile(cam_dir, sprintf('CoreView_%d_%s_%04d.raw', run_num, cam, n));
d = dir(fullfile(cam_dir, sprintf('CoreView_%d_%s_*.raw', run_num, cam)));
d = d([d.bytes] >= 4096*3072*2);                 % skip files still being copied
avail = sort(cellfun(@(s) str2double(regexp(s, '_(\d+)\.raw$', 'tokens', 'once')), {d.name}));
if isempty(avail), error('No complete frames in %s', cam_dir); end
if isempty(frame_end), frame_end = max(avail); end
frames  = intersect(frame_start:frame_step:frame_end, avail);
nFrames = numel(frames);
fprintf('%s: %d complete frames (%d..%d); writing %d -> %s\n', ...
        cam, numel(avail), avail(1), avail(end), nFrames, out_file);

%% ---------------- Scale, axes, background, time ----------------
fplate = fullfile(raw_root, sprintf('CoreView_%d', plate_num), cam, ...
                  sprintf('CoreView_%d_%s_%02d.raw', plate_num, cam, plate_frm));
[mmpp, c0] = lif_plate_scale(fplate, sc);

ny = 3072/downsample;  nx = 4096/downsample;
x_cm =  ((1:nx)*downsample - c0(1)) * mmpp / 10;
z_cm = -((1:ny)*downsample - c0(2)) * mmpp / 10;

zo.bg = [];
if ~isempty(bg_frame)
    zo.bg = lif_load_raw(fname(bg_frame));
    zo.bg = zo.bg(1:downsample:end, 1:downsample:end);
    fprintf('Background: frame %d (dye-free)\n', bg_frame);
    [z0, z_cm, dy] = lif_surface_datum(zo.bg, z_cm, mmpp*downsample, dy);   % z = 0 at the bg-frame water level
else
    warning('No bg_frame: z stays measured from the plate (board centre).');
end

H = load(hw_file);
if isfield(H, 'camStartElapsed') && ~isnan(H.camStartElapsed)
    camDelay = H.camStartElapsed - H.fanStartElapsed;
else
    camDelay = H.runConfig.DELAY_BEFORE_TRIG;
    warning('No camStartElapsed in %s -- using DELAY_BEFORE_TRIG = %g s.', hw_file, camDelay);
end
t_wind = camDelay + (frames - 1) / fs;

%% ---------------- Figure + writer ----------------


v = VideoWriter(out_file, 'MPEG-4');
v.FrameRate = video_fps;
v.Quality   = video_quality;
open(v);
cleanupObj = onCleanup(@() close(v));

z_front_all = nan(nFrames, nx);
z_surf_all  = nan(nFrames, nx);
z99_all     = nan(nFrames, 1);              % all-x average, cm below the surface
z99_x_all   = [];                           % [nFrames x nbin] local z99(x)
x99_cm      = [];

%% ---------------- Frames ----------------
tic;
for i = 1:nFrames
    n = frames(i);
    I = lif_load_raw(fname(n));
    I = I(1:downsample:end, 1:downsample:end);

    [z_front_all(i,:), z_surf_all(i,:)] = lif_dye_front(I, z_cm, mmpp*downsample, dy);
    [z99_x, jc, z99_all(i), zeta, theta_bar, theta_img, T_img] = ...
        lif_z99(I, z_cm, z_surf_all(i,:), z_front_all(i,:), mmpp*downsample, zo);
    if isempty(z99_x_all), z99_x_all = nan(nFrames, numel(jc)); x99_cm = x_cm(jc); end
    z99_x_all(i,:) = z99_x;

    if strcmp(fo.left, 'transmission'), L = T_img; else, L = I; end
    lif_surf_update(h, L, theta_img, x_cm, z_front_all(i,:), z_surf_all(i,:), ...
        z99_x, jc, zeta, theta_bar, z99_all(i), ...
        sprintf('CoreView %d, Water\\_SURF, frame %d, $t = %.2f$ s since wind start', ...
                run_num, n, t_wind(i)), t_wind(i));

    rgb = print(h.fig, '-RGBImage', '-r0');                 % works with an invisible figure
    rgb = rgb(1:2*floor(end/2), 1:2*floor(end/2), :);       % H.264 needs even size
    writeVideo(v, rgb);

    if mod(i, 25) == 0 || i == nFrames
        el = toc;
        fprintf('  %d/%d  frame %d  t = %.2f s  z99 = %.2f cm  --  %.2f fps  ETA %.0f s\n', ...
                i, nFrames, n, t_wind(i), z99_all(i), i/el, (nFrames-i)/(i/el));
    end
end
clear cleanupObj;   % closes the writer
fprintf('\nWrote %s  (%d frames, %.1f s)\n', out_file, nFrames, toc);

%% ---------------- Save the z99 time series ----------------
% Frame bg_frame is the dye-free reference itself: its z99 is noise.
depth_all = z_surf_all - z_front_all;       % [nFrames x nx] cm below surface
zo = rmfield(zo, intersect({'bg', 'dark'}, fieldnames(zo)));   % don't save images
zo.bg_frame = bg_frame;  zo.dark_file = dark_file;
[p, b] = fileparts(out_file);
ts_file = fullfile(p, [b '_z99.mat']);
save(ts_file, 'frames', 't_wind', 'x_cm', 'z99_all', 'z99_x_all', 'x99_cm', ...
     'z_front_all', 'z_surf_all', 'depth_all', 'mmpp', 'c0', 'downsample', ...
     'dy', 'zo', 'sc', 'hw_file');
fprintf('Saved z99 time series to %s\n', ts_file);
