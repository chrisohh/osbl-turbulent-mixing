%% Water_SURF video: theta over the hot-wire velocity series
% Two panels, one video frame per camera frame:
%   top    theta = ln(I_bg/I) (dye concentration, lighting removed), full
%          width, with the detected surface.  No z99, no raw panel, no
%          horizontally averaged profile.
%   bottom hot-wire velocity vs time since wind start, red line at the
%          current frame.
% Same util/lif_surf_figure as surf_theta_video.m, with opts.panels
% 'theta' and opts.show_z99 false, so the theta image itself is identical
% to the one in that video.  Use surf_theta_video.m when you want the raw
% panel, the averaged profile and the saved z99 time series; this script
% saves no time series.
clear; close all;
addpath(fullfile(fileparts(mfilename('fullpath')), 'util'));

%% ---------------- User settings ----------------
raw_root  = 'D:\HLAB_2026\Core';
cam       = 'Water_Surf';
run_num   = 95;

frame_start = 1;
frame_step  = 10;          % every Nth acquired frame
frame_end   = [];          % [] = last available

plate_num = 94;  plate_frm = 1;
sc.plate_scale = 1.0;
sc.square_mm   = 11;       % small plate: coarse 11x10 layout with 11 mm squares
sc.offset_mm   = 25;       % plate was this far TOWARD the camera
sc.n_medium    = 1.33;
sc.lens_f_mm   = 35;       % <-- CHECK

hw_file = 'D:\HLAB_2026\hotwire\hotwire_20260923_113939.mat';
fs      = 50;

hw_avg_s = 0.02;           % block average (s); 0.02 = one camera frame
hw_xlim  = [];             % s since wind start; [] = whole record

downsample = 2;

fo.cmap            = 'bone';
fo.img_xlim        = [-25 22];  % image panels, cm
fo.img_ylim        = [-22 7];
fo.z99_floor_cm    = -2.3;      % no z99(x) line while the frame z99 is shallower (side-view floor, dye still at the surface)
fo.theta_clims     = [0 0.25];
fo.theta_log       = false;
fo.theta_clims_log = [0.005 0.3];
% The image keeps its true aspect (axis image) and the theta panel is now
% full width, so this figure wants to be TALLER than the three-panel one:
% image height ~ width/1.33 for a 4096x3072 frame, and the image gets 2/3
% of the figure height above the hot-wire panel.  Too wide and you get
% white margins either side of the image.
fo.fig_pos         = [40 40 1400 1200];
fo.visible         = 'on';
fo.panels          = 'theta';   % theta only, full width (see lif_surf_figure)
fo.show_z99        = false;     % no z99 overlay or legend entry

dy.dye_frac       = 0.90;  % see lif_dye_front
dy.dye_thr        = [];
dy.smooth_mm      = 0;      % (only for filter 'gauss')
dy.filter         = 'median3';  % image noise filter: 3x3 median keeps edges and amplitude (lif_img_filter)
dy.surf_method    = 'edge';     % surface = top of the dark band, as the side view shows it
dy.surf_median_cm = 0;      % no running median on the surface line
dy.surf_window_cm = [4 10];
dy.surf_gap_cm    = 0.4;

% theta still comes from lif_z99 (it returns theta_img), so its settings
% are still needed even though z99 itself is not plotted.
zo.frac          = 0.99;   % see lif_z99
zo.smooth_mm     = dy.smooth_mm;
zo.surf_gap_cm   = 0.15;   % top band left out of z99 (meniscus); small since bg_shift aligns the surface,
                           % and early on the whole dye layer is only ~1 cm thick
zo.fit_order     = 2;
zo.fit_margin_cm = 0.5;
zo.bin_cm        = 0;      % z99(x) per pixel column: no window
zo.step_cm       = 0;
zo.xmed          = 11;     % median along x over ~2.5 mm: removes z99(x) spikes
zo.filter        = dy.filter;
zo.surf_method   = dy.surf_method;
zo.snr_min       = 4;
zo.tail_cm       = 1;
bg_frames        = 1:10:100;   % dye-free frames, averaged -> theta = ln(I_bg/I).
                               % ALL must be dye-free: check when the dye arrives.
dark_file        = [];     % lens-capped frame, same exposure (.raw path); [] = none

out_file      = fullfile(raw_root, sprintf('CoreView_%d', run_num), ...
                         sprintf('CoreView_%d_theta_hw.mp4', run_num));
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
cal_dir = fullfile(fileparts(mfilename('fullpath')), 'calibration');   % shared caches
[mmpp, c0] = lif_scale_cached(fplate, sc, fullfile(cal_dir, sprintf('scale_%s_Core%d.mat', cam, plate_num)));

ny = 3072/downsample;  nx = 4096/downsample;
x_cm =  ((1:nx)*downsample - c0(1)) * mmpp / 10;
z_cm = -((1:ny)*downsample - c0(2)) * mmpp / 10;

zo.bg = [];
if ~isempty(bg_frames)
    bg_frames = intersect(bg_frames, avail);
    if isempty(bg_frames), error('None of the requested bg_frames exist in %s', cam_dir); end
    zo.bg = lif_load_bg(fname, bg_frames, downsample);
    fprintf('Background: mean of %d dye-free frames (%d..%d)\n', ...
            numel(bg_frames), min(bg_frames), max(bg_frames));
    zo.dark = [];
    if ~isempty(dark_file)
        zo.dark = lif_load_raw(dark_file);
        zo.dark = zo.dark(1:downsample:end, 1:downsample:end);
        fprintf('Dark frame: %s\n', dark_file);
    end
    [z0, z_cm, dy] = lif_surface_datum(zo.bg, z_cm, mmpp*downsample, dy);   % z = 0 at the bg water level
else
    warning('No bg_frames: z stays measured from the plate (board centre).');
end

[~, hw_name] = fileparts(hw_file);
[camDelay, hw] = lif_hw_cached(hw_file, hw_avg_s, fullfile(cal_dir, ['hwsync_' hw_name '.mat']));
t_wind = camDelay + (frames - 1) / fs;

%% ---------------- Figure + writer ----------------
fo.hw = hw;
fo.hw.xlim = hw_xlim;
h = lif_surf_figure(x_cm, z_cm, fo);

v = VideoWriter(out_file, 'MPEG-4');
v.FrameRate = video_fps;
v.Quality   = video_quality;
open(v);
cleanupObj = onCleanup(@() close(v));

%% ---------------- Frames ----------------
tic;
for i = 1:nFrames
    n = frames(i);
    I = lif_load_raw(fname(n));
    I = I(1:downsample:end, 1:downsample:end);

    [z_front, z_surf] = lif_dye_front(I, z_cm, mmpp*downsample, dy);
    [z99_x, jc, z99, zeta, theta_bar, theta_img, ~, st] = ...
        lif_z99(I, z_cm, z_surf, z_front, mmpp*downsample, zo);

    % The raw / profile / z99 arguments are ignored by a 'theta' figure,
    % but passing the real values keeps the call identical to the one in
    % surf_theta_video.m, so adding a panel back needs no other change.
    lif_surf_update(h, I, theta_img, x_cm, z_front, z_surf, ...
        z99_x, jc, zeta, theta_bar, z99, ...
        sprintf('CoreView %d, Water\\_SURF, frame %d, $t = %.2f$ s since wind start', ...
                run_num, n, t_wind(i)), t_wind(i), st);

    rgb = print(h.fig, '-RGBImage', '-r0');                 % works with an invisible figure
    if i == 1
        vsz = 2*floor([size(rgb,1) size(rgb,2)]/2);         % H.264 needs even size
    end
    rgb = lif_fit_frame(rgb, vsz);                          % print size can drift by a pixel
    writeVideo(v, rgb);

    if mod(i, 25) == 0 || i == nFrames
        el = toc;
        fprintf('  %d/%d  frame %d  t = %.2f s  --  %.2f fps  ETA %.0f s\n', ...
                i, nFrames, n, t_wind(i), i/el, (nFrames-i)/(i/el));
    end
end
clear cleanupObj;   % closes the writer
fprintf('\nWrote %s  (%d frames, %.1f s)\n', out_file, nFrames, toc);
