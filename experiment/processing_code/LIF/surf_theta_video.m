%% Water_SURF video: raw | theta | profile, with z99 time series
% Same three-panel figure as surf_frame_calibrated.m (util/lif_surf_figure
% + lif_surf_update), one video frame per camera frame:
%   1. raw counts as recorded, dye front / surface / local z99(x)
%   2. theta = ln(I_bg/I) (dye concentration, lighting removed)
%   3. all-x averaged theta vs depth, z99 and the 5-95% spread of z99(x)
% For the stripped theta-over-hot-wire video, see theta_hw_video.m.
% Title: time since wind start, t = camDelay + (n-1)/fs, from the hot-wire
% file.  The z99 / dye-front time series is saved next to the video
% (<out_file>_z99.mat) for the z99(t) plot.
clearvars -except run_num run_queue; close all;   % set run_num before calling to pick the run
addpath(fullfile(fileparts(mfilename('fullpath')), 'util'));

%% ---------------- User settings ----------------
raw_root  = 'D:\HLAB_2026\Core';
if ~exist('run_num', 'var'), run_num = 95; end   % run settings: util/lif_run_config.m
R = lif_run_config(run_num);
cam = R.cam;

frame_start = 1;
frame_step  = 10;           % every Nth acquired frame
frame_end   = [];          % [] = last available

plate_num = R.plate_num;
sc.plate_scale = 1.0;
sc.square_mm   = R.square_mm;  % small plate: coarse 11x10 layout with 11 mm squares
sc.offset_mm   = R.offset_mm;  % plate this far TOWARD the camera from the tank centre
sc.n_medium    = 1.33;
sc.lens_f_mm   = R.lens_f_mm;  % Nikon 60 mm

hw_file = R.hw_file;
fs      = 50;

% Hot-wire panel under the images, vertical line at the frame time
show_hw  = true;
hw_avg_s = 0.02;           % block average (s); 0.02 = one camera frame
hw_xlim  = [];             % s since wind start; [] = whole record

downsample = 2;

fo.cmap            = 'bone';    % see lif_surf_figure
fo.img_xlim        = [-25 22];  % image panels, cm
fo.img_ylim        = [-22 7];
zg.floor_win       = [20 40];   % no z99(x) line until the run's z99 leaves its own surface-layer
zg.floor_fac       = 1.25;      % plateau (median over floor_win) by floor_fac / 3 sigma / 0.25 cm;
zg.fixed           = [];        % see util/lif_z99_gate.  A number (e.g. -2.3, run 95) = old fixed floor
fo.left            = 'raw';%'transmission';   % left panel: 'transmission' I/I_bg (lighting removed) | 'raw'
fo.trans_clims     = [0.75 1.05];
fo.theta_clims     = [0 0.25];
fo.theta_log       = false;
fo.theta_clims_log = [0.005 0.3];
fo.profile_xlim    = [-0.01 0.3];
fo.profile_zlim    = [];        % set by lif_surf_figure to the theta image's z range
fo.z99_band        = [5 95];   % z99(x) spread shaded in the profile: [25 75] steadier, [] = none
fo.fig_pos         = [40 40 2000 650];
fo.visible         = 'on';

dy.dye_frac       = 0.90;  % see lif_dye_front
dy.dye_thr        = [];
dy.smooth_mm      = 0;      % (only for filter 'gauss')
dy.filter         = 'median3';  % image noise filter: 3x3 median keeps edges and amplitude (lif_img_filter)
dy.surf_method    = 'edge';     % surface = top of the dark band, as the side view shows it
dy.surf_median_cm = 0;      % no running median on the surface line
dy.surf_window_cm = [4 10];
dy.surf_gap_cm    = 0.4;

zo.frac          = 0.99;   % see lif_z99
zo.smooth_mm     = dy.smooth_mm;
zo.surf_gap_cm   = 0.15;   % top band left out of z99 (meniscus); small since bg_shift aligns the surface,
                           % and early on the whole dye layer is only ~1 cm thick
zo.fit_order     = 2;
zo.fit_margin_cm = 0.5;
line_bin         = 'auto';  % z99(x) line window (cm): 'auto' = line_k / dye SNR of the run (worst of
                           % 3 late frames, util/lif_dye_snr), 0.5 cm steps in line_lim; 'registry' =
                           % R.line_bin_cm; or a number.  SNR 118 ~17 -> 1 cm, 120 ~6 -> 3 cm, 96 ~2 -> 4 cm
line_k           = 18;
line_lim         = [1 4];
zo.bin_cm        = R.line_bin_cm;   % (replaced below unless line_bin = 'registry'; 0 = per column)
zo.step_cm       = 0;      % (0 with bin 0: every column; with a window: one value per column too)
zo.xmed          = 11;     % median along x over ~2.5 mm: removes z99(x) spikes
zo.band_fill     = 'none';  % DISPLAY only: 'below' paints the surface band (the dark line) with
                           % the mixed-layer theta under it; 'black' = black like the dark
                           % line; 'none' = as measured.  z99 is unaffected (see lif_z99)
zo.filter        = dy.filter;
zo.surf_method   = dy.surf_method;
zo.z99_all       = 'profile';  % frame z99 from the all-x averaged profile (paper eq. 3.14):
                           % the only estimate that stays monotonic for faint runs (96)
zo.snr_min       = 4;      % no z99 unless theta clears this many sigma of the noise
zo.tail_cm       = 1;      % integral stops after this much quiet water (see lif_z99)
bg_frames        = R.bg_frames;   % dye-free frames, averaged -> theta = ln(I_bg/I).
                           % ALL must be dye-free: check when the dye arrives (per run).
dark_file        = [];     % lens-capped frame, same exposure (.raw path): subtracted from I and I_bg; [] = none
datum_tol_cm     = 1;      % move z = 0 if the early-run surface is further than this from the bg level

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
fplate = R.plate_file;
cal_dir = fullfile(fileparts(mfilename('fullpath')), 'calibration');   % caches for the other computer
[mmpp, c0] = lif_scale_cached(fplate, sc, fullfile(cal_dir, sprintf('scale_%s_Core%d.mat', cam, plate_num)));

ny = 3072/downsample;  nx = 4096/downsample;
x_cm =  ((1:nx)*downsample - c0(1)) * mmpp / 10;
z_cm = -((1:ny)*downsample - c0(2)) * mmpp / 10;

zo.bg = [];  relight = struct('applied', false);  datum_shift = 0;
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
    [z0, z_cm, dy] = lif_surface_datum(zo.bg, z_cm, mmpp*downsample, dy);   % z = 0 at the bg-frame water level
    [zo.bg, relight] = lif_bg_relight(zo.bg, fname, avail, bg_frames, downsample, z_cm, struct('fs', fs));
    if relight.applied                       % lighting stepped after the bg frames (102: t ~ 14.5 s)
        fprintf('Lighting step at frame %d: bg re-lit from frames %d..%d (gain %.3f, 5-95%% %.3f..%.3f)\n', ...
                relight.step_frame, relight.frames_used(1), relight.frames_used(end), relight.gain, relight.gain_range);
    end
    % Datum check: the water surface 1-4 s after the bg frames (wind on,
    % waves still mm-sized) must sit at the bg water level.  102: the bg
    % frames' edge is a feature 2.6 cm above the water that vanishes with
    % the first ripples (12.4 s; the camera did not move) -> z = 0 moves to
    % the early-run surface.  Depths below the local surface are unaffected.
    ck = intersect(max(bg_frames) + round((1:0.5:4) * fs), avail);
    zck = nan(size(ck));
    for k = 1:numel(ck)
        I = lif_load_raw(fname(ck(k)));  I = I(1:downsample:end, 1:downsample:end);
        [~, zs_k] = lif_dye_front(I, z_cm, mmpp*downsample, dy);
        zck(k) = median(zs_k, 'omitnan');
    end
    datum_shift = median(zck, 'omitnan');
    if abs(datum_shift) > datum_tol_cm
        z_cm = z_cm - datum_shift;
        dy.surf_window_cm = dy.surf_window_cm - datum_shift;
        % lif_z99 windows that are heights above the BG water level stay
        % there (102: the air reference at 2-6 cm would otherwise take in
        % the bg's own dark edge and add theta ~ +0.05 everywhere)
        zw = struct('air_window_cm', [2 6], 'bg_window_cm', [-1 1]);   % lif_z99 defaults
        for fz = fieldnames(zw).'
            if ~isfield(zo, fz{1}), zo.(fz{1}) = zw.(fz{1}); end
            zo.(fz{1}) = zo.(fz{1}) - datum_shift;
        end
        warning('surf_theta_video:datum', ['Run %d: surface %.2f cm from the bg water level 1-4 s after ' ...
                'the bg frames -> z = 0 moved to the early-run surface.'], run_num, datum_shift);
    end
else
    warning('No bg_frames: z stays measured from the plate (board centre).');
end

[camDelay, hw, wave] = lif_wind_clock(R, hw_avg_s, cal_dir);   % hot-wire / ref probe / wave file, or assumed start + stand-in; wave gauge (own or stand-in)
t_wind = camDelay + (frames - 1) / fs;

%% ---------------- Figure + writer ----------------
if show_hw && ~isempty(hw)
    fo.hw = hw;
    fo.hw.xlim = hw_xlim;
    if ~isempty(wave), fo.wave = wave; end   % wave-gauge panel under the wind
end
fo.ct = lif_ct_profile(R.density_file);   % stratified runs: pre-run CT profile in the profile panel ([] otherwise)
if ~isempty(fo.ct)
    fprintf('CT profile %s: pycnocline %.1f..%.1f cm\n', fo.ct.file, fo.ct.pyc);
end
fo.z99_floor_cm = NaN;          % per frame from the gate below
h = lif_surf_figure(x_cm, z_cm, fo);
G = lif_z99_gate(zg);

vid_file = lif_video_file(out_file);   % time-stamped name if the old video is open somewhere
v = VideoWriter(vid_file, 'MPEG-4');
v.FrameRate = video_fps;
v.Quality   = video_quality;
open(v);
cleanupObj = onCleanup(@() close(v));

z_front_all = nan(nFrames, nx);
z_surf_all  = nan(nFrames, nx);
z99_all     = nan(nFrames, 1);              % all-x average, cm below the surface
z99_x_all   = [];                           % [nFrames x nbin] local z99(x)
x99_cm      = [];
z90_all     = nan(nFrames, 1);              % 90% depth (median of bins), cm
zrms_all    = nan(nFrames, 1);              % RMS dye depth (median of bins), cm
peak_all    = nan(nFrames, 1);              % max of the all-x theta profile
deep_all    = nan(nFrames, 1);              % fraction of x with dye at the frame bottom (z99 beyond view)

%% ---------------- z99(x) line window from the dye SNR ----------------
if isnumeric(line_bin)
    zo.bin_cm = line_bin;
elseif strcmp(line_bin, 'auto')
    % late frames: the dye is most diluted there, so they set the worst SNR
    pf = frames(max(1, round(numel(frames) * [0.6 0.8 1])));
    snr = nan(size(pf));
    for k = 1:numel(pf)
        I = lif_load_raw(fname(pf(k)));  I = I(1:downsample:end, 1:downsample:end);
        [zf_k, zs_k] = lif_dye_front(I, z_cm, mmpp*downsample, dy);
        [~, ~, ~, ~, ~, th_k] = lif_z99(I, z_cm, zs_k, zf_k, mmpp*downsample, zo);
        snr(k) = lif_dye_snr(th_k, z_cm, zs_k);
    end
    zo.bin_cm = min(max(round(2 * line_k / min(snr)) / 2, line_lim(1)), line_lim(2));
    fprintf('z99(x) line window %.1f cm (dye SNR %s at frames %s; registry %.1f cm)\n', ...
            zo.bin_cm, mat2str(round(snr, 1)), mat2str(pf), R.line_bin_cm);
end

%% ---------------- Frames ----------------
tic;
for i = 1:nFrames
    n = frames(i);
    I = lif_load_raw(fname(n));
    I = I(1:downsample:end, 1:downsample:end);

    [z_front_all(i,:), z_surf_all(i,:)] = lif_dye_front(I, z_cm, mmpp*downsample, dy);
    [z99_x, jc, z99_all(i), zeta, theta_bar, theta_img, T_img, st] = ...
        lif_z99(I, z_cm, z_surf_all(i,:), z_front_all(i,:), mmpp*downsample, zo);
    if isempty(z99_x_all), z99_x_all = nan(nFrames, numel(jc)); x99_cm = x_cm(jc); end
    z99_x_all(i,:) = z99_x;
    z90_all(i) = st.z90;  zrms_all(i) = st.zrms;  peak_all(i) = st.peak;  deep_all(i) = st.deep_frac;

    [h.opts.z99_floor_cm, G] = lif_z99_gate(G, t_wind(i), z99_all(i));
    if strcmp(fo.left, 'transmission'), L = T_img; else, L = I; end
    lif_surf_update(h, L, theta_img, x_cm, z_front_all(i,:), z_surf_all(i,:), ...
        z99_x, jc, zeta, theta_bar, z99_all(i), ...
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
        fprintf('  %d/%d  frame %d  t = %.2f s  z99 = %.2f cm  --  %.2f fps  ETA %.0f s\n', ...
                i, nFrames, n, t_wind(i), z99_all(i), i/el, (nFrames-i)/(i/el));
    end
end
clear cleanupObj;   % closes the writer
fprintf('\nWrote %s  (%d frames, %.1f s)\n', vid_file, nFrames, toc);

%% ---------------- Save the z99 time series ----------------
% z99 is NaN while there is no detectable dye (bg_frames included: they are
% the reference itself), so the time series starts where the dye arrives.
depth_all = z_surf_all - z_front_all;       % [nFrames x nx] cm below surface
zo = rmfield(zo, intersect({'bg', 'dark'}, fieldnames(zo)));   % don't save images
zo.bg_frames = bg_frames;  zo.dark_file = dark_file;
[p, b] = fileparts(out_file);
ts_file = fullfile(p, [b '_z99.mat']);
save(ts_file, 'frames', 't_wind', 'x_cm', 'z99_all', 'z99_x_all', 'x99_cm', ...
     'z90_all', 'zrms_all', 'peak_all', 'deep_all', ...
     'z_front_all', 'z_surf_all', 'depth_all', 'mmpp', 'c0', 'downsample', ...
     'dy', 'zo', 'sc', 'hw_file', 'run_num', 'R', 'G', 'relight', 'datum_shift');   % G.t_on: when the z99 line came on
fprintf('Saved z99 time series to %s\n', ts_file);
