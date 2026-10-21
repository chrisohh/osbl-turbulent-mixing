%% CoreView 105 transverse LIF video: raw | C | profile, with a z99 time series
% The time-series counterpart of trans_frame_calibrated.m -- same helpers,
% same figure (util/lif_trans_figure + lif_trans_update), one video frame
% per camera frame, so a still pulled from the movie looks exactly like the
% snapshot.  See trans_frame_calibrated.m for what each setting does and
% why; the ones that matter only for a sequence are:
%
%   * the surface datum and the sheet profile S(y) are both taken ONCE,
%     from ref_frame, and held for the whole run.  They are properties of
%     the optics and of the undisturbed initial dye layer, not of the
%     frame, and re-deriving them per frame would make C drift as the dye
%     disperses -- exactly the thing the video is meant to show.
%   * frames are read in order and only the CData is replaced, so the
%     figure is built once.
%   * the run fades out rather than stopping: the dye disperses and
%     bleaches, and by frame ~2600 the sheet signal is a few percent of
%     its early level (frame 3249 is pure sensor noise).  frame_end is set
%     accordingly -- check it against the sheet-level printout below.
%
% The z99 / penetration time series is saved next to the video
% (<out_file>_z99.mat) for the z99(t) plot.
clear; close all;
addpath(fullfile(fileparts(mfilename('fullpath')), 'util'));

%% ---------------- User settings ----------------
raw_root  = 'D:\HLAB_2026\Core';
cam       = 'Water_PIV';
run_num   = 105;
ref_frame = 1;             % undisturbed dye sheet: surface datum AND S(y)

frame_start = 1;
frame_step  = 5;           % every Nth acquired frame
frame_end   = 2600;        % [] = last available; the dye has faded well before then

% ---- Scale (see trans_frame_calibrated.m: the plate is yawed, so only the
%      vertical scale is trustworthy) --------------------------------------
plate_num = 91;  plate_frm = 1;
sc.square_mm = 11;
sc.offset_mm = 0;          % plate-to-sheet offset toward the camera (CHECK)
sc.n_medium  = 1.33;
sc.lens_f_mm = 35;         % <-- CHECK
sc.mmpp      = [];         % [] = from the plate; a number overrides it

hw_file = '';              % '' = no hot-wire file for this run yet
fs      = 50;              % camera trigger rate (Hz) -- CHECK for run 105

show_hw  = ~isempty(hw_file);
hw_avg_s = 0.02;
hw_xlim  = [];

downsample = 2;
col_crop = [1 3900];       % drop the bright reflection at the right edge
row_crop = [1 3072];       % keep all rows: lif_trans_conc needs the clear water

% ---- Figure ----------------------------------------------------------
fo.cmap        = 'hot';
fo.raw_log     = true;
fo.clims_log   = [0.5 300];
fo.C_log       = true;
fo.C_clims_log = [3e-3 1];
fo.C_label     = '$C/C_{ref}$';
fo.profile_xlim = [0 0.45];   % FIXED for a video: an autoscaled x axis makes
                              % the profile jump about from frame to frame
fo.profile_zlim = [-4 0.5];
fo.zlim_cm      = [-4 1.5];
fo.ylim_cm      = [];
fo.fig_pos      = [40 40 1700 760];
fo.visible      = 'off';

% ---- Surface ---------------------------------------------------------
su.smooth_mm      = 1.0;
su.surf_smooth_cm = 3;
su.search_cm      = [];
su.min_counts     = 15;
su.fit_order      = 1;
su.track          = false;  % hold the frame-1 line; see trans_frame_calibrated.m
su.offset_cm      = 0;

% ---- Concentration + z99 ---------------------------------------------
co.frac        = 0.99;
co.smooth_mm   = su.smooth_mm;
co.surf_gap_cm = 0;
co.top_cm      = 0.3;
co.deep_cm     = 6;
co.deep_min_rows = 200;
co.min_peak    = 15;
co.int_depth_cm = 5;
co.normalize   = 'ref';
co.bin_cm      = 1;
co.min_bin_frac = 0.8;
co.pen_frac    = 0.02;
co.pen_min_counts = 5;
co.pen_smooth_cm  = 1.0;
co.dark        = [];

out_file      = fullfile(raw_root, sprintf('CoreView_%d', run_num), ...
                         sprintf('CoreView_%d_trans_conc.mp4', run_num));
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

%% ---------------- Scale, axes ----------------
cal_dir = fullfile(fileparts(mfilename('fullpath')), 'calibration');
fplate  = fullfile(raw_root, sprintf('CoreView_%d', plate_num), cam, ...
                   sprintf('CoreView_%d_%s_%02d.raw', plate_num, cam, plate_frm));
[mmpp, c0] = lif_trans_scale_cached(fplate, sc, ...
                 fullfile(cal_dir, sprintf('scale_%s_Core%d_trans.mat', cam, plate_num)));

crop = @(A) A(row_crop(1):downsample:row_crop(2), col_crop(1):downsample:col_crop(2));
Iref = crop(lif_load_raw(fname(ref_frame)));
[ny, nx] = size(Iref);
dx_mm = mmpp * downsample;
y_cm =  ((col_crop(1) + (0:nx-1)*downsample) - c0(1)) * mmpp / 10;
z_cm = -((row_crop(1) + (0:ny-1)*downsample) - c0(2)) * mmpp / 10;
if ~isempty(co.dark) && ischar(co.dark), co.dark = crop(lif_load_raw(co.dark)); end

%% ---------------- Surface datum + sheet profile, once ----------------
[z_ref, Sref] = lif_trans_surface(Iref, z_cm, dx_mm, su);
z0    = mean(z_ref, 'omitnan');
z_cm  = z_cm - z0;
z_ref = z_ref - z0;
su.prior = Sref;
co.ref   = Iref;
fprintf('Datum: frame %d water level -> z = 0 (line RMS %.1f px, camera roll %.2f deg)\n', ...
        ref_frame, Sref.rms_px, Sref.roll_deg);

%% ---------------- Time ----------------
if ~isempty(hw_file)
    [~, hw_name] = fileparts(hw_file);
    [camDelay, hw] = lif_hw_cached(hw_file, hw_avg_s, fullfile(cal_dir, ['hwsync_' hw_name '.mat']));
    t_wind = camDelay + (frames - 1) / fs;
    t_what = 'since wind start';
else
    hw = [];  camDelay = NaN;
    t_wind = (frames - 1) / fs;
    t_what = 'after frame 1';
    warning(['No hot-wire file for CoreView %d: times are measured from frame 1, ' ...
             'not from wind start.  Set hw_file once run_experiment.m has one.'], run_num);
end

%% ---------------- Figure + writer ----------------
if show_hw
    fo.hw = hw;
    fo.hw.xlim = hw_xlim;
end
h = lif_trans_figure(y_cm, z_cm, fo);

v = VideoWriter(out_file, 'MPEG-4');
v.FrameRate = video_fps;
v.Quality   = video_quality;
open(v);
cleanupObj = onCleanup(@() close(v));

z_surf_all = nan(nFrames, nx);
z_pen_all  = nan(nFrames, nx);
surf_all   = nan(nFrames, nx);      % surface-band counts: how much dye is left
z99_all    = nan(nFrames, 1);
z99_y_all  = [];
y99_cm     = [];
ref_col    = [];

%% ---------------- Frames ----------------
tic;
for i = 1:nFrames
    n = frames(i);
    I = crop(lif_load_raw(fname(n)));

    [z_surf_all(i,:), ~] = lif_trans_surface(I, z_cm, dx_mm, su);
    R = lif_trans_conc(I, z_cm, z_surf_all(i,:), dx_mm, co);
    z99_all(i)    = R.z99;
    z_pen_all(i,:) = R.z_pen;
    surf_all(i,:)  = R.surf_col;
    if isempty(z99_y_all)
        z99_y_all = nan(nFrames, numel(R.jc));
        y99_cm    = y_cm(R.jc);
        ref_col   = R.ref_col;
    end
    z99_y_all(i,:) = R.z99_y;

    lif_trans_update(h, I, R, y_cm, z_surf_all(i,:), ...
        sprintf('CoreView %d, transverse LIF (%s), frame %d, $t = %.2f$ s %s', ...
                run_num, strrep(cam, '_', '\_'), n, t_wind(i), t_what), t_wind(i));

    rgb = print(h.fig, '-RGBImage', '-r0');                 % works with an invisible figure
    rgb = rgb(1:2*floor(end/2), 1:2*floor(end/2), :);       % H.264 needs even size
    writeVideo(v, rgb);

    if mod(i, 25) == 0 || i == nFrames
        el = toc;
        fprintf('  %d/%d  frame %d  t = %.2f s  z99 = %5.2f cm  surface dye %3.0f%% of frame %d  --  %.2f fps  ETA %.0f s\n', ...
                i, nFrames, n, t_wind(i), z99_all(i), ...
                100*median(R.surf_col ./ R.ref_col, 'omitnan'), ref_frame, ...
                i/el, (nFrames-i)/(i/el));
    end
end
clear cleanupObj;   % closes the writer
fprintf('\nWrote %s  (%d frames, %.1f s)\n', out_file, nFrames, toc);

%% ---------------- Save the z99 time series ----------------
% Frame ref_frame is the reference itself: C there is 1 by construction.
co = rmfield(co, intersect({'ref', 'dark'}, fieldnames(co)));   % don't save images
co.ref_frame = ref_frame;
depth_all = z_surf_all - z_pen_all;        % [nFrames x nx] dye penetration, cm
[p, b] = fileparts(out_file);
ts_file = fullfile(p, [b '_z99.mat']);
save(ts_file, 'frames', 't_wind', 't_what', 'camDelay', 'fs', 'y_cm', ...
     'z99_all', 'z99_y_all', 'y99_cm', 'z_surf_all', 'z_pen_all', 'depth_all', ...
     'surf_all', 'ref_col', 'mmpp', 'c0', 'downsample', 'row_crop', 'col_crop', ...
     'su', 'co', 'sc', 'hw_file');
fprintf('Saved z99 time series to %s\n', ts_file);
