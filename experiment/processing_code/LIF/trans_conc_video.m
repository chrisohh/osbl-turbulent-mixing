%% Transverse LIF video: raw | C | profile (+ wind, wave, CT), with a z99 time series
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
clearvars -except run_num run_queue; close all;   % set run_num before calling to pick the run
addpath(fullfile(fileparts(mfilename('fullpath')), 'util'));

%% ---------------- User settings ----------------
% Run settings (camera, plate, reference frame, wind record) come from
% processing_code/common/run_registry.csv via util/lif_trans_config.m.
% Laser-sheet runs: 104, 105 (unstrat), 119, 121 (strat).
raw_root  = 'D:\HLAB_2026\Core';
if ~exist('run_num', 'var'), run_num = 105; end
Rt  = lif_trans_config(run_num);
cam = Rt.cam;
ref_frame = Rt.ref_frame;  % undisturbed dye sheet: surface datum AND S(y); NaN = pick automatically

frame_start = 1;
frame_step  = 10;          % every Nth acquired frame (0.2 s, as the side view)
frame_end   = [];          % [] = last available; z99 is gated off once the dye fades

% ---- Flattening onto the laser-sheet plane --------------------------------
% The camera sees the sheet ~58 deg off its axis.  The plate stood IN the
% sheet, so its corners give the full plane -> image map (util/
% lif_trans_rectify: the IR fisheye-correction model, homography only for
% this 60 mm lens), and every frame is resampled on a regular (y, z) grid
% in that plane -- as seen along the plate normal, true cm on both axes.
% The vertical-only scale used before squeezed y to about half.
plate_num = Rt.plate_num;
rc.square_mm = 11;
rc.dx_mm     = 0.125;      % output grid: the old 2x-downsampled vertical pitch;
                           % about the raw pitch across the tank (~0.12 mm/px)
rc.cols      = [1 3900];   % drop the bright reflection at the right edge

fs      = 50;              % camera trigger rate (Hz)

hw_avg_s = 0.02;
hw_xlim  = [];

% ---- Figure ----------------------------------------------------------
fo.cmap        = 'hot';
fo.raw_log     = true;
fo.clims_log   = [0.5 300];
fo.C_log       = true;
fo.C_clims_log = [3e-3 1];
fo.C_label     = '$C/C_{ref}$';
fo.profile_xlim = [0 0.45];   % FIXED for a video: an autoscaled x axis makes
                              % the profile jump about from frame to frame
fo.zlim_cm      = [-16 1.5];
fo.ylim_cm      = [];
fo.visible      = 'off';
fo.show_ds      = 2;        % display every 2nd pixel (~0.25 mm): finer than the panels show
fo.z99_band     = [5 95];   % z99(y) spread shaded in the profile ([] = none)
zg.floor_win    = [20 40];  % no z99(y) line until z99 leaves the run's own surface-layer plateau
zg.floor_fac    = 1.25;     % (util/lif_z99_gate, as surf_theta_video); a number in zg.fixed = fixed floor
zg.fixed        = [];

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
co.int_depth_cm = 15;      % integral stops this far below the surface: 5 cm (run 105's early
                           % frames) capped the deepening; 10-18 cm agree within ~0.1 cm
co.normalize   = 'ref';
co.bin_cm      = 1;
co.min_bin_frac = 0.8;
co.pen_frac    = 0.02;
co.pen_min_counts = 5;
co.pen_smooth_cm  = 1.0;
co.dark        = [];
sig_gate       = 0.10;     % z99 blank when the sheet signal (p99.9 counts) is below this
                           % fraction of the reference frame's: no dye yet (119/121 before
                           % ~34 s) or diluted/bleached late in the run -- z99 there is noise

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

%% ---------------- Flattening map, axes ----------------
cal_dir = fullfile(fileparts(mfilename('fullpath')), 'calibration');
fplate  = Rt.plate_file;   % plate frames _01.._05, first usable is taken
rc.cache_file = fullfile(cal_dir, sprintf('rectify_%s_Core%d_trans.mat', ...
                         regexprep(cam, '\W', ''), plate_num));
Tr = lif_trans_rectify(fplate, rc);

crop = @(A) lif_trans_apply(Tr, A);        % flatten (name kept from the cropping it replaced)
if isempty(ref_frame) || isnan(ref_frame)
    % Not set in the registry (119, 121: frame 1 is dark -- the dye is not
    % in the sheet yet).  Take the first frame whose bright tail reaches half
    % of its maximum over the first 30 s: the dye layer has arrived but the
    % wind has not mixed it yet.  Check it, then write it into the registry.
    cand = intersect(1:10:1500, avail);  b = nan(size(cand));
    for k = 1:numel(cand), J = crop(lif_load_raw(fname(cand(k)))); b(k) = prctile(J(:), 99.9); end
    ref_frame = cand(find(b >= 0.5 * max(b), 1, 'first'));
    fprintf('ref_frame not in the registry: using frame %d (first with p99.9 >= half its max over frames 1-1500).\n', ref_frame);
end
Iref = crop(lif_load_raw(fname(ref_frame)));
lev_ref = prctile(Iref(:), 99.9);
dx_mm = Tr.dx_mm;
y_cm  = Tr.y_mm / 10;      % cross-tank, 0 = plate centre
z_cm  = Tr.z_mm / 10;      % re-zeroed on the reference water level below
co.valid = Tr.footprint;   % zero-filled outside the camera's view
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
if isempty(fo.ylim_cm)     % zoom on the columns the sheet lights (the footprint is wider)
    lit = find(max(medfilt2(Iref, [5 5]), [], 1) >= su.min_counts);   % median: not hot pixels
    fo.ylim_cm = y_cm(lit([1 end])) + [-1 1];
    fprintf('Sheet lights y = %.1f .. %.1f cm (%.1f cm across)\n', y_cm(lit([1 end])), diff(y_cm(lit([1 end]))));
end

%% ---------------- Time ----------------
% Wind clock from the registry: hot-wire / wave-gauge file (104, 105, 121),
% or the assumed camera start (cam_delay_s, 10 s) with another run's wind
% as a labelled stand-in (119)
[camDelay, hw, wave] = lif_wind_clock(Rt, hw_avg_s, cal_dir);   % + wave gauge (own or stand-in)
t_wind = camDelay + (frames - 1) / fs;
t_what = 'since wind start';
if isempty(Rt.hw_file), t_what = 'since wind start (camera start assumed)'; end
show_hw = ~isempty(hw);
hw_file = Rt.hw_file;

%% ---------------- Figure + writer ----------------
if show_hw
    fo.hw = hw;
    fo.hw.xlim = hw_xlim;
    if ~isempty(wave), fo.wave = wave; end   % wave-gauge panel under the wind
end
fo.ct = lif_ct_profile(Rt.density_file);   % stratified runs: pre-run CT profile in the profile panel
if ~isempty(fo.ct), fprintf('CT profile %s: pycnocline %.1f..%.1f cm\n', fo.ct.file, fo.ct.pyc); end
fo.z99_floor_cm = NaN;                     % per frame from the gate
h = lif_trans_figure(y_cm, z_cm, fo);
G = lif_z99_gate(zg);

vid_file = lif_video_file(out_file);   % time-stamped name if the old video is open somewhere
v = VideoWriter(vid_file, 'MPEG-4');
v.FrameRate = video_fps;
v.Quality   = video_quality;
open(v);
cleanupObj = onCleanup(@() close(v));

nx = numel(y_cm);
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

    try
        [z_surf_all(i,:), ~] = lif_trans_surface(I, z_cm, dx_mm, su);
    catch
        % No dye in the sheet yet (119, 121: the dye reaches the sheet
        % ~20 s after the cameras start) -- hold the reference line; z99
        % comes out NaN from lif_trans_conc's min_peak gate
        z_surf_all(i,:) = z_ref;
    end
    R = lif_trans_conc(I, z_cm, z_surf_all(i,:), dx_mm, co);
    if prctile(I(:), 99.9) < sig_gate * lev_ref        % too little dye in the sheet
        R.z99 = NaN;  R.z99_y(:) = NaN;
    end
    z99_all(i)    = R.z99;
    z_pen_all(i,:) = R.z_pen;
    surf_all(i,:)  = R.surf_col;
    if isempty(z99_y_all)
        z99_y_all = nan(nFrames, numel(R.jc));
        y99_cm    = y_cm(R.jc);
        ref_col   = R.ref_col;
    end
    z99_y_all(i,:) = R.z99_y;
    [h.opts.z99_floor_cm, G] = lif_z99_gate(G, t_wind(i), z99_all(i));

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
fprintf('\nWrote %s  (%d frames, %.1f s)\n', vid_file, nFrames, toc);

%% ---------------- Save the z99 time series ----------------
% Frame ref_frame is the reference itself: C there is 1 by construction.
co = rmfield(co, intersect({'ref', 'dark', 'valid'}, fieldnames(co)));   % don't save images
co.ref_frame = ref_frame;
depth_all = z_surf_all - z_pen_all;        % [nFrames x nx] dye penetration, cm
[p, b] = fileparts(out_file);
ts_file = fullfile(p, [b '_z99.mat']);
rect = rmfield(Tr, {'U', 'V', 'footprint'});   % the flattening model, without the big maps
save(ts_file, 'frames', 't_wind', 't_what', 'camDelay', 'fs', 'y_cm', ...
     'z99_all', 'z99_y_all', 'y99_cm', 'z_surf_all', 'z_pen_all', 'depth_all', ...
     'surf_all', 'ref_col', 'rect', 'rc', ...
     'su', 'co', 'hw_file', 'run_num', 'Rt', 'ref_frame', 'G');   % G.t_on: when the z99 line came on
fprintf('Saved z99 time series to %s\n', ts_file);
