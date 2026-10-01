%% distortion_correction_movie.m
% Side-by-side movie of an IR recording: raw camera frame (left) and the
% distortion-corrected frame on the water surface in tank cm (right), using
% the single-view board calibration (calibration/calibrate_surface_single.m)
% instead of the side-wall polyfits of distortion_correction.m.
%
% Default: Rec-000094 with the coarse (40 mm, full-size) board snapshot
% Snap-000096. The camera must not have moved between the two files.
%
% Cells:
%   1. config
%   2. surface calibration from the snapshot (or load the saved one)
%   3. distortion-correction map + wall cross-check on the recording itself
%   4. pass 1: time-mean + colour limits from a strided subset
%   5. pass 2: write the movie

clear; clc;

%% ---- 1. Config ----
here = fileparts(mfilename('fullpath'));
if isempty(here), here = pwd; end
addpath(here, fullfile(here, 'calibration'));
addpath('C:\Program Files\FLIR Systems\sdks\file\bin\Release');

cfg.rec_path   = 'D:\HLAB_2026\IR_camera\Rec-000094.ats';
cfg.snap_path  = 'D:\HLAB_2026\IR_camera\Snap-000096.ats';
cfg.cal_path   = fullfile(here, 'calibration', 'surface_cal_Snap-000096.mat');

cfg.recalibrate = true;
cfg.cal_opts = struct('repick', {{'pattern bottom edge', 'left wall'}});

cfg.out_path   = 'D:\HLAB_2026\IR_camera\Processed\Rec-000094_distortion_correction.mp4';
cfg.dx_mm      = 1;               % corrected grid spacing
cfg.wall_y_mm  = 250;             % clip cross-wind extent at the side walls
cfg.x_lim_mm   = [];              % along-wind range [min max]; [] = board corners +/- margin
cfg.x_margin_mm = 60;             % extrapolation beyond the outermost corners (1.5 squares)
cfg.frames     = [];              % [] = all, else e.g. 1:2000
cfg.mean_stride = 10;             % every Nth frame for the colour limits (and time-mean)
cfg.ref_frames = 1:10;            % reference subtracted in 'anomaly' mode:
                                  %   1:10 = initial quiescent frames, as in
                                  %          distortion_correction.m (img_avg)
                                  %   []   = mean over the whole record (strided)
cfg.mode       = 'anomaly';       % 'anomaly' (T - time mean) or 'temperature'
cfg.remove_frame_offset = false;   % subtract each frame's spatial median (anomaly mode):
                                  % removes uniform drift (surface cooling, camera
                                  % drift, NUC jumps) so only spatial structure shows
cfg.clim       = [];              % [] = 1st/99th percentile of the subset
cfg.clim_symmetric = false;       % false: limits follow the data (e.g. -0.8..0.1);
                                  % the colormap is rebuilt so white stays at 0
cfg.mask_color = [0.6 0.6 0.6];   % gray outside the side walls / outside the image
cfg.cmap       = hot(256);        % used as-is over cfg.clim, e.g. clim = [-0.8 0]
cfg.zero_white = false;           % true with cmap = diverging_cmap(256): white kept at 0
                                  % symmetric (Langmuir cells). 'hot' suits
                                  % absolute T, not a signed anomaly.
cfg.fps_out    = [];              % [] = recording frame rate
cfg.quality    = 95;

%% ---- 2. Surface calibration ----
if cfg.recalibrate || ~isfile(cfg.cal_path)
    cfg.cal_opts.out_path = cfg.cal_path;
    % Walls are clicked on the recording (water/wall contrast shows the
    % waterline); pattern edges on the snapshot. Mean of a few frames.
    if ~isfield(cfg.cal_opts, 'wall_image')
        vr = FlirMovieReader(cfg.rec_path);
        vr.unit = 'temperatureFactory';
        acc = 0;
        for i = 1:30, acc = acc + double(step(vr)); end
        cfg.cal_opts.wall_image = acc / 30;
        clear vr
    end
    cal = calibrate_surface_single(cfg.snap_path, cfg.cal_opts);
else
    S = load(cfg.cal_path, 'cal');
    cal = S.cal;
    fprintf('Loaded %s (rms %.2f px)\n', cfg.cal_path, cal.rms_px);
end

%% ---- 3. Distortion-correction map ----
v = FlirMovieReader(cfg.rec_path);
v.unit = 'temperatureFactory';
nF  = double(v.sourceInfo.presetInfo(1).numFrames);
fps = double(v.sourceInfo.presetInfo(1).frameRate);
img_size = double([v.sourceInfo.imageHeight v.sourceInfo.imageWidth]);   % SDK returns ints
assert(isequal(img_size, cal.image_size), ...
    'Recording is %dx%d but calibration is %dx%d.', img_size, cal.image_size);
if isempty(cfg.frames), cfg.frames = 1:nF; end
if isempty(cfg.fps_out), cfg.fps_out = fps; end
fprintf('%s: %d frames at %.1f Hz\n', cfg.rec_path, nF, fps);

bb = cal.bbox_mm;
% Along-wind extent from the detected board corners plus a margin. Pixels
% outside the board are unconstrained: near the top of the image the
% homography approaches its horizon and beyond the last corners the radial
% polynomial extrapolates, so mapping image-edge pixels gave x ranges of
% hundreds of metres.
if isempty(cfg.x_lim_mm)
    xw = cal.world_points(:,1);
    cfg.x_lim_mm = [min(xw) max(xw)] + cfg.x_margin_mm * [-1 1];
end
xg = ceil(cfg.x_lim_mm(1)) : cfg.dx_mm : floor(cfg.x_lim_mm(2));         % along-wind
yg = max(ceil(bb(3)), -cfg.wall_y_mm) : cfg.dx_mm : min(floor(bb(4)), cfg.wall_y_mm);  % cross-wind
[Yg, Xg] = meshgrid(yg, xg);
p = surf_world2px(cal, [Xg(:) Yg(:)]);
U = reshape(p(:,1), size(Xg));
V = reshape(p(:,2), size(Xg));
outside = U < 1 | U > img_size(2) | V < 1 | V > img_size(1);
fprintf('Corrected grid %d x %d at %.1f mm (x %.0f..%.0f, y %.0f..%.0f mm)\n', ...
    numel(xg), numel(yg), cfg.dx_mm, xg([1 end]), yg([1 end]));

F = griddedInterpolant({1:img_size(1), 1:img_size(2)}, zeros(img_size), 'linear', 'none');

% Cross-check against the side walls: the predicted wall lines (cyan)
% should sit on the walls in the recording itself. This is the
% independent check the old wall-polyfit method provided.
reset(v);
f1 = double(step(v));
figure('Name', 'Wall check on recording', 'Position', [100 100 1300 550]);
subplot(1,2,1);
imagesc(f1); colormap(gca, 'hot'); axis image; hold on; colorbar;   % absolute T
xs = linspace(cfg.x_lim_mm(1), cfg.x_lim_mm(2), 300)';
for yw = cfg.wall_y_mm * [-1 1]
    pw = surf_world2px(cal, [xs, yw*ones(size(xs))]);
    plot(pw(:,1), pw(:,2), 'c-', 'LineWidth', 1.5);
end
title('Frame 1 with predicted walls y = \pm250 mm (cyan)');
subplot(1,2,2);
imagesc(yg, xg, correct_frame(F, f1, V, U, outside)); colormap(gca, 'hot');
axis image; colorbar;
xlabel('y cross-wind (mm)'); ylabel('x along-wind (mm)');
title('Frame 1 distortion corrected');

%% ---- 4. Pass 1: time-mean and colour limits ----
idx_mean = cfg.frames(1:cfg.mean_stride:end);
if isempty(cfg.ref_frames), idx_ref = idx_mean; else, idx_ref = cfg.ref_frames; end
acc = zeros(img_size);
sub = zeros([numel(xg) numel(yg) numel(idx_mean)], 'single');
reset(v);
n = 0; n_ref = 0;
for fi = 1:max([idx_mean idx_ref])
    fr = step(v);
    if any(fi == idx_ref)
        n_ref = n_ref + 1;
        acc = acc + double(fr);
    end
    if any(fi == idx_mean)
        n = n + 1;
        sub(:,:,n) = correct_frame(F, double(fr), V, U, outside);
    end
    if mod(fi, 500) == 0, fprintf('  pass 1: frame %d / %d\n', fi, max(idx_mean)); end
end
T_mean = acc / n_ref;                  % reference field (initial frames or record mean)
fprintf('Reference: mean of %d frames (%d..%d)\n', n_ref, min(idx_ref), max(idx_ref));
R_mean = correct_frame(F, T_mean, V, U, outside);

if strcmpi(cfg.mode, 'anomaly')
    sub = sub - single(R_mean);
    if cfg.remove_frame_offset
        for k = 1:n
            sk = sub(:,:,k);
            sub(:,:,k) = sk - median(sk(isfinite(sk)));
        end
    end
end
if isempty(cfg.clim)
    s = sort(sub(isfinite(sub)));
    cfg.clim = double(s(round([0.01 0.99] * numel(s))))';
    if strcmpi(cfg.mode, 'anomaly') && cfg.clim_symmetric
        cfg.clim = max(abs(cfg.clim)) * [-1 1];
    end
end
% Colormap matched to the limits: white at 0, equal |dT| -> equal colour
% strength, whether or not the limits are symmetric
% cfg.cmap is used as-is over cfg.clim, unless cfg.zero_white: then a
% diverging map is resampled so white stays at 0 for skewed limits.
if ischar(cfg.cmap), cmap_movie = feval(cfg.cmap, 256); else, cmap_movie = cfg.cmap; end
if isfield(cfg, 'zero_white') && cfg.zero_white
    cmap_movie = zero_centred_cmap(cmap_movie, cfg.clim);
end
fprintf('Colour limits [%.3f %.3f] (%s)\n', cfg.clim, cfg.mode);
clear sub

%% ---- 5. Pass 2: write movie ----
out_dir = fileparts(cfg.out_path);
if ~isfolder(out_dir), mkdir(out_dir); end
vw = VideoWriter(cfg.out_path, 'MPEG-4');
vw.FrameRate = cfg.fps_out;
vw.Quality   = cfg.quality;
open(vw);

fig = figure('Color', 'w', 'Position', [50 100 1400 640], 'Visible', 'off');
[~, rec_name] = fileparts(cfg.rec_path);
ht = sgtitle(fig, '', 'Interpreter', 'none', 'FontSize', 12);

% Pixels of the raw frame outside the side walls (|y| > wall_y_mm) are
% grayed out, matching the clipped extent of the corrected panel
[uu, vv] = meshgrid(1:img_size(2), 1:img_size(1));
xy_px = surf_px2world(cal, [uu(:) vv(:)]);
raw_mask = reshape(~all(isfinite(xy_px), 2) | abs(xy_px(:,2)) > cfg.wall_y_mm, img_size);

% Left: raw camera frame, pixel coordinates
ax1 = subplot(1, 2, 1, 'Parent', fig);
h1  = imagesc(ax1, nan(img_size));
set(h1, 'AlphaData', ~raw_mask);
set(ax1, 'Color', cfg.mask_color);
axis(ax1, 'image'); colormap(ax1, cmap_movie); caxis(ax1, cfg.clim);
xlabel(ax1, 'u (px)'); ylabel(ax1, 'v (px)');
title(ax1, 'Raw');

% Right: distortion corrected, tank coordinates
ax2 = subplot(1, 2, 2, 'Parent', fig);
h2  = imagesc(ax2, yg/10, xg/10, nan(numel(xg), numel(yg)));
set(h2, 'AlphaData', ~outside);
set(ax2, 'Color', cfg.mask_color);
axis(ax2, 'image'); colormap(ax2, cmap_movie); caxis(ax2, cfg.clim);
xlabel(ax2, 'cross-wind y (cm)'); ylabel(ax2, 'along-wind x (cm)');
title(ax2, 'Distortion corrected');

cb = colorbar(ax2);
if strcmpi(cfg.mode, 'anomaly')
    cb.Label.String = 'T - \langle T \rangle (\circC)';
else
    cb.Label.String = 'T (\circC)';
end

reset(v);
last = max(cfg.frames);
offset = nan(1, last);                 % spatial median of each frame's anomaly
for fi = 1:last
    fr = double(step(v));
    if ~any(fi == cfg.frames), continue; end
    R = correct_frame(F, fr, V, U, outside);
    if strcmpi(cfg.mode, 'anomaly')
        fr = fr - T_mean;
        R  = R - R_mean;
        offset(fi) = median(R(isfinite(R)));
        if cfg.remove_frame_offset      % same scalar on both panels
            fr = fr - offset(fi);
            R  = R  - offset(fi);
        end
    end
    set(h1, 'CData', fr);
    set(h2, 'CData', R);
    set(ht, 'String', sprintf('%s   frame %d   t = %.3f s', rec_name, fi, (fi - 1) / fps));
    writeVideo(vw, getframe(fig));
    if mod(fi, 200) == 0, fprintf('  pass 2: frame %d / %d\n', fi, last); end
end
close(vw); close(fig);
fprintf('Wrote %s\n', cfg.out_path);

% Whole-frame temperature offset vs time: a smooth trend is real surface
% cooling/warming or camera body drift; sudden steps are NUC events.
if strcmpi(cfg.mode, 'anomaly')
    t = ((1:last) - 1) / fps;
    figure('Name', 'Frame offset');
    plot(t, offset, 'k.-');
    xlabel('t (s)'); ylabel('median(T - \langle T \rangle) over surface (\circC)');
    grid on;
    if cfg.remove_frame_offset, s_ = 'removed from movie'; else, s_ = 'kept in movie'; end
    title(sprintf('%s whole-frame offset (%s)', rec_name, s_), 'Interpreter', 'none');
end

%% ---- local functions ----
function R = correct_frame(F, fr, V, U, outside)
F.Values = fr;
R = F(V, U);
R(outside) = NaN;
end

function cmap = diverging_cmap(n)
% Blue - white - red, linear in each half (no toolbox needed)
h = floor(n/2);
b = [0.02 0.19 0.38];  r = [0.40 0.00 0.05];  w = [1 1 1];
mid_b = [0.26 0.58 0.76]; mid_r = [0.84 0.38 0.30];
lo = interp1([0 0.5 1], [b; mid_b; w], linspace(0, 1, h));
hi = interp1([0 0.5 1], [w; mid_r; r], linspace(0, 1, n - h));
cmap = [lo; hi];
end

function cmap = zero_centred_cmap(base, clim, n)
% Resample a symmetric diverging map so that, over [clim(1) clim(2)], white
% sits at 0 and colour strength scales with |value| / max(|clim|). With
% skewed limits (e.g. -0.8..0.1) the map is mostly blue with a little red,
% and a -0.1 and a +0.1 pixel still look equally strong.
if nargin < 3, n = 256; end
t = linspace(clim(1), clim(2), n)' / max(abs(clim));      % in [-1, 1]
s = linspace(-1, 1, size(base, 1))';
cmap = interp1(s, base, t);
end
