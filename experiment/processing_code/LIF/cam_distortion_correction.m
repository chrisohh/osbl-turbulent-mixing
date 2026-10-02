%% cam_distortion_correction.m
% Raw | plate-corrected view for the CoreView cameras, following
% IR Camera/distortion_correction_movie.m.  Calibration from ONE plate
% image per camera (util/cam_plate_calibrate.m: homography; the lenses are
% not fisheye, so radial distortion is off by default and only reported).
%
% Cells:
%   1. config (pick the camera)
%   2. calibrate from the plate image + check figure:
%        plate raw with detected corners, reprojection residual x20 and the
%        back-projected world grid | plate corrected onto the plane
%   3. a data frame: raw | corrected in cm (z = 0 at the frame-1 water
%      level for the side views)
%   4. optional movie of raw | corrected over frames
%
% Cameras (CoreView 95):
%   Water_SURF  side view,    plate Core 94, small plate (11 mm), 2.5 cm toward camera
%   Water_PIV   angled view,  plate Core 91, small plate (11 mm)
%   nadir       CISG colour,  plate Core 97, big plate (40 mm), colocated with IR snap_96
clear; clc;
here = fileparts(mfilename('fullpath'));
if isempty(here), here = pwd; end
addpath(fullfile(here, 'util'), fullfile(here, '..', 'Slope_Gauge'), ...   % cisg_load_coreview
        fullfile(here, '..', 'IR Camera', 'calibration'));             % surf_world2px / surf_px2world

%% ---- 1. Config ----
cam_name = 'Water_SURF';          % 'Water_SURF' | 'Water_PIV' | 'nadir'

raw_root = 'D:\HLAB_2026\Core';
run_num  = 95;
frame    = 1875;                  % data frame for cell 3
datum_frame = 1;                  % dye-free first frame: z = 0 at its water level (side views)

cams.Water_SURF = struct('folder', 'Water_Surf', 'loader', 'mono', ...
    'plate', fullfile(raw_root, 'CoreView_94', 'Water_Surf', 'CoreView_94_Water_Surf_01.raw'), ...
    'square_mm', 11, 'offset_mm', 25, 'y_up', true, 'side', true, ...
    'xlab', 'x along tank', 'ylab', 'z');
cams.Water_PIV  = struct('folder', 'Water_PIV', 'loader', 'mono', ...
    'plate', fullfile(raw_root, 'CoreView_91', 'Water_PIV', 'CoreView_91_Water_PIV_01.raw'), ...
    'square_mm', 11, 'offset_mm', 0, 'y_up', true, 'side', false, ...
    'xlab', 'x (plate plane)', 'ylab', 'y (plate plane)');
cams.nadir      = struct('folder', 'Flare 12M125 CCL C1112A00', 'loader', 'color', ...
    'plate', fullfile(raw_root, 'CoreView_97', 'Flare 12M125 CCL C1112A00', ...
                      'CoreView_97_Flare 12M125 CCL C1112A00_001.raw'), ...
    'square_mm', 40, 'offset_mm', 0, 'y_up', false, 'side', false, ...
    'xlab', 'x along wind', 'ylab', 'y cross wind');
cam = cams.(cam_name);

cal_opts = struct('square_mm', cam.square_mm, 'offset_mm', cam.offset_mm, ...
                  'n_medium', 1.33, 'y_up', cam.y_up, 'n_radial', 0, ...
                  'lens_f_mm', 35, 'pixel_um', 5.5);     % lens_f_mm: CHECK per camera
cal_path = fullfile(here, 'calibration', sprintf('cal_%s.mat', cam_name));
recalibrate = true;

dx_mm      = 0.5;                 % corrected grid spacing
extent     = 'image';             % 'image' = whole image footprint, 'board' = corners + margin
grid_mm    = 20;                  % world grid overlaid on the raw plate (check)
cmap       = 'bone';
clims      = [];                  % counts; [] = 1st/99.5th percentile

make_movie  = false;              % cell 4
movie_frames = 1:5:1875;
movie_path  = fullfile(raw_root, sprintf('CoreView_%d', run_num), ...
                       sprintf('CoreView_%d_%s_corrected.mp4', run_num, cam_name));
movie_fps   = 10;

load_img = @(f) load_frame(f, cam.loader);
fname = @(n) fullfile(raw_root, sprintf('CoreView_%d', run_num), cam.folder, ...
                      sprintf('CoreView_%d_%s_%04d.raw', run_num, cam.folder, n));

%% ---- 2. Calibrate from the plate ----
Pl = load_img(cam.plate);
if recalibrate || ~isfile(cal_path)
    cal = cam_plate_calibrate(Pl, cal_opts);
    if ~isfolder(fileparts(cal_path)), mkdir(fileparts(cal_path)); end
    save(cal_path, 'cal');
    fprintf('Saved %s\n', cal_path);
else
    S = load(cal_path, 'cal');  cal = S.cal;
    fprintf('Loaded %s (rms %.2f px)\n', cal_path, cal.rms_px);
end

% Corrected grid (row 1 = top = largest y)
img_size = cal.image_size;
switch extent
    case 'image'
        [uu, vv] = meshgrid(linspace(1, img_size(2), 60), linspace(1, img_size(1), 60));
        xy = surf_px2world(cal, [uu(:) vv(:)]);
        bb = [min(xy(:,1)) max(xy(:,1)) min(xy(:,2)) max(xy(:,2))];
    case 'board'
        bb = cal.bbox_mm;
end
xg = ceil(bb(1)) : dx_mm : floor(bb(2));
yg = floor(bb(4)) : -dx_mm : ceil(bb(3));
[Xg, Yg] = meshgrid(xg, yg);
p = surf_world2px(cal, [Xg(:) Yg(:)]);
U = reshape(p(:,1), size(Xg));  V = reshape(p(:,2), size(Xg));
outside = U < 1 | U > img_size(2) | V < 1 | V > img_size(1);
correct = @(img) mask_out(interp2(double(img), U, V, 'linear', NaN), outside);
fprintf('Corrected grid %d x %d at %.2f mm (x %.0f..%.0f, y %.0f..%.0f mm)\n', ...
        numel(yg), numel(xg), dx_mm, xg([1 end]), yg([end 1]));

% Check figure
figure('Name', ['Plate correction ' cam_name], 'Color', 'w', 'Position', [60 80 1600 650]);
subplot(1, 2, 1);
imagesc(Pl); colormap(gca, gray); axis image; hold on;
pts = cal.image_points;
pr  = surf_world2px(cal, cal.world_points);
plot(pts(:,1), pts(:,2), 'g+', 'MarkerSize', 6);
quiver(pts(:,1), pts(:,2), 20*(pr(:,1)-pts(:,1)), 20*(pr(:,2)-pts(:,2)), 0, 'r');
overlay_grid(cal, bb, grid_mm);
title(sprintf('Plate raw: corners (+), residual x20 (red), %d mm grid (cyan); rms %.2f px', ...
      grid_mm, cal.rms_px));
xlabel('u (px)'); ylabel('v (px)');
subplot(1, 2, 2);
imagesc(xg/10, yg/10, correct(Pl)); colormap(gca, gray); axis image; set(gca, 'YDir', 'normal');
hold on; plot(cal.world_points(:,1)/10, cal.world_points(:,2)/10, 'g+', 'MarkerSize', 6);
xlabel([cam.xlab ' (cm)']); ylabel([cam.ylab ' (cm, board centre)']);
title(sprintf('Plate corrected (%.2f mm grid), %.4f mm/px at board', dx_mm, cal.mm_per_px_center));

%% ---- 3. Data frame: raw | corrected ----
I  = load_img(fname(frame));
Ic = correct(I);
y0 = 0;  ylab = [cam.ylab ' (cm, board centre)'];
if cam.side && ~isempty(datum_frame)
    % z = 0 at the datum frame's water level: darkest row of the corrected
    % datum frame (as lif_dye_front), averaged along x
    Bc = correct(load_img(fname(datum_frame)));
    dy = struct('surf_window_cm', [4 10], 'smooth_mm', 1.5);
    [~, zs] = lif_dye_front(fillmissing(Bc, 'nearest', 2), yg/10, dx_mm, dy);
    y0 = 10 * mean(zs, 'omitnan');
    ylab = [cam.ylab sprintf(' (cm, 0 = frame-%d water level)', datum_frame)];
    fprintf('Datum: frame-%d water level at %.2f cm above the board centre\n', datum_frame, y0/10);
end
if isempty(clims), clims = prctile(I(:), [1 99.5]); end

figure('Name', sprintf('%s frame %d raw | corrected', cam_name, frame), 'Color', 'w', ...
       'Position', [60 80 1700 650]);
subplot(1, 2, 1);
imagesc(I); colormap(gca, cmap); caxis(clims); axis image; colorbar;
xlabel('u (px)'); ylabel('v (px)'); title(sprintf('Raw, frame %d', frame));
subplot(1, 2, 2);
h = imagesc(xg/10, (yg - y0)/10, Ic); set(h, 'AlphaData', ~isnan(Ic));
colormap(gca, cmap); caxis(clims); axis image; set(gca, 'YDir', 'normal'); colorbar;
xlabel([cam.xlab ' (cm)']); ylabel(ylab);
title('Plate corrected');

%% ---- 4. Movie: raw | corrected ----
if make_movie
    vw = VideoWriter(movie_path, 'MPEG-4');
    vw.FrameRate = movie_fps;  vw.Quality = 95;
    open(vw);
    fig = figure('Color', 'w', 'Position', [50 100 1700 650], 'Visible', 'off');
    ax1 = subplot(1, 2, 1, 'Parent', fig);
    h1 = imagesc(ax1, nan(img_size)); axis(ax1, 'image'); colormap(ax1, cmap); caxis(ax1, clims);
    xlabel(ax1, 'u (px)'); ylabel(ax1, 'v (px)'); title(ax1, 'Raw');
    ax2 = subplot(1, 2, 2, 'Parent', fig);
    h2 = imagesc(ax2, xg/10, (yg - y0)/10, nan(size(Xg))); set(h2, 'AlphaData', ~outside);
    axis(ax2, 'image'); set(ax2, 'YDir', 'normal'); colormap(ax2, cmap); caxis(ax2, clims);
    xlabel(ax2, [cam.xlab ' (cm)']); ylabel(ax2, ylab); title(ax2, 'Plate corrected');
    ht = sgtitle(fig, '', 'Interpreter', 'none');
    for n = movie_frames
        f = fname(n);
        d = dir(f);
        if isempty(d) || d.bytes < 1e6, continue; end         % missing / still copying
        I = load_img(f);
        set(h1, 'CData', I);  set(h2, 'CData', correct(I));
        set(ht, 'String', sprintf('CoreView %d %s  frame %d', run_num, cam_name, n));
        rgb = print(fig, '-RGBImage', '-r0');
        writeVideo(vw, rgb(1:2*floor(end/2), 1:2*floor(end/2), :));
    end
    close(vw); close(fig);
    fprintf('Wrote %s\n', movie_path);
end

%% ---- local functions ----
function I = load_frame(f, loader)
    switch loader
        case 'mono',  I = lif_load_raw(f);
        case 'color', I = double(rgb2gray(cisg_load_coreview(f)));
    end
end

function R = mask_out(R, outside)
    R(outside) = NaN;
end

function overlay_grid(cal, bb, spacing)
    gx = ceil(bb(1)/spacing)*spacing : spacing : floor(bb(2)/spacing)*spacing;
    gy = ceil(bb(3)/spacing)*spacing : spacing : floor(bb(4)/spacing)*spacing;
    xs = linspace(bb(1), bb(2), 200)';  ys = linspace(bb(3), bb(4), 200)';
    for y = gy
        p = surf_world2px(cal, [xs, y*ones(size(xs))]);
        plot(p(:,1), p(:,2), 'c-', 'LineWidth', 0.4);
    end
    for x = gx
        p = surf_world2px(cal, [x*ones(size(ys)), ys]);
        plot(p(:,1), p(:,2), 'c-', 'LineWidth', 0.4);
    end
end
