%% cisg_check_colocation.m
% IR <-> CISG colocation check.  Runs cisg_plate_calibrate, maps the IR-detected
% plate corners (cal_ir.image_points) into the CISG image through the shared
% tank plane and compares them with the corners detected in the CISG frame.
% Only the ~42 IR corners that fall in the CISG field of view are compared.
%
% Because both cals use the same lattice, the residual is the combined
% measurement error of the two corner sets (IR corners: 0.5 px ~ 0.6 mm), not
% a test of the CISG fit alone.  The block placement on the board (centre
% assumption, parity) is NOT tested here: it cannot be seen from corners.
%
% Writes calibration/cisg_colocation_check.png.
clearvars -except run_num; clc;                    % run_num = 95 (default), 96 or 101: uses that run's plate cal and IR snapshot cal
if ~exist('run_num', 'var'), run_num = 95; end
here = fileparts(mfilename('fullpath'));
if isempty(here), here = pwd; end
sg = fileparts(here);
addpath(here, sg, fullfile(sg, 'util'), fullfile(fileparts(sg), 'IR Camera', 'calibration'));

%% ---- 1. Calibrate ----
rd = cisg_run_defaults(run_num);
S = load(fullfile(here, sprintf('cisg_cal_CoreView_%d.mat', rd.plate)), 'cal');   % the saved cal (not recalibrated here)
cal_cisg = S.cal;
S = load(cal_cisg.ir_cal_path, 'cal');
cal_ir = S.cal;

%% ---- 2. IR corners -> CISG ----
[uv_map, xy_ir] = cisg_ir_to_cisg(cal_ir, cal_cisg, cal_ir.image_points);   % IR px -> plate mm -> tank mm -> CISG px
sz = cal_cisg.image_size;
in_fov = uv_map(:,1) >= 1 & uv_map(:,1) <= sz(2) & uv_map(:,2) >= 1 & uv_map(:,2) <= sz(1);

% pair with the CISG corners by ideal plate-frame position (shared lattice)
[found, idx] = ismembertol(cal_ir.world_points, cal_cisg.world_points_plate, 1e-6, 'ByRows', true);
pair = found & in_fov;
d_px = uv_map(pair,:) - cal_cisg.image_points(idx(pair),:);
xy_cisg = surf_px2world(cal_cisg, cal_cisg.image_points(idx(pair),:));
d_mm = xy_ir(pair,:) - xy_cisg;                   % IR-measured vs CISG-measured tank position

r_px = hypot(d_px(:,1), d_px(:,2));
r_mm = hypot(d_mm(:,1), d_mm(:,2));
fprintf('\n=== IR -> CISG corner colocation (%d paired of %d IR corners; %d in the CISG FOV) ===\n', ...
        nnz(pair), size(cal_ir.image_points,1), nnz(in_fov));
fprintf('  residual: mean %.2f px (%.3f mm), max %.2f px (%.3f mm), rms %.2f px\n', ...
        mean(r_px), mean(r_mm), max(r_px), max(r_mm), sqrt(mean(r_px.^2)));
fprintf('  mean offset (IR - CISG): dx %+.3f mm, dy %+.3f mm (%+.2f, %+.2f px)\n', ...
        mean(d_mm), mean(d_px));
fprintf('  CISG fit %.3f px rms; IR fit %.3f px rms; CISG %.5f mm/px\n', ...
        cal_cisg.rms_px, cal_ir.rms_px, cal_cisg.mm_per_px_center);

%% ---- 3. Overlay ----
I = double(rgb2gray(cisg_load_coreview(cal_cisg.frame_path, [3072 4100])));
fig = figure('Visible', 'off', 'Position', [100 100 1100 850]);
imagesc(I); colormap(gray); axis image; hold on;
h1 = plot(cal_cisg.image_points(:,1), cal_cisg.image_points(:,2), 'g+', 'MarkerSize', 10, 'LineWidth', 1.2);
h2 = plot(uv_map(in_fov,1), uv_map(in_fov,2), 'ro', 'MarkerSize', 9, 'LineWidth', 1);
% camera-aligned tank axes (cyan, x right / y up) and plate axes (yellow, rotated by theta), 100 mm long
p0 = cal_cisg.plate_centre_px;
px = surf_world2px(cal_cisg, [100 0]) - p0;
py = surf_world2px(cal_cisg, [0 100]) - p0;
Rp = cal_cisg.plate_to_tank;                      % plate axes expressed in tank coordinates
qx = surf_world2px(cal_cisg, 100 * Rp(:,1)') - p0;
qy = surf_world2px(cal_cisg, 100 * Rp(:,2)') - p0;
h3 = quiver(p0(1), p0(2), px(1), px(2), 0, 'c', 'LineWidth', 3, 'MaxHeadSize', 0.6);
quiver(p0(1), p0(2), py(1), py(2), 0, 'c', 'LineWidth', 3, 'MaxHeadSize', 0.6);
h4 = quiver(p0(1), p0(2), qx(1), qx(2), 0, 'y', 'LineWidth', 2, 'MaxHeadSize', 0.6);
quiver(p0(1), p0(2), qy(1), qy(2), 0, 'y', 'LineWidth', 2, 'MaxHeadSize', 0.6);
text(p0(1) + 1.15*px(1), p0(2) + 1.15*px(2), 'x tank', 'Color', 'c', 'FontSize', 14, 'FontWeight', 'bold');
text(p0(1) + 1.15*py(1), p0(2) + 1.15*py(2), 'y tank', 'Color', 'c', 'FontSize', 14, 'FontWeight', 'bold');
text(p0(1) + 1.15*qx(1), p0(2) + 1.15*qx(2) + 80, 'x plate', 'Color', 'y', 'FontSize', 14, 'FontWeight', 'bold');
text(p0(1) + 1.15*qy(1) + 80, p0(2) + 1.15*qy(2), 'y plate', 'Color', 'y', 'FontSize', 14, 'FontWeight', 'bold');
h5 = plot(p0(1), p0(2), 'cx', 'MarkerSize', 14, 'LineWidth', 2);
h6 = plot(cal_cisg.image_centre_px(1), cal_cisg.image_centre_px(2), 'm+', 'MarkerSize', 18, 'LineWidth', 2);
legend([h1 h2 h3 h4 h5 h6], {'CISG corners', 'IR corners mapped', 'tank axes (camera-aligned)', ...
       sprintf('plate axes (theta %+.2f deg)', cal_cisg.theta_deg), 'plate centre', 'image centre'}, ...
       'Color', [0.2 0.2 0.2], 'TextColor', 'w');
title(sprintf('CISG plate frame: CISG corners (g+) vs IR corners mapped (ro); mean %.2f px, max %.2f px', mean(r_px), max(r_px)));
xlabel('u (px)'); ylabel('v (px)');
png = fullfile(here, 'cisg_colocation_check.png');  if rd.plate ~= 97, png = fullfile(here, sprintf('cisg_colocation_check_%d.png', rd.plate)); end
exportgraphics(fig, png, 'Resolution', 100);
close(fig);
fprintf('Saved %s\n', png);
