%% Camera frame + hot-wire time series, with the frame's instant marked
% Top:    one frame from a CoreView camera (default: Water_PIV, the angled
%         view -- "PIV" is just the camera name, these are plain images).
% Bottom: hot-wire velocity from the matching run_experiment.m .mat file,
%         with a vertical bar at the time that frame was exposed.
%
% Timing: all cameras are triggered together at 50 Hz by run_experiment.m.
% Frame n is taken to be trigger n, so on the wind clock (0 = wind start,
% same as hwT_wind) it sits at
%     t_frame = (camStartElapsed - fanStartElapsed) + (n - 1 + frame_offset) / fs
% frame_offset absorbs any mismatch between the camera software's frame
% numbering and the first trigger (e.g. if recording was armed late).
clear; close all;
addpath(fullfile(fileparts(mfilename('fullpath')), 'util'));

%% ---------------- User settings ----------------
coreview_num = 95;
raw_root     = 'D:\HLAB_2026\Core';
cam          = 'Water_PIV';     % camera folder name (also in the file name)
frame        = 1875;            % same frame as the single Water_SURF image
frame_offset = 0;               % frames; see header note

hw_file = 'D:\HLAB_2026\hotwire\hotwire_20260923_113939.mat';

fs         = 50;                % Hz, camera trigger rate
downsample = 2;                 % plot every Nth pixel
pct_lims   = [1 99.5];          % colour limits from these percentiles
cmap       = 'gray';

%% ---------------- Load frame ----------------
cam_dir = fullfile(raw_root, sprintf('CoreView_%d', coreview_num), cam);
fname   = fullfile(cam_dir, sprintf('CoreView_%d_%s_%04d.raw', coreview_num, cam, frame));
if ~isfile(fname)
    lst = dir(fullfile(cam_dir, sprintf('CoreView_%d_%s_*.raw', coreview_num, cam)));
    n   = sort(cellfun(@(s) str2double(regexp(s, '_(\d+)\.raw$', 'tokens', 'once')), {lst.name}));
    error('Frame %d not found in %s (available: %d..%d, %d files).', ...
          frame, cam_dir, n(1), n(end), numel(n));
end
I = lif_load_raw(fname);
I = I(1:downsample:end, 1:downsample:end);

%% ---------------- Load hot-wire ----------------
H = load(hw_file);

if isfield(H, 'camStartElapsed') && ~isnan(H.camStartElapsed)
    camDelay = H.camStartElapsed - H.fanStartElapsed;
    camStop  = H.camStopElapsed  - H.fanStartElapsed;
else
    % Older files: fall back to the requested (not measured) delay
    camDelay = H.runConfig.DELAY_BEFORE_TRIG;
    camStop  = NaN;
    warning('No camStartElapsed in %s -- using DELAY_BEFORE_TRIG = %g s.', hw_file, camDelay);
end
t_frame = camDelay + (frame - 1 + frame_offset) / fs;
fprintf('Camera start %.3f s after wind start; frame %d at t = %.3f s\n', camDelay, frame, t_frame);

if t_frame < H.hwT_wind(1) || t_frame > H.hwT_wind(end)
    warning('t_frame = %.2f s is outside the hot-wire record (%.2f..%.2f s).', ...
            t_frame, H.hwT_wind(1), H.hwT_wind(end));
end

%% ---------------- Plot ----------------
figure('Name', sprintf('CoreView_%d %s frame %d + hot-wire', coreview_num, cam, frame), ...
       'Position', [60 60 1100 900], 'Color', 'w');
tl = tiledlayout(2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');

% Frame
ax1 = nexttile(tl);
[ny, nx] = size(I);
imagesc((1:nx)*downsample, (1:ny)*downsample, I);
axis image; colormap(ax1, cmap);
caxis(prctile(I(:), pct_lims));
c = colorbar; c.Label.String = 'counts';
xlabel('$x$ (px)', 'Interpreter', 'latex');
ylabel('$y$ (px)', 'Interpreter', 'latex');
title(sprintf('CoreView %d, %s, frame %d, $t = %.2f$ s', coreview_num, ...
      strrep(cam, '_', '\_'), frame, t_frame), 'Interpreter', 'latex');
set(ax1, 'FontSize', 13, 'FontName', 'times');

% Hot-wire
ax2 = nexttile(tl);
hold on;
leg = {};
if ~isempty(H.U)
    plot(H.hwT_wind, H.U, 'LineWidth', 0.75); leg{end+1} = '$U$ (hot-wire)';
end
if isfield(H, 'U_ref') && ~isempty(H.U_ref)
    plot(H.hwT_wind, H.U_ref, 'LineWidth', 0.75); leg{end+1} = '$U_{\rm ref}$';
end
yl = ylim;
if ~isnan(camStop)   % camera-recording window
    p = patch([camDelay camStop camStop camDelay], yl([1 1 2 2]), [0.85 0.85 0.85], ...
              'EdgeColor', 'none', 'HandleVisibility', 'off');
    uistack(p, 'bottom');
end
xline(t_frame, 'r-', 'LineWidth', 2);
leg{end+1} = sprintf('frame %d', frame);
ylim(yl); hold off; box on; grid on;
xlabel('Time since wind start (s)');
ylabel('Velocity (m/s)');
legend(leg, 'Interpreter', 'latex', 'Location', 'best');
[~, hw_name] = fileparts(hw_file);
title(strrep(hw_name, '_', '\_'));
set(ax2, 'FontSize', 13, 'FontName', 'times');
