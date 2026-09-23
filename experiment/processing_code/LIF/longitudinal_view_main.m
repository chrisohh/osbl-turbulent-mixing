%% Longitudinal view - two-camera snapshot viewer (dye visualisation check)
% Loads ONE frame from each of the two MCL cameras in a CoreView run and
% displays them side by side so the contrast/contour levels can be tuned
% until the dye is visible.  No spatial calibration is applied: axes are in
% pixels.
clear; close all;
addpath(fullfile(fileparts(mfilename('fullpath')), 'util'));

%% ---------------- User settings ----------------
coreview_num = 89;
raw_root     = 'D:\HLAB_2026\Core';

cams = { 'Flare 12M125 MCL C0100800', ...   % camera 1 (longitudinal/side view)
         'Flare 12M125 MCL C0100801' };     % camera 2

snapshot_frame = 1800;      % <-- the single frame to plot
ref_frame      = 1;         % background frame (dye-free) for normalisation

fs = 50;                    % Hz (only used for the time label)

% Display mode:
%   'raw'   - plot raw counts
%   'norm'  - plot I/I_ref (transmission).  Dye absorbs -> values < 1.
%   'od'    - plot -log10(I/I_ref) (optical depth).  Dye -> values > 0.
display_mode = 'od';



%% ---------------- Load ----------------
dt = 1/fs;
nCam = numel(cams);
I    = cell(1, nCam);
Iref = cell(1, nCam);

for k = 1:nCam
    cam_dir = fullfile(raw_root, sprintf('CoreView_%d', coreview_num), cams{k});
    fname = @(n) fullfile(cam_dir, sprintf('CoreView_%d_%s_%04d.raw', ...
                                           coreview_num, cams{k}, n));

    fprintf('Camera %d: %s\n', k, cams{k});
    fprintf('  snapshot: %s\n', fname(snapshot_frame));
    I{k} = lif_load_raw(fname(snapshot_frame));

    if ~strcmp(display_mode, 'raw')
        Iref{k} = lif_load_raw(fname(ref_frame));
    end

    fprintf('  counts: min=%g  median=%g  p99.9=%g  max=%g\n', ...
            min(I{k}(:)), median(I{k}(:)), prctile(I{k}(:), 99.9), max(I{k}(:)));
end

%% ---------------- Build the displayed field ----------------
F = cell(1, nCam);
for k = 1:nCam
    switch display_mode
        case 'raw'
            F{k} = I{k};
            cbl  = 'counts';
        case 'norm'
            F{k} = I{k} ./ max(Iref{k}, 1);
            cbl  = '$I/I_{\rm ref}$';
        case 'od'
            F{k} = -log10(max(I{k}, 1) ./ max(Iref{k}, 1));
            cbl  = '$-\log_{10}(I/I_{\rm ref})$';
        otherwise
            error('Unknown display_mode: %s', display_mode);
    end
    F{k} = F{k}(1:downsample:end, 1:downsample:end);
end

%% ---------------- Plot ----------------
% Colour limits.  Leave empty to auto-set from percentiles of the frame.
clim_raw  = [];             % e.g. [0 700]
clim_norm = [0.7 1.05];
clim_od   = [0 0.05];
pct_lims  = [1 99.5];       % percentiles used when the clim above is empty

% Contour overlay (set to [] to switch off).  Levels are in the units of
% whatever display_mode produces.
contour_levels = [];%[0.02 0.05 0.10 0.20];   % optical-depth contours
contour_smooth = 0;%9;         % box-filter width in px before contouring (0 = none)

downsample = 2;             % plot every Nth pixel (2 keeps 4096x3072 responsive)
cmap       = 'bone';      % try 'gray', 'parula', 'turbo', 'hot'
y_dir      = 'reverse';     % 'reverse' = image as the camera sees it,
                            % 'normal'  = flip so row 1 is at the bottom

figure('Name', sprintf('CoreView_%d frame %d', coreview_num, snapshot_frame), ...
       'Position', [60 60 1500 650], 'Color', 'w');
tl = tiledlayout(1, nCam, 'TileSpacing', 'compact', 'Padding', 'compact');

k=1;
    nexttile;
    [ny, nx] = size(F{k});
    xpx = (1:nx) * downsample;
    ypx = (1:ny) * downsample;

    imagesc(xpx, ypx, F{k},100);
    axis image; set(gca, 'YDir', y_dir);
    colormap(gca, cmap);

    switch display_mode
        case 'raw',  cl = clim_raw;
        case 'norm', cl = clim_norm;
        case 'od',   cl = clim_od;
    end
    if isempty(cl)
        cl = prctile(F{k}(:), pct_lims);
    end
    caxis(cl);

    if ~isempty(contour_levels)
        G = F{k};
        if contour_smooth > 1
            G = conv2(G, ones(contour_smooth)/contour_smooth^2, 'same');
        end
        hold on;
        contour(xpx, ypx, G, contour_levels, 'k', 'LineWidth', 0.75);
        hold off;
    end

    c = colorbar;
    c.Label.Interpreter = 'latex';
    c.Label.String      = cbl;
    c.Label.FontSize    = 13;

    title(sprintf('cam %d (%s)', k, cams{k}(end-8:end)), 'Interpreter', 'none');
    xlabel('$x$ (px)', 'Interpreter', 'latex');
    ylabel('$z$ (px)', 'Interpreter', 'latex');
    set(gca, 'FontSize', 13, 'FontName', 'times');
%%
k=2;

    nexttile;
    [ny, nx] = size(F{k});
    xpx = (1:nx) * downsample;
    ypx = (1:ny) * downsample;

    imagesc(xpx, ypx, F{k});
    axis image; set(gca, 'YDir', y_dir);
    colormap(gca, cmap);

    switch display_mode
        case 'raw',  cl = clim_raw;
        case 'norm', cl = clim_norm;
        case 'od',   cl = [0 0.3]
    end
    if isempty(cl)
        cl = prctile(F{k}(:), pct_lims);
    end
    caxis(cl);

    if ~isempty(contour_levels)
        G = F{k};
        if contour_smooth > 1
            G = conv2(G, ones(contour_smooth)/contour_smooth^2, 'same');
        end
        hold on;
        contour(xpx, ypx, G, contour_levels, 'k', 'LineWidth', 0.75);
        hold off;
    end

    c = colorbar;
    c.Label.Interpreter = 'latex';
    c.Label.String      = cbl;
    c.Label.FontSize    = 13;

    title(sprintf('cam %d (%s)', k, cams{k}(end-8:end)), 'Interpreter', 'none');
    xlabel('$x$ (px)', 'Interpreter', 'latex');
    ylabel('$z$ (px)', 'Interpreter', 'latex');
    set(gca, 'FontSize', 13, 'FontName', 'times');

title(tl, sprintf('CoreView %d, frame %d, $t = %.2f$ s, mode: %s', ...
      coreview_num, snapshot_frame, (snapshot_frame-1)*dt, display_mode), ...
      'Interpreter', 'latex', 'FontSize', 15);

%% ---------------- Tuning help ----------------
% Intensity histogram of the displayed field - use it to pick clim/contours.
figure('Name', 'Displayed-field histogram', 'Position', [80 80 900 350], 'Color', 'w');
for k = 1:nCam
    subplot(1, nCam, k);
    histogram(F{k}(:), 200, 'EdgeColor', 'none');
    set(gca, 'YScale', 'log'); grid on;
    xlabel(cbl, 'Interpreter', 'latex'); ylabel('count');
    title(sprintf('cam %d', k));
end

fprintf('\nDone. To hunt for the dye: change snapshot_frame, then tighten\n');
fprintf('clim_%s / contour_levels using the histogram above.\n', display_mode);
