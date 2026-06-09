% plot_check2_not_swaying.m
%
% "Is it not swaying?" — Check 2 (visual scale-separation sanity check).
%
% Loads intrp_u_res from PIVMat_TURB (built by run_decomposition_loop.m),
% computes the horizontal-mean <u>(z0, t) at 3-4 depths, then applies a
% sliding-window average over a swept range of window lengths Delta_t.
% All curves for all window lengths land on one panel per depth.
%
% Interpretation:
%   Short Delta_t curves  → wiggly (wave-period oscillations contaminate mean)
%   Long  Delta_t curves  → flat   (plateau = valid scale-separation window)
%
% Optional: set make_video = true to export an MP4 that sweeps through
% window lengths so the curve visibly settles — useful for presentations.

clear; clc;

%% =========================================================================
%% PARAMETERS  (keep in sync with run_decomposition_loop.m)
%% =========================================================================
LONG = 'D:\DelawareDataBackup\Longitudinal\PIV\';
ii   = 4;               % experiment index

DT = 10e-3;             % sec per image pair
DX = 1/17697.69;        % m per pixel  (for depth labels)

% Window lengths to sweep (seconds).
% Cover sub-wave-period → several wave periods so the plateau is visible.
win_lengths_s = [0.05, 0.1, 0.2, 0.5, 1.0, 2.0, 4.0];

% Row indices into the wave-following grid (1 = just below surface).
% Adjust if your grid has a different number of rows.
depth_idx = [3, 8, 15, 25];

% Set true to also export an MP4 (sweeps Delta_t from short to long)
make_video = false;

%% =========================================================================
%% PATHS
%% =========================================================================
DIRS      = dir(LONG);
DIRS      = DIRS(3:end);
exp_name  = DIRS(ii).name;
turb_save = [LONG exp_name '/Chris_recompute/PIVMat_TURB/'];

turb_files = dir([turb_save exp_name '_compVel_*.mat']);
N_frames   = length(turb_files);

if N_frames == 0
    error('No files found in %s\nRun run_decomposition_loop.m first.', turb_save);
end
fprintf('Experiment: %s   Frames in PIVMat_TURB: %d\n', exp_name, N_frames);

%% =========================================================================
%% BUILD u_bar(z, t)  — horizontal mean of intrp_u_res for every frame
%% =========================================================================
% Determine grid depth from frame 1
tmp0 = load(fullfile(turb_save, turb_files(1).name), 'decomposedVel', 'pivRes');
Nz   = size(tmp0.decomposedVel.compVel.intrp_u_res, 1);

u_bar = NaN(Nz, N_frames);

fprintf('Loading %d frames...', N_frames);
for ff = 1:N_frames
    d = load(fullfile(turb_save, turb_files(ff).name), 'decomposedVel');
    u_bar(:, ff) = nanmean(double(d.decomposedVel.compVel.intrp_u_res), 2);
end
fprintf(' done.\n');

t = (0:N_frames-1) * DT;   % seconds

%% =========================================================================
%% SLIDING-WINDOW AVERAGES  (movmean, centred)
%% =========================================================================
depth_idx = depth_idx(depth_idx <= Nz);   % guard against out-of-range
Nd        = numel(depth_idx);
Nw        = numel(win_lengths_s);
win_frames = max(1, round(win_lengths_s / DT));

% u_avg{iw}(id, ff) = sliding average at depth id, frame ff, window iw
u_avg = cell(Nw, 1);
for iw = 1:Nw
    ua = NaN(Nd, N_frames);
    for id = 1:Nd
        ua(id, :) = movmean(u_bar(depth_idx(id), :), win_frames(iw), ...
                            'omitnan', 'Endpoints','shrink');
    end
    u_avg{iw} = ua;
end

%% =========================================================================
%% DEPTH LABELS  (from last loaded pivRes)
%% =========================================================================
pR    = tmp0.pivRes;
z_ax  = (pR.zPIV - pR.GS/2) * DX * 1e3;          % mm, length = Nz (approx)
z_ax  = z_ax(1:min(numel(z_ax), Nz));
dlbl  = arrayfun(@(d) sprintf('z \\approx %.1f mm', z_ax(min(d,numel(z_ax)))), ...
                 depth_idx, 'UniformOutput', false);

%% =========================================================================
%% PRIMARY FIGURE — tiledlayout, one tile per depth
%% =========================================================================
cmap   = parula(Nw);
lw_vec = linspace(0.8, 2.5, Nw);

fig = figure('Color','white','Position',[80 80 1100 min(900, 200*Nd)]);
tlo = tiledlayout(Nd, 1, 'TileSpacing','compact','Padding','compact');

for id = 1:Nd
    ax = nexttile;
    hold on; box on; grid on;

    for iw = 1:Nw
        plot(t, u_avg{iw}(id,:), 'Color', cmap(iw,:), 'LineWidth', lw_vec(iw));
    end

    ylabel('\langle u \rangle (m/s)', 'FontSize', 9)
    title(dlbl{id}, 'FontSize', 9, 'FontWeight', 'normal')
    xlim([t(1), t(end)])

    if id == 1
        leg_str = arrayfun(@(w) sprintf('\\Deltat = %.2f s', w), ...
                           win_lengths_s, 'UniformOutput', false);
        legend(leg_str, 'Location','eastoutside', 'FontSize', 7.5)
    end
    if id < Nd
        set(ax, 'XTickLabel', [])
    end
end

xlabel(tlo, 't (s)', 'FontSize', 10)
drawnow

%% =========================================================================
%% OPTIONAL VIDEO — sweep Delta_t from short to long (single depth panel)
%% =========================================================================
if make_video
    vid_file = fullfile(turb_save, sprintf('check2_not_swaying_%s.mp4', exp_name));
    vw = VideoWriter(vid_file, 'MPEG-4');
    vw.FrameRate = 4;   % 4 fps → viewer has time to read each curve
    open(vw);

    figV = figure('Color','white','Position',[200 200 900 350]);

    for iw = 1:Nw
        clf(figV); hold on; box on; grid on;

        % ghost of all shorter windows (light grey)
        for jw = 1:iw-1
            plot(t, u_avg{jw}(1,:), 'Color', [0.8 0.8 0.8], 'LineWidth', 0.6)
        end
        % current window in colour (use first depth for the video)
        plot(t, u_avg{iw}(1,:), 'Color', cmap(iw,:), 'LineWidth', 2.5)

        xlabel('t (s)'); ylabel('\langle u \rangle (m/s)')
        title(sprintf('\\Deltat = %.2f s  |  %s  |  %s', ...
              win_lengths_s(iw), dlbl{1}, strrep(exp_name,'_','\_')))
        xlim([t(1), t(end)])
        drawnow

        writeVideo(vw, getframe(figV));
    end

    close(vw);
    close(figV);
    fprintf('Video saved to %s\n', vid_file);
end
