% Help you select mean profile time window
% Loads intrp_u_res from PIVMat_TURB (built by run_decomposition_loop.m),
% computes the horizontal-mean <u>(z, t) and plots as a Hovmoller diagram.

clear; clc;

%% =========================================================================
%% PARAMETERS
%% =========================================================================
p = get_piv_params('ExpLCL_2_01');

exp_name  = p.exp_name;
turb_save = strcat([p.turb_save, '\']);
DX        = p.DX;           % 1/17697.69 m/pixel
DT        = p.DT;           % 10e-3 s, inter-frame (A→B), for velocity scaling
Fs_cam    = p.Fs_cam;       % 14.4 Hz — camera acquisition rate (images)
Fs_PIV    = p.Fs_PIV;       % 7.2 Hz  — velocity map rate (one pair = two images)
DT_pair   = p.DT_pair;      % 1/7.2 s — time between consecutive pairs (for time axis)

turb_files = dir([turb_save exp_name '_compVel_*.mat']);
N_frames   = length(turb_files);

if N_frames == 0
    error('No files found in %s\nRun run_decomposition_loop.m first.', turb_save);
end
fprintf('Experiment: %s   Frames: %d\n', exp_name, N_frames);

%% =========================================================================
%% BUILD u_bar(z, t)
%% =========================================================================
tmp0 = load(fullfile(turb_save, turb_files(1).name), 'decomposedVel', 'pivRes');
Nz   = size(tmp0.decomposedVel.compVel.intrp_u_res, 1);

u_bar = NaN(Nz, N_frames);

fprintf('Loading %d frames...', N_frames);
for ff = 1:N_frames
    d = load(fullfile(turb_save, turb_files(ff).name), 'decomposedVel');
    u_bar(:, ff) = nanmean(double(d.decomposedVel.compVel.intrp_u_res), 2);
end
fprintf(' done.\n');

t  = p.t_imaging_wrt_ramp + (0:N_frames-1) * DT_pair;   % s since wind ramp-up (t=0 = wind on)

% z_ax in curvilinear (wave-following) coordinates: ζ = 0 is the surface,
% increasing downward. Row i of intrp_u_res is i*GS pixels below the surface.
GS   = tmp0.pivRes.GS;
z_ax = (1:Nz) * GS * DX;   % m below surface (ζ)

%% =========================================================================
%% HOVMOLLER — raw u_bar
%% =========================================================================
fig = figure('Color','white','Position',[80 80 1100 420]);
set(fig, 'DefaultAxesColor', [0.85 0.85 0.85]);   % gray for NaN/masked
h=imagesc(t, z_ax, u_bar);
set(h, 'AlphaData', ~isnan(u_bar));
axis xy;
colormap(brewermap([],'Spectral'));
cb = colorbar;
cb.Label.String = '\langle u_{raw}-u_{orb} \rangle (m/s)';
clim([-0.01, 0.12]);
xlabel('t (s)');
ylabel('\zeta (m)');
title(sprintf('Horizontal-mean u_{raw}-u_{orb}  |  %s', strrep(exp_name,'_','\_')));
ylim([0,0.1])
set(gca, 'YDir', 'reverse');

%% =========================================================================
%% ENSEMBLE AVERAGE across runs
%% =========================================================================
ens_runs = {'ExpLCL_2_01', 'ExpLCL_2_02', 'ExpLCL_2_03'};
N_runs   = numel(ens_runs);

ub = cell(N_runs, 1);   % u_bar(z, t)    per run (horizontal mean)
uc = cell(N_runs, 1);   % u_center(z, t) per run (center-x column)
tt = cell(N_runs, 1);   % time axis (s since wind ramp) per run
for ir = 1:N_runs
    [ub{ir}, tt{ir}, ~, uc{ir}] = load_u_bar(ens_runs{ir});
    fprintf('  %s: %d frames\n', ens_runs{ir}, size(ub{ir}, 2));
end

%%
fig = figure('Color','white','Position',[80 80 1100 1100]);
set(fig, 'DefaultAxesColor', [0.85 0.85 0.85]);   % gray for NaN/masked
for ir=1:N_runs
subplot(3,1,ir);imagesc(ub{ir,1})
h=imagesc(tt{ir,1}, z_ax*1e3, ub{ir,1});
set(h, 'AlphaData', ~isnan(u_bar));
axis xy;
colormap(brewermap([],'Spectral'));
cb = colorbar;
cb.Label.String = '\langle u_{raw}-u_{orb} \rangle (m/s)';
clim([-0.01, 0.12]);

ylabel('\zeta (mm)');
ylim([0,0.1*1e3])
set(gca, 'YDir', 'reverse');
title(sprintf('Run %d', ir));
end
xlabel('t (s)');

%%
% Common time grid = overlap window of all runs (each run has its own
% t_imaging_wrt_ramp offset, so align on physical time, not frame index).
t0    = max(cellfun(@(x) x(1),   tt));
t1    = min(cellfun(@(x) x(end), tt));
t_ens = t0 : DT_pair : t1;

Nz_common   = min(cellfun(@(x) size(x, 1), ub));
u_stack     = NaN(Nz_common, numel(t_ens), N_runs);
u_cen_stack = NaN(Nz_common, numel(t_ens), N_runs);
for ir = 1:N_runs
    % interp each depth onto the common time grid
    u_stack(:, :, ir)     = interp1(tt{ir}, ub{ir}(1:Nz_common, :).', t_ens).';
    u_cen_stack(:, :, ir) = interp1(tt{ir}, uc{ir}(1:Nz_common, :).', t_ens).';
end
u_bar_ens = mean(u_stack,     3, 'omitnan');   % ensemble mean <u>(z, t)  (horizontal mean)
u_cen_ens = mean(u_cen_stack, 3, 'omitnan');   % ensemble mean u(z, t)    (center-x column)

figure('Color','white','Position',[80 80 900 min(900, 200*Nd)]);
tlo = tiledlayout(Nd, 1, 'TileSpacing','compact','Padding','compact');
for id = 1:Nd
    di = depth_idx(id);
    if di > Nz_common, continue; end
    ax = nexttile; hold on; box on; grid on

    % individual runs (thin, faded)
    hr = gobjects(N_runs, 1);
    for ir = 1:N_runs
        hr(ir) = plot(t_ens, squeeze(u_stack(di, :, ir)), ...
                      'LineWidth', 0.8, 'Color', [0.6 0.6 0.6]);
    end
    % ensemble mean (bold)
    he = plot(t_ens, u_bar_ens(di, :), 'k', 'LineWidth', 2);

    ylabel('\langle u \rangle (m/s)', 'FontSize', 9)
    title(dlbl{id}, 'FontSize', 9, 'FontWeight', 'normal')
    xlim([t_ens(1), t_ens(end)])
    if id == 1
        legend([hr(1), he], {'individual runs', 'ensemble mean'}, ...
               'Location','eastoutside', 'FontSize', 8)
    end
    if id < Nd
        set(ax, 'XTickLabel', [])
    end
end
xlabel(tlo, 't (s)', 'FontSize', 10)
title(tlo, sprintf('Ensemble-mean u over %d runs (%s)', N_runs, ...
      strrep(strjoin(ens_runs, ', '), '_', '\_')), 'FontSize', 10)

%% =========================================================================
%% HOVMOLLER — ensemble mean
%% =========================================================================
z_ens = z_ax(1:Nz_common);
fig = figure('Color','white','Position',[80 80 1100 420]);
set(fig, 'DefaultAxesColor', [0.85 0.85 0.85]);   % gray for NaN/masked
h = imagesc(t_ens, z_ens*1e3, u_bar_ens);
set(h, 'AlphaData', ~isnan(u_bar_ens));
axis xy;
colormap(brewermap([],'Spectral'));
cb = colorbar;
cb.Label.String = '\langle u_{raw}-u_{orb} \rangle (m/s)';
clim([-0.01, 0.12]);
xlabel('t (s)');
ylabel('\zeta (mm)');
title(sprintf('Ensemble-average horizontal-average u  |  %d runs', N_runs));
ylim([0, 0.1*1e3])
set(gca, 'YDir', 'reverse');

%% =========================================================================
%% TIME SERIES — ensemble mean at selected depths (overlaid)
%% =========================================================================
figure('Color','white'); hold on;
for id = 1:Nd
    di = depth_idx(id);
    if di > Nz_common, continue; end
    plot(t_ens, u_bar_ens(di, :), 'LineWidth', 1.4, 'DisplayName', dlbl{id})
end
xlabel('t (s)'); ylabel('\langle u \rangle (m/s)')
xlim([t_ens(1), t_ens(end)])
title(sprintf('Ensemble-mean u at selected depths  |  %d runs', N_runs))
legend('Location','best')


%% =========================================================================
%% LINE PLOT — sliding-window averages at selected depths
%% =========================================================================
win_lengths_s = DT_pair*[1,2,5,10,20,40];
depth_mm = [2, 5, 10, 20,50];
depth_idx = arrayfun(@(d) find(min(abs(z_ax*1e3 - d)) == abs(z_ax*1e3 - d), 1), depth_mm);
% depth_idx     = [2, 5, 10, 20,50];

depth_idx  = depth_idx(depth_idx <= Nz);
Nd         = numel(depth_idx);
Nw         = numel(win_lengths_s);
win_frames = max(1, round(win_lengths_s / DT_pair));

u_avg = cell(Nw, 1);
for iw = 1:Nw
    ua = NaN(Nd, N_frames);
    for id = 1:Nd
        ua(id, :) = movmean(u_bar(depth_idx(id), :), win_frames(iw), ...
                            'omitnan', 'Endpoints', 'shrink');
    end
    u_avg{iw} = ua;
end

dlbl = arrayfun(@(d) sprintf('z \\approx %.1f mm', 1e3*z_ax(min(d, numel(z_ax*1e3)))), ...
                depth_idx, 'UniformOutput', false);

cmap   = parula(Nw);
lw_vec = linspace(0.8, 2.5, Nw);

figure('Color','white','Position',[80 560 1100 min(900, 200*Nd)]);
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
        leg_str = arrayfun(@(w) sprintf('t_{window} = %.2f s', w), ...
                           win_lengths_s, 'UniformOutput', false);
        legend(leg_str, 'Location','eastoutside', 'FontSize', 7.5)
    end
    if id < Nd
        set(ax, 'XTickLabel', [])
    end
end
xlabel(tlo, 't (s)', 'FontSize', 10)
title(tlo, sprintf('Sliding-window mean u  |  %s', strrep(exp_name,'_','\_')), 'FontSize', 10)
drawnow

%% =========================================================================
%% SLIDING-WINDOW HOVMOLLERS (ensemble) — continuous z
%% =========================================================================
% Smooth the ensemble+horizontal mean u_bar_ens with each candidate time
% window and show the full z-t field, so the temporal wobble is visible at
% every depth (not just the few selected ones).
win_frames_ens = max(1, round(win_lengths_s / DT_pair));

figure('Color','white','Position',[60 60 1200 700]);
tloH = tiledlayout(2, ceil(Nw/2), 'TileSpacing','compact','Padding','compact');
for iw = 1:Nw
    uw = movmean(u_bar_ens, win_frames_ens(iw), 2, 'omitnan', 'Endpoints','shrink');
    ax = nexttile;
    set(ax, 'Color', [0.85 0.85 0.85]);   % gray for NaN/masked
    hh = imagesc(t_ens, z_ens*1e3, uw);
    set(hh, 'AlphaData', ~isnan(uw));
    axis xy; set(ax, 'YDir', 'reverse'); ylim([0, 0.1*1e3]);
    colormap(brewermap([],'Spectral')); clim([-0.01, 0.12]);
    title(sprintf('t_{win} = %.2f s', win_lengths_s(iw)), 'FontWeight','normal');
end
cb = colorbar; cb.Layout.Tile = 'east';
cb.Label.String = '\langle u \rangle (m/s)';
xlabel(tloH, 't (s)'); ylabel(tloH, '\zeta (mm)');
title(tloH, sprintf('Sliding-window smoothing of ensemble mean  |  %d runs', N_runs));

%% =========================================================================
%% WOBBLE vs WINDOW — pick a window from the knee
%% =========================================================================
% removed(z, W) = RMS over t of (u - smoothed). Rises steeply while the window
% is eating turbulent/wave wobble, then flattens once the wobble is gone.
% The knee (per depth) is a good window length.
removed = NaN(Nz_common, Nw);
for iw = 1:Nw
    uw = movmean(u_bar_ens, win_frames_ens(iw), 2, 'omitnan', 'Endpoints','shrink');
    removed(:, iw) = sqrt(mean((u_bar_ens - uw).^2, 2, 'omitnan'));
end

figure('Color','white','Position',[80 80 1050 430]);
tloW = tiledlayout(1, 2, 'TileSpacing','compact','Padding','compact');

% (a) z-W heatmap of removed fluctuation
nexttile;
imagesc(win_lengths_s, z_ens*1e3, removed);
axis xy; set(gca, 'YDir', 'reverse'); ylim([0, 0.1*1e3]);
colormap(gca, parula); cbw = colorbar; cbw.Label.String = 'removed RMS (m/s)';
xlabel('t_{window} (s)'); ylabel('\zeta (mm)');
title('Fluctuation removed vs window & depth');

% (b) knee curves at the selected depths
nexttile; hold on; box on; grid on
for id = 1:Nd
    di = depth_idx(id);
    if di > Nz_common, continue; end
    plot(win_lengths_s, removed(di, :), '-o', 'LineWidth', 1.2, 'DisplayName', dlbl{id});
end
xlabel('t_{window} (s)'); ylabel('removed RMS (m/s)');
title('Knee \approx good window'); legend('Location','best', 'FontSize', 8);
drawnow

%% =========================================================================
%% PROFILE EVOLUTION — u(z) profiles colored by time
%% =========================================================================
% Each line is the ensemble-mean profile u(z) at one instant; color = time.
% Lets you see the profile shape and its temporal wobble directly.
t_prof_win = [t_ens(1), t_ens(end)];   % restrict here to zoom a window of interest
win_frames_prof = 1;                   % set >1 to lightly time-smooth first
nprof      = 40;                        % number of profiles to draw

u_prof = movmean(u_bar_ens, win_frames_prof, 2, 'omitnan', 'Endpoints','shrink');
tsel   = linspace(t_prof_win(1), t_prof_win(2), nprof);
[~, ksel] = arrayfun(@(ts) min(abs(t_ens - ts)), tsel);

cmap = parula(nprof);
figure('Color','white','Position',[80 80 560 640]); hold on; box on
for ii = 1:nprof
    plot(u_prof(:, ksel(ii)), z_ens*1e3, 'Color', cmap(ii,:), 'LineWidth', 1.0);
end
set(gca, 'YDir', 'reverse'); ylim([0, 0.1*1e3]);
xlabel('\langle u \rangle (m/s)'); ylabel('\zeta (mm)');
title(sprintf('Ensemble-mean u(z) profiles over time  |  %d runs', N_runs));
colormap(gca, cmap); clim(t_prof_win);
cb = colorbar; cb.Label.String = 't (s)';
drawnow

%% =========================================================================
%% MEAN PROFILE — u(z), time-averaged over a window
%% =========================================================================
t_avg_win = [t_ens(1), t_ens(end)];   % <-- set the averaging window, e.g. [50 65]
kwin  = t_ens >= t_avg_win(1) & t_ens <= t_avg_win(2);

u_prof_mean = mean(u_bar_ens(:, kwin), 2, 'omitnan');   % <u>(z) over the window
u_prof_std  = std(u_bar_ens(:, kwin), 0, 2, 'omitnan'); % temporal wobble at each z

figure('Color','white','Position',[80 80 480 640]); hold on; box on; grid on
% ± temporal std band
zmm = z_ens(:)*1e3;
fill([u_prof_mean-u_prof_std; flipud(u_prof_mean+u_prof_std)], ...
     [zmm; flipud(zmm)], [0.85 0.88 0.95], 'EdgeColor','none', ...
     'DisplayName','\pm std over time');
plot(u_prof_mean, zmm, 'k', 'LineWidth', 2, 'DisplayName','mean');
set(gca, 'YDir', 'reverse'); ylim([0, 0.1*1e3]);
xlabel('\langle u \rangle (m/s)'); ylabel('\zeta (mm)');
title(sprintf('Ensemble-mean u(z),  t = %.1f-%.1f s  |  %d runs', ...
      t_avg_win(1), t_avg_win(2), N_runs));
legend('Location','southeast', 'FontSize', 8);
drawnow

%% =========================================================================
%% PROFILES AT SELECTED TIMES — u(z) at chosen instants, overlaid
%% =========================================================================
inst_times_s = [38 50 55 57 59 65];     % s since wind ramp
[~, inst_k]  = arrayfun(@(ts) min(abs(t_ens - ts)), inst_times_s);

cmap_i = parula(numel(inst_times_s));
figure('Color','white','Position',[80 80 500 660]); hold on; box on; grid on
for ii = 1:numel(inst_k)
    plot(u_bar_ens(:, inst_k(ii)), z_ens*1e3, 'Color', cmap_i(ii,:), ...
         'LineWidth', 1.6, 'DisplayName', sprintf('t \\approx %.1f s', t_ens(inst_k(ii))));
end
plot(0.13, 0, 'rp', 'MarkerSize', 14, 'MarkerFaceColor', 'r', ...
     'DisplayName', '13 cm/s @ surface');
set(gca, 'YDir', 'reverse'); ylim([0, 0.1*1e3]);
xlabel('\langle u \rangle (m/s)'); ylabel('\zeta (mm)');
title(sprintf('Ensemble-mean u(z) at selected times  |  %d runs', N_runs));
legend('Location','southeast', 'FontSize', 8);
drawnow

%% =========================================================================
%% LOCAL FUNCTIONS
%% =========================================================================
function [u_bar, t, z_ax, u_center] = load_u_bar(exp_name)
% Build horizontal-mean <u>(z, t) and center-column u(z, t) for one run.
%   u_bar    : horizontal (x) mean at each depth/frame
%   u_center : single center-x column at each depth/frame (instantaneous)
    p = get_piv_params(exp_name);
    turb_save  = strcat([p.turb_save, '\']);
    turb_files = dir([turb_save exp_name '_compVel_*.mat']);
    N_frames   = length(turb_files);
    if N_frames == 0
        error('No files found in %s\nRun run_decomposition_loop.m first.', turb_save);
    end

    tmp0     = load(fullfile(turb_save, turb_files(1).name), 'decomposedVel', 'pivRes');
    Nz       = size(tmp0.decomposedVel.compVel.intrp_u_res, 1);
    u_bar    = NaN(Nz, N_frames);
    u_center = NaN(Nz, N_frames);
    for ff = 1:N_frames
        d   = load(fullfile(turb_save, turb_files(ff).name), 'decomposedVel');
        uxz = double(d.decomposedVel.compVel.intrp_u_res);
        u_bar(:, ff)    = nanmean(uxz, 2);
        cen             = round(size(uxz, 2) / 2);   % center column in x
        u_center(:, ff) = uxz(:, cen);
    end

    t    = p.t_imaging_wrt_ramp + (0:N_frames-1) * p.DT_pair;
    z_ax = (1:Nz) * tmp0.pivRes.GS * p.DX;
end
