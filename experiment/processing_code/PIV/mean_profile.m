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
%% LINE PLOT — sliding-window averages at selected depths
%% =========================================================================
win_lengths_s = DT_pair*[1,2,5,10,20,40];
depth_idx     = [3, 8, 15, 25];

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
