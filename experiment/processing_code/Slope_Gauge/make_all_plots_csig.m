%% MAKE_ALL_PLOTS_CSIG  Veron-style wave-field plotting suite from CSIG eta(x,y,t).
%
% Loads CISG_slopes_CoreView_22.mat (produced by cisg_main.m), reconstructs
% eta from (Sx, Sy), and renders every CSIG-derivable figure from
% D:\Scripps\OSBL-notes\plotting_plan.tex categories A, B, C, E.3.
%
% Stage windows are deferred — first cut runs on the full record. Hardcode
% stage onset times later and add overlays.

clear; clc;

addpath(fullfile(pwd, 'util'));
addpath(fullfile(pwd, 'plotting'));

mat_file = 'Y:\HLAB_2026\SlopeGauge\Processed\CISG_slopes_CoreView_22.mat';
fprintf('Loading %s ...\n', mat_file);
S = load(mat_file);

Sx    = S.Sx;
Sy    = S.Sy;
time  = S.time;
setup = S.setup;
clear S;
%%
dx = setup.cm_per_pixel_x / 100;     % m
dy = setup.cm_per_pixel_y / 100;     % m

% Zero outside ROI for FFT cleanliness
mask = setup.valid_mask;
Sx(~mask) = 0;
Sy(~mask) = 0;

fprintf('Reconstructing eta from slopes ...\n');
eta = cisg_reconstruct_eta(Sx, Sy, dx, dy);

if ~exist('figures', 'dir'), mkdir figures; end

data = struct('Sx', Sx, 'Sy', Sy, 'eta', eta, 'time', time, ...
              'setup', setup, 'dx', dx, 'dy', dy);

fprintf('\n=== Generating Veron-style CSIG plots ===\n');
%% plot_A1_snapshots_csig(data);

%%

plot_A2_spectrograms_csig(data);
plot_A3_kspectra_csig(data);
plot_A4_slope_variance(data);
plot_A6_directional_B(data);
plot_B1_kx_omega_csig(data);
plot_B3_phase_speed_dev(data);
plot_B4_angular_dist(data);
plot_C1_local_k_map(data);
plot_C4_cross_stream_coh(data);
plot_C5_heterogeneity_idx(data);
plot_C6_local_stokes_map(data);
plot_E3_bandwidth_correction(data);

fprintf('\nAll CSIG plots written to figures/\n');
