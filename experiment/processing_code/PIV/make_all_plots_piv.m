%% MAKE_ALL_PLOTS_PIV  Veron-style wave-field plots from PIV-extracted eta(x,t).
%
% Loads EXTRACTED_SURFACES/A.mat (variable A1) and the matching x-grid file,
% then renders every PIV-eta-derivable figure from
% D:\Scripps\OSBL-notes\plotting_plan.tex (subset of categories A, B, C, E).
%
% Stage windows are deferred — first cut runs on the full record. Hardcode
% stage onset times later and add overlays.

clear; clc;

addpath(fullfile(pwd, 'util'));
addpath(fullfile(pwd, 'plotting'));

root      = get_server_root();
dirin_eta = [root 'Longitudinal\PIV\ExpLCL_2_01\PIVRaw\EXTRACTED_SURFACES\'];

fprintf('Loading eta from %s ...\n', dirin_eta);
load([dirin_eta 'A.mat']);                                % A1 [Nx_eta x Nt_eta]  m
load([dirin_eta 'ExpLCL_2_01_Surface_000.mat']);          % surface_x1

x_eta        = surface_x1(:);
Fs           = 7.2;          % Hz, PIV pair rate
delay        = 73 + 5;       % s  (imaging starts delay s after trigger)
t_wind_onset = 40;           % s  trigger -> wind onset (JFM time origin)
t_eta = delay + (0:size(A1,2)-1)/Fs - t_wind_onset;
eta   = A1;
clear A1 surface_x1;

fprintf('eta: Nx=%d  Nt=%d  t = %.1f .. %.1f s\n', ...
        size(eta,1), size(eta,2), t_eta(1), t_eta(end));

if ~exist('figures', 'dir'), mkdir figures; end

data = struct('eta', eta, 'x_eta', x_eta, 't_eta', t_eta, 'Fs', Fs);

fprintf('\n=== Generating Veron-style PIV eta(x,t) plots ===\n');

plot_A1_eta_snapshots_piv(data);
plot_A2_spectrogram_piv(data);
plot_A3_kspectrum_piv(data);
plot_A5_eps_tilde(data);
plot_B2_kx_omega_piv(data);
plot_C2_local_k_hovmoller(data);
plot_C3_envelope_hovmoller(data);
plot_E1_stokes_profile(data);

fprintf('\nAll PIV plots written to figures/\n');
