clear; clc;

%% --- PATHS ---
addpath('\\Airseaserver41\D\DelawareDataBackup\codes\202204 - IR PIV to compute surface velocity\');
addpath('\\Airseaserver41\D\DelawareDataBackup\codes\202204 - Transverse structures PIV\');

%% --- PARAMETERS ---
params.g     = 9.81;
params.gamma = 7.2e-5;      % m^3/s^2 (surface tension / density)
params.Fs    = 7.2;          % Hz PIV frame rate
params.x_piv = 12.0;        % m fetch of cross-wind PIV plane
t_image_offset = 38;         % s (first frame in JFM coordinates)

%% --- LOAD ETA (Ramp 2, Exp 1) ---
dirin_eta = '\\Airseaserver41\D\DelawareDataBackup\Longitudinal\PIV\ExpLCL_2_01\PIVRaw\EXTRACTED_SURFACES\';
Fs_eta    = 7.2;
delay_eta = 73 + 5;

load([dirin_eta, 'A.mat']);
load([dirin_eta, 'ExpLCL_2_01_Surface_000.mat']);

x_eta = surface_x1(:);                                   % [Nx x 1] m
t_eta = delay_eta + (0:size(A1,2)-1) / Fs_eta - 40;     % JFM time
eta   = A1;                                               % [Nx x Nt]
clear A1 surface_x1;

%% --- LOAD LONGITUDINAL PIV DATA (Ramp 2, Exp 1) ---
% Adjust this path if the longitudinal velocity data lives elsewhere.
dir_piv_long = '\\Airseaserver41\D\DelawareDataBackup\Longitudinal\PIV\ExpLCL_2_01\PIVMat\';

D_tmp = dir([dir_piv_long, '*.mat']);
tmp   = load([dir_piv_long, D_tmp(1).name], 'compVel');
DX    = tmp.compVel.DX;
GS    = tmp.compVel.GS;
clear tmp;

% Load all longitudinal frames.
% If process_one_directory_transverse_PIV works for longitudinal too, use it.
% Otherwise replace with the appropriate loader.
STAT = process_one_directory_transverse_PIV(dir_piv_long);

% In the longitudinal (x-z) PIV plane:
%   compVel.delx → along-wind (x)  → horizontal velocity u
%   compVel.dely → vertical   (z)  → vertical velocity   w
u_long = STAT.u;   % [Nt x Nx x Nz]  along-wind velocity
w_long = STAT.v;   % [Nt x Nx x Nz]  vertical velocity

[Nt, Nx_piv, Nz_piv] = size(w_long);

% --- PIV coordinate vectors ---
t_piv = t_image_offset + (0:Nt-1) / params.Fs;   % JFM time (s)

% z: 0 at surface, negative downward
z_piv = -(0:Nz_piv-1)' * DX;   % [Nz x 1] m

% x: along-wind (relative to image, convert to physical if DX known)
x_piv = (0:Nx_piv-1)' * DX;    % [Nx x 1] m  (relative to left edge of FOV)

fprintf('Longitudinal PIV: Nt=%d, Nx=%d, Nz=%d, DX=%.4f m\n', Nt, Nx_piv, Nz_piv, DX);

%% --- SELECT A SNAPSHOT ---
% Choose a time in the developed LC turbulence phase.
% t ~ 55 s is well into the Langmuir turbulence regime.
t_snap = 55;
[~, it_snap] = min(abs(t_piv - t_snap));
fprintf('Selected snapshot: t = %.2f s (frame %d)\n', t_piv(it_snap), it_snap);

% Extract single-frame Cartesian fields [Nz x Nx] (transpose from [Nx x Nz])
u_cart = squeeze(u_long(it_snap, :, :))';   % [Nz x Nx]
w_cart = squeeze(w_long(it_snap, :, :))';   % [Nz x Nx]

% Surface elevation at this time step (interpolated from eta data)
[~, it_eta] = min(abs(t_eta - t_piv(it_snap)));
eta_snap = interp1(x_eta, eta(:, it_eta), x_piv, 'linear', 'extrap');
eta_snap = eta_snap(:)';  % [1 x Nx]

%% === STEP 1: WAVE-FOLLOWING COORDINATE TRANSFORM (WATER SIDE) ===
% zeta: depth below local surface, positive downward.
% Grid spacing matches PIV resolution, starting at surface.
dz = abs(z_piv(2) - z_piv(1));
zeta_vec = (0:dz:abs(z_piv(end)))';   % [Nz x 1] m

transfo = generate_transfo_water(eta_snap, x_piv, zeta_vec);

% Transform velocity fields to wave-following coordinates
u_wf = transform_field_to_wavefollowing(u_cart, z_piv, transfo);
w_wf = transform_field_to_wavefollowing(w_cart, z_piv, transfo);

fprintf('Wave-following transform done: %d x %d grid\n', size(u_wf));

%% === STEP 2: PHASE & ENSEMBLE AVERAGING (all snapshots) ===
% Bin each x-column of each frame by its wave phase, then average.
% This is the Fabio approach (phiAvNensAv_fab) adapted for our data format.

Nbins = 180;  % phase bins (same as Fabio default)
phase_edges = linspace(-pi, pi, Nbins + 1);

% Accumulators for phase average and ensemble average
u_phaseSum = zeros(length(zeta_vec), Nbins);
w_phaseSum = zeros(length(zeta_vec), Nbins);
phaseCounts = zeros(length(zeta_vec), Nbins);

u_ensembleSum = zeros(length(zeta_vec), 1);
w_ensembleSum = zeros(length(zeta_vec), 1);
ensembleCounts = zeros(length(zeta_vec), 1);

fprintf('Computing phase & ensemble averages over %d frames...\n', Nt);
for it = 1:Nt
    % Extract Cartesian fields for this frame
    u_t = squeeze(u_long(it, :, :))';   % [Nz x Nx]
    w_t = squeeze(w_long(it, :, :))';

    % Surface elevation for this frame
    [~, it_e] = min(abs(t_eta - t_piv(it)));
    eta_t = interp1(x_eta, eta(:, it_e), x_piv, 'linear', 'extrap');
    eta_t = eta_t(:)';

    % Wave-following transform
    transfo_t = generate_transfo_water(eta_t, x_piv, zeta_vec);
    u_wf_t = transform_field_to_wavefollowing(u_t, z_piv, transfo_t);
    w_wf_t = transform_field_to_wavefollowing(w_t, z_piv, transfo_t);

    % Phase at each x from Hilbert of this snapshot's eta
    phase_t = transfo_t.phase;   % [1 x Nx]

    % Bin by phase
    bin_idx = discretize(phase_t, phase_edges);

    for ib = 1:Nbins
        cols = find(bin_idx == ib);
        if isempty(cols), continue; end

        u_chunk = u_wf_t(:, cols);
        w_chunk = w_wf_t(:, cols);

        u_phaseSum(:, ib) = nansum([u_phaseSum(:, ib), nansum(u_chunk, 2)], 2);
        w_phaseSum(:, ib) = nansum([w_phaseSum(:, ib), nansum(w_chunk, 2)], 2);
        phaseCounts(:, ib) = phaseCounts(:, ib) + sum(isfinite(u_chunk), 2);
    end

    % Ensemble sum (all columns regardless of phase)
    u_ensembleSum = u_ensembleSum + nansum(u_wf_t, 2);
    w_ensembleSum = w_ensembleSum + nansum(w_wf_t, 2);
    ensembleCounts = ensembleCounts + sum(isfinite(u_wf_t), 2);

    if mod(it, round(Nt/10)) == 0
        fprintf('  %d%% done\n', round(100*it/Nt));
    end
end

% --- Averages ---
u_phaseAvg    = u_phaseSum ./ phaseCounts;         % [Nz x Nbins]
w_phaseAvg    = w_phaseSum ./ phaseCounts;
u_ensembleAvg = u_ensembleSum ./ ensembleCounts;   % [Nz x 1]
w_ensembleAvg = w_ensembleSum ./ ensembleCounts;

fprintf('Phase average: %d bins, min %d max %d counts per bin\n', ...
    Nbins, min(phaseCounts(:)), max(phaseCounts(:)));

%% === STEP 3: DECOMPOSE SNAPSHOT INTO TURBULENT & WAVE-COHERENT ===
% For the selected snapshot:
%   turbulent:      u' = u - <u>(phase)      (subtract phase average)
%   wave-coherent:  u~ = <u>(phase) - <u>    (phase average minus ensemble)

phase_snap = transfo.phase;
bin_snap = discretize(phase_snap, phase_edges);

u_turb = NaN(size(u_wf));
w_turb = NaN(size(w_wf));
u_wave = NaN(size(u_wf));
w_wave = NaN(size(w_wf));

for ix = 1:Nx_piv
    ib = bin_snap(ix);
    if isnan(ib), continue; end

    u_turb(:, ix) = u_wf(:, ix) - u_phaseAvg(:, ib);
    w_turb(:, ix) = w_wf(:, ix) - w_phaseAvg(:, ib);

    u_wave(:, ix) = u_phaseAvg(:, ib) - u_ensembleAvg;
    w_wave(:, ix) = w_phaseAvg(:, ib) - w_ensembleAvg;
end

%% === STEP 4: INVERSE TRANSFORM TURBULENT FIELD BACK TO CARTESIAN ===
u_turb_cart = inverse_transform_to_cartesian(u_turb, z_piv, transfo);
w_turb_cart = inverse_transform_to_cartesian(w_turb, z_piv, transfo);
u_wave_cart = inverse_transform_to_cartesian(u_wave, z_piv, transfo);
w_wave_cart = inverse_transform_to_cartesian(w_wave, z_piv, transfo);

%% --- SAVE ---
save('longitudinal_wave_decomposition.mat', ...
    'u_cart', 'w_cart', 'u_wf', 'w_wf', ...
    'u_phaseAvg', 'w_phaseAvg', 'u_ensembleAvg', 'w_ensembleAvg', ...
    'u_turb', 'w_turb', 'u_wave', 'w_wave', ...
    'u_turb_cart', 'w_turb_cart', 'u_wave_cart', 'w_wave_cart', ...
    'transfo', 'zeta_vec', 'z_piv', 'x_piv', 't_piv', 'params', ...
    'it_snap', 'Nbins', 'phase_edges', 'phaseCounts');
disp('Saved longitudinal_wave_decomposition.mat');

%% === PLOTS ===
plot_wave_decomposition(u_cart, w_cart, u_wf, w_wf, ...
    u_phaseAvg, w_phaseAvg, u_ensembleAvg, w_ensembleAvg, ...
    u_turb, w_turb, u_wave, w_wave, ...
    transfo, zeta_vec, z_piv, x_piv, t_piv(it_snap), ...
    Nbins, phase_edges);
