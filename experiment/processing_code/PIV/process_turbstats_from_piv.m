clear; clc;

%% --- PATHS ---
% Smoothn utility (N-D robust spline smoothing)
addpath('\\Airseaserver41\D\DelawareDataBackup\codes\202204 - IR PIV to compute surface velocity\');

% process_one_directory_transverse_PIV
addpath('\\Airseaserver41\D\DelawareDataBackup\codes\202204 - Transverse structures PIV\');

%% --- PARAMETERS ---
params.g     = 9.81;        % m/s^2
params.gamma = 7.2e-5;      % m^3/s^2  (surface tension / density, JFM Eq 3.6)
params.H     = 0.71;        % m  water depth (JFM Section 2)
params.Fs    = 7.2;         % Hz PIV frame rate
params.x_piv = 12.0;        % m  along-wind fetch of cross-wind PIV plane

% --- Ramp 2 / JFM 2023 timing ---
% Trigger at t=0. Wind starts ramping at t=40 s (trigger delay).
% Imaging begins at t = 5 + IID = 5 + 73 = 78 s after trigger.
% JFM paper t=0 is defined as wind onset (t=40 s absolute).
% → first PIV frame is at t_JFM = 78 - 40 = 38 s
t_image_offset = 38;  % s  (time of first PIV frame in JFM coordinates)

%% --- LOAD ETA (Ramp 2, Exp 1) ---
% run_all_R2.m uses hardcoded E:\ paths, so inline the loading with server paths.
% Source: \\Airseaserver41\D\DelawareDataBackup\Longitudinal\PIV\ExpLCL_2_01\

dirin_eta = '\\Airseaserver41\D\DelawareDataBackup\Longitudinal\PIV\ExpLCL_2_01\PIVRaw\EXTRACTED_SURFACES\';
Fs_eta    = 7.2;   % Hz
delay_eta = 73 + 5;  % s  (= 78 s after trigger = 38 s in JFM time)

load([dirin_eta, 'A.mat']);                                    % loads A1 [Nx × Nt_eta]
load([dirin_eta, 'ExpLCL_2_01_Surface_000.mat']);              % loads surface_x1 [Nx × 1]

x_eta = surface_x1(:);                                        % [Nx × 1]  m, along-wind
t_eta = delay_eta + (0:size(A1,2)-1) / Fs_eta - 40;          % JFM time (s since wind onset)
eta   = A1;                                                    % [Nx × Nt_eta]

clear A1 surface_x1;

%% --- LOAD PIV DATA (Ramp 2, Repeat 01) ---
dir_piv = '\\Airseaserver41\D\DelawareDataBackup\Transverse\PIV\ExpLCTB_2_01\PIVMat2\';

% Get spatial scale (DX, GS) from the first compVel file without loading all frames.
D_tmp = dir([dir_piv, '*.mat']);
tmp   = load([dir_piv, D_tmp(1).name], 'compVel');
DX    = tmp.compVel.DX;   % m/pixel
GS    = tmp.compVel.GS;   % grid spacing (pixels); grid cell size = DX*GS
clear tmp D_tmp;

% Load all frames via existing function
STAT = process_one_directory_transverse_PIV(dir_piv);

% --- Velocity orientation in the cross-wind (y–z) PIV plane ---
%   compVel.delx → image horizontal = cross-wind direction (y)  → STAT.u
%   compVel.dely → image vertical   = depth direction    (z)    → STAT.v
v_cross = STAT.u;   % cross-wind velocity  [Nt × Ny × Nz]
w_vert  = STAT.v;   % vertical velocity    [Nt × Ny × Nz]

[Nt, Ny, Nz] = size(w_vert);

% --- PIV time vector (JFM: seconds since wind onset) ---
t_piv = t_image_offset + (0:Nt-1) / params.Fs;

% --- PIV z-grid: 0 at water surface, negative downward ---
% The PIV image top row is at the surface; each pixel row is DX metres deeper.
% Verify: Nz * DX should equal the FOV height (~0.139 m per JFM Section 2).
z_piv = -(0:Nz-1)' * DX;   % [Nz × 1]  m

fprintf('PIV: Nt=%d frames, Ny=%d, Nz=%d pixels, DX=%.4f m\n', Nt, Ny, Nz, DX);
fprintf('     z spans %.1f mm to %.1f mm (expect ~139 mm)\n', ...
        1e3*z_piv(1), 1e3*z_piv(end));

%% --- COMPUTE WAVE ORBITAL VELOCITIES ---
[w_wave, ~, wave_info] = compute_wave_orbital_velocities(eta, x_eta, t_eta, z_piv, params);
% w_wave: [Nx_eta × Nz × Nt_eta]

% Extract at PIV along-wind location (x = 12 m fetch)
[~, ix_piv] = min(abs(x_eta - params.x_piv));
w_wave_at_piv = squeeze(w_wave(ix_piv, :, :));   % [Nz × Nt_eta]

% Interpolate onto PIV time vector.
% Both signals are at 7.2 Hz but may start at slightly different absolute times.
% interp1 extrapolates to 0 outside the eta time range.
w_wave_profile = interp1(t_eta, w_wave_at_piv', t_piv, 'linear', 0)';
% [Nz × Nt_piv]

%% --- METHOD A: HILBERT WAVE SUBTRACTION ---
[v_turb_A, w_turb_A] = remove_wave_from_piv(v_cross, w_vert, 'hilbert', w_wave_profile);
stats_A = compute_horizontal_avg_turbstats(v_turb_A, w_turb_A, z_piv, t_piv);

%% --- METHOD B: Y-MEAN SUBTRACTION ---
[v_turb_B, w_turb_B] = remove_wave_from_piv(v_cross, w_vert, 'hmean');
stats_B = compute_horizontal_avg_turbstats(v_turb_B, w_turb_B, z_piv, t_piv);

%% --- SAVE ---
save('turbstats_JFM2023.mat', 'stats_A', 'stats_B', 'wave_info', 'params');
disp('Done. Saved turbstats_JFM2023.mat');

%% --- PLOT ---
plot_turbstats_profiles(stats_A);
