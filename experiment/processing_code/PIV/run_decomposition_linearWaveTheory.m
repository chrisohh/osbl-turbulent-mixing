% run_decomposition_linearWaveTheory.m
% Step-by-step Method A (linear wave theory) wave-turbulence decomposition
% for longitudinal PIV. Plots at every step for sanity checking.
%
% Uses existing Fabrice functions: FindSurfaceCapillary,
% generateTransfo_LC_noLFV_2023, transformVelField_decay_forFab,
% reverseTransformVelField_decay_forFab.

clear; clc;

%% =========================================================================
%% PARAMETERS
%% =========================================================================
LONG = 'D:\DelawareDataBackup\Longitudinal\PIV\';

% Add paths to Fabrice's functions
addpath('D:\Scripps\GC-Wave-Gen\M-Files_FabMarcNovDec2014\');
addpath('D:\Scripps\GC-Wave-Gen\M-Files_FabMarcNovDec2014\FabriceScripts\');
addpath('D:\Scripps\GC-Wave-Gen\M-Files_FabMarcNovDec2014\CrapperOptimizedFindSurface\');

ii = 4;                     % experiment index
image_pair_number = 216;    % sanity-check frame

DX = 1/17697.69;           % m per pixel
DT = 10e-3;                % sec per image pair

% Pixel offset: maps surface-image coords to PIV-image coords
% (experiment-specific; from Main_LC.m)
SU_OFFSET = -1716 + 287;

%% =========================================================================
%% LOAD DATA
%% =========================================================================
DIRS = dir(LONG);
DIRS = DIRS(3:end);
exp_name = DIRS(ii).name;
num_of_digits = 3;
load_path = [LONG exp_name];
pair_str = sprintf(['%0' num2str(num_of_digits) 'd'], image_pair_number);

fprintf('Experiment: %s   Pair: %s\n', exp_name, pair_str);

% --- Raw PIV image (for Step 0/1 plotting) ---
load([load_path '/PIVRaw/PIV/' exp_name '_Piv_' pair_str '_a.mat']);
IM_a = imgPiv;
[h, w] = size(IM_a);

% --- Surface detection ---
imSurfa = FindSurfaceCapillary( ...
    [load_path '/PIVRaw/PIVSURF/' exp_name '_Pivsurf_' pair_str '_a.mat'], ...
    findMask = true);

% --- Load or compute velocity ---
compvel_file = [load_path '/PIVMat_2023/' exp_name '_compVel_' pair_str '.mat'];
if exist(compvel_file, 'file')
    fprintf('Loading compVel from %s\n', compvel_file);
    load(compvel_file);
else
    fprintf('compVel not found — computing from raw images...\n');
    load([load_path '/PIVRaw/PIV/' exp_name '_Piv_' pair_str '_b.mat']);
    IM_b = imgPiv;
    imSurfb = FindSurfaceCapillary( ...
        [load_path '/PIVRaw/PIVSURF/' exp_name '_Pivsurf_' pair_str '_b.mat'], ...
        findMask = true);

    IntrWndw = [256 128 64 24 16 8];
    GrdSpc   = [128  64 32 12  8 4];
    compVel = ComputeVelocities_Quick_NoFilt_Deform_Water( ...
        IM_a, IM_b, imSurfa.mask, imSurfb.mask, IntrWndw, GrdSpc);
    compVel.DX = DX;
    compVel.DT = DT;
end

% --- Extract fields ---
u = compVel.delta_x .* compVel.Mask;   % [Nz_piv x Nx_piv] pixels
w_vel = compVel.delta_z .* compVel.Mask;
[Nz_piv, Nx_piv] = size(u);
Surface_PIV = imSurfa.surfacePIVImg;    % [1 x w] pixel row coords

fprintf('PIV grid: %d x %d   Image: %d x %d\n', Nz_piv, Nx_piv, h, w);

%% =========================================================================
%% STEP 0 — RAW SURFACE IMAGE + DETECTED SURFACE
%% =========================================================================
cl = [1, 0.4, 0.4];

figure(1); clf;
imagesc(imSurfa.ImgScaledToPIVSmallCrop, [0, 300])
colormap bone
hold on
plot(imSurfa.surfaceSurfImgScaled, '-', 'Color', cl, 'LineWidth', 1.5)
plot(imSurfa.surface_raw, '-r', 'LineWidth', 1.5)
daspect([1, 1, 1])
legend('Filtered surface', 'Raw surface', 'Location', 'best')
title(sprintf('Step 0: Surface detection — %s pair %s', exp_name, pair_str), ...
    'Interpreter', 'none')
drawnow

%% =========================================================================
%% STEP 1 — CARTESIAN VELOCITY + SURFACE ON PIV IMAGE
%% =========================================================================
figure(2); clf;

subplot(2,2,1)
imagesc(IM_a, [0, 300])
colormap(gca, gray)
hold on
plot(Surface_PIV, '-', 'Color', cl, 'LineWidth', 1.5)
daspect([1,1,1])
title('PIV image A + surface')

subplot(2,2,2)
imagesc((1:w)*DX*1e3, (1:h)*DX*1e3, IM_a, [0, 300])
colormap(gca, gray)
hold on
plot((1:length(Surface_PIV))*DX*1e3, Surface_PIV*DX*1e3, '-', 'Color', cl, 'LineWidth', 1.5)
daspect([1,1,1])
xlabel('x (mm)'); ylabel('z (mm)')
title('PIV image A (physical coords)')

subplot(2,2,3)
imagesc((1:w)*DX*1e3, (1:h)*DX*1e3, u * DX/DT)
colorbar; colormap(gca, gray)
xlabel('x (mm)'); ylabel('z (mm)')
title('u (m/s)')
daspect([1,1,1])

subplot(2,2,4)
imagesc((1:w)*DX*1e3, (1:h)*DX*1e3, w_vel * DX/DT)
colorbar; colormap(gca, gray)
xlabel('x (mm)'); ylabel('z (mm)')
title('w (m/s)')
daspect([1,1,1])

sgtitle(sprintf('Step 1: Cartesian velocity — %s pair %s', exp_name, pair_str), ...
    'Interpreter', 'none')
drawnow

%% =========================================================================
%% STEP 2 — COORDINATE TRANSFORM + METHOD A ORBITALS
%% =========================================================================
pivRes.zPIV = compVel.zPIV;
pivRes.xPIV = compVel.xPIV;
pivRes.GS   = compVel.GS;
pivRes.mask  = compVel.Mask;

transfo = generateTransfo_LC_noLFV_2023(compVel, Surface_PIV, pivRes);

SU   = transfo.SU(2:end,:);        % remove zeta=0 row (the surface itself)
SU   = SU + SU_OFFSET;             % surface-image → PIV-image pixel coords
ORBX = transfo.ORBX(2:end,:);      % horizontal orbital velocity
ORBZ = transfo.ORBZ(2:end,:);      % vertical orbital velocity

pivRes.pf_surf = SU(1,:);

fprintf('SU size: %d x %d   ORBX size: %d x %d\n', size(SU), size(ORBX));

% --- Plot SU grid lines on Cartesian velocity ---
figure(3); clf;
imagesc((1:w)*DX*1e3, (1:h)*DX*1e3, u * DX/DT)
colorbar; colormap gray; hold on

% Plot every 20th constant-zeta line
n_su_lines = size(SU, 1);
line_skip = max(1, round(n_su_lines / 15));
for iz = 1:line_skip:n_su_lines
    plot((1:size(SU,2)) * compVel.GS * DX * 1e3, SU(iz,:) * DX * 1e3, '-r', 'LineWidth', 0.5)
end
% Surface line
plot((1:length(Surface_PIV))*DX*1e3, Surface_PIV*DX*1e3, '-', 'Color', cl, 'LineWidth', 2)

xlabel('x (mm)'); ylabel('z (mm)')
daspect([1,1,1])
title(sprintf('Step 2: Constant-\\zeta lines on u (m/s) — pair %s', pair_str))
drawnow

%% =========================================================================
%% STEP 3 — TRANSFORM VELOCITY TO WAVE-FOLLOWING
%% =========================================================================
intrp_u = transformVelField_decay_forFab(u, pivRes, SU);
intrp_w = transformVelField_decay_forFab(w_vel, pivRes, SU);

figure(4); clf;
subplot(1,2,1)
imagesc(intrp_u * DX/DT)
colorbar; colormap gray
xlabel('x (PIV columns)'); ylabel('\zeta (levels)')
title('u in wave-following (m/s)')

subplot(1,2,2)
imagesc(intrp_w * DX/DT)
colorbar; colormap gray
xlabel('x (PIV columns)'); ylabel('\zeta (levels)')
title('w in wave-following (m/s)')

sgtitle(sprintf('Step 3: Wave-following velocity — pair %s', pair_str))
drawnow

%% =========================================================================
%% STEP 4 — METHOD A ORBITALS
%% =========================================================================
figure(5); clf;
subplot(1,2,1)
imagesc(ORBX * DX/DT)
colorbar; colormap gray
xlabel('x (PIV columns)'); ylabel('\zeta (levels)')
title('ORBX — u orbital (m/s)')

subplot(1,2,2)
imagesc(ORBZ * DX/DT)
colorbar; colormap gray
xlabel('x (PIV columns)'); ylabel('\zeta (levels)')
title('ORBZ — w orbital (m/s)')

sgtitle(sprintf('Step 4: Method A (linear theory) orbital velocities — pair %s', pair_str))
drawnow

%% =========================================================================
%% STEP 5 — SUBTRACT ORBITALS → RESIDUAL
%% =========================================================================
u_res = intrp_u - ORBX;
w_res = intrp_w - ORBZ;

figure(6); clf;
% Row 1: measured (wave-following)
subplot(3,2,1)
imagesc(intrp_u * DX/DT); colorbar; colormap gray
title('Measured u (wave-following)')
ylabel('Measured')

subplot(3,2,2)
imagesc(intrp_w * DX/DT); colorbar; colormap gray
title('Measured w (wave-following)')

% Row 2: orbitals
subplot(3,2,3)
imagesc(ORBX * DX/DT); colorbar; colormap gray
title('ORBX (Method A)')
ylabel('Orbital')

subplot(3,2,4)
imagesc(ORBZ * DX/DT); colorbar; colormap gray
title('ORBZ (Method A)')

% Row 3: residual
subplot(3,2,5)
imagesc(u_res * DX/DT); colorbar; colormap gray
title('u residual (mean + turb)')
ylabel('Residual')
xlabel('x');

subplot(3,2,6)
imagesc(w_res * DX/DT); colorbar; colormap gray
title('w residual (mean + turb)')
xlabel('x');

sgtitle(sprintf('Step 5: Measured - Orbital = Residual — pair %s', pair_str))
drawnow

%% =========================================================================
%% STEP 6 — TRANSFORM RESIDUAL BACK TO CARTESIAN
%% =========================================================================
u_turb_cart = reverseTransformVelField_decay_forFab(u_res, pivRes, SU);
w_turb_cart = reverseTransformVelField_decay_forFab(w_res, pivRes, SU);

figure(7); clf;
subplot(2,2,1)
imagesc(u * DX/DT); colorbar; colormap gray
title('Original u (Cartesian)')

subplot(2,2,2)
imagesc(w_vel * DX/DT); colorbar; colormap gray
title('Original w (Cartesian)')

subplot(2,2,3)
imagesc(u_turb_cart * DX/DT); colorbar; colormap gray
title('u residual (back to Cartesian)')

subplot(2,2,4)
imagesc(w_turb_cart * DX/DT); colorbar; colormap gray
title('w residual (back to Cartesian)')

sgtitle(sprintf('Step 6: Cartesian comparison — pair %s', pair_str))
drawnow

fprintf('\nDone. Inspect figures 1-7 for sanity checks.\n');
