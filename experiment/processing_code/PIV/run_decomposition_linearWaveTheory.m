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

% rootpath = 'C:\Users\airsealab\Documents\GitHub';
rootpath = 'D:\Scripps';

% Add paths to Fabrice's functions
addpath(strcat(rootpath,'\GC-Wave-Gen\M-Files_FabMarcNovDec2014\'));
addpath(strcat(rootpath,'\GC-Wave-Gen\M-Files_FabMarcNovDec2014\FabriceScripts\'));
addpath(strcat(rootpath,'\GC-Wave-Gen\M-Files_FabMarcNovDec2014\CrapperOptimizedFindSurface\'));

ii = 4;                     % experiment index
image_pair_number = 203;    % sanity-check frame

DX = 1/17697.69;           % m per pixel
DT = 10e-3;                % sec per image pair

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

chris_path = [load_path '/Chris_recompute/'];
piv_save   = [chris_path 'PIVMat/'];
turb_save  = [chris_path 'PIVMat_TURB/'];

compvel_file = [piv_save exp_name '_compVel_' pair_str '.mat'];
if exist(compvel_file, 'file')
    fprintf('Loading compVel from %s\n', compvel_file);
    load(compvel_file,'imSurfa','imSurfb');
else
    % --- Surface detection ---
imSurfa = FindSurfaceCapillary( ...
    [load_path '/PIVRaw/PIVSURF/' exp_name '_Pivsurf_' pair_str '_a.mat'], ...
    findMask = true);
imSurfb = FindSurfaceCapillary( ...
    [load_path '/PIVRaw/PIVSURF/' exp_name '_Pivsurf_' pair_str '_b.mat'], ...
    findMask = true);
end

%% --- Raw PIV image ---
load([load_path '/PIVRaw/PIV/' exp_name '_Piv_' pair_str '_a.mat']);
IM_a = imgPiv;
load([load_path '/PIVRaw/PIV/' exp_name '_Piv_' pair_str '_b.mat']);
IM_b = imgPiv;
[height, width] = size(IM_a);

%% =========================================================================
%% STEP 0 — RAW SURFACE IMAGE + DETECTED SURFACE
%% =========================================================================
%Raw surface camera
figure;
imagesc(imSurfa.ImgScaledToPIVSmallCrop,[0,300]);
colormap bone
axis tight;axis equal
title(sprintf('Surface Camera: %s frame %s', exp_name, pair_str), ...
    'Interpreter', 'none')

%Raw PIV camera
figure;
subplot(1,2,1)
imagesc(IM_a,[0,300]);
hold on;plot(imSurfa.surfacePIVImg, '-r', 'LineWidth', 1)
colormap gray
axis tight;axis equal
title(sprintf('PIV %s frame %s - A', exp_name, pair_str), ...
    'Interpreter', 'none')

subplot(1,2,2)
imagesc(IM_b,[0,300]);
hold on;plot(imSurfb.surfacePIVImg, '-r', 'LineWidth', 1)
colormap gray
axis tight;axis equal
title(sprintf('PIV %s frame %s - B', exp_name, pair_str), ...
    'Interpreter', 'none')

set(gcf,'Color','white')
drawnow

%% =========================================================================
%% STEP 1 — CARTESIAN VELOCITY + SURFACE ON PIV IMAGE
%% =========================================================================
%% --- Load or compute velocity ---
rerun = 1;

compvel_file = [piv_save exp_name '_compVel_' pair_str '.mat'];
if exist(compvel_file, 'file') && rerun==0
    fprintf('Loading compVel from %s\n', compvel_file);
    load(compvel_file,'compVel');
else
    fprintf('compVel not found — computing from raw images...\n');
%     IntrWndw = [256 128 64 32]; %[256 128 64 32 16 8];
%     GrdSpc   = IntrWndw/2;
%     compVel = ComputeVelocities_Quick_NoFilt_Deform_Water( ...
%         IM_a, IM_b, imSurfa.mask, imSurfb.mask, IntrWndw, GrdSpc);
    %%
% Each row = one pyramid level: [IW_x, IW_z]
IntrWndw = [% wide x for large horizontal displacement,
            256 64;
            128 32;
             64 32;    %   shallow z to stay below surface
             32 16];
GrdSpc = [IntrWndw(:,1)/2,   IntrWndw(:,2)/2];   % 50% both, all passes
GrdSpc(end,1) = IntrWndw(end,1)/4;                % 75% x at final pass only
% compVel = ComputeVelocities_Quick_Filt_Deform_Water_dcorFilt( ...
%          IM_a, IM_b, imSurfa.mask, imSurfb.mask, IntrWndw, GrdSpc,1);
dcorGate=0;
%allow inter-pass UOD
iuod.enabled  = true;
iuod.remove   = 2.0;
iuod.reinsert = 3.0;
iuod.minvec   = 5;     % 5 is slightly stricter
compVel = ComputeVelocities_Quick_Filt_Deform_Water_dcorFilt( ...
    IM_a, IM_b, imSurfa.mask, imSurfb.mask, IntrWndw, GrdSpc, dcorGate, iuod);

    compVel.DX = DX;
    compVel.DT = DT;
end


%% Validate
clearvars opts
% Bare minimum — just gap fill, no removal:
% opts.do_uod    = false;
% opts.do_groups = false;
% opts.denoise   = false;
% opts.fill_gaps = true;
% compVel = validatePIV(compVel, opts);
% defaults: UOD + group removal + denoising + gap fill
% compVel = validatePIV(compVel);

% % Step 1 — UOD only (no gap fill between passes, no groups):
% opts.do_uod    = true;
% opts.fill_gaps = true;
% compVel = validatePIV(compVel, opts);
% 
% % Step 2 — add group removal:
opts.do_uod    = true;
opts.do_groups = true;
opts.fill_gaps = true;
compVel = validatePIV(compVel, opts);
% 
% % Step 3 — add smoothing:
% opts.do_uod    = true;
% opts.do_groups = true;
% opts.fill_gaps = true;
% opts.denoise   = true;   % or smoothn externally with tunable S
% compVel = validatePIV(compVel, opts);

% % Step 4 — if still noisy, add dcor gate:
% opts.do_dcor  = true;
% opts.dcor_min = 0.3;
% compVel = validatePIV(compVel, opts);

% remove outlier
% compVel.delta_x = removeOutliers(compVel.delta_x, compVel.dcor);
% % compVel.delta_z = removeOutliers(compVel.delta_z, compVel.dcor);
% optional smoothing pass
compVel.delta_x = smoothn(compVel.delta_x, 0.1, 'robust'); %Fabrice had 0.4
compVel.delta_z = smoothn(compVel.delta_z, 0.1, 'robust');

%% Extract fields
delta_x = compVel.delta_x .* compVel.Mask;   % in pixel per frame
delta_z = compVel.delta_z .* compVel.Mask;
x = compVel.xPIV.*DX; %in m
z = compVel.zPIV.*DX; %in m
[Nz_piv, Nx_piv] = size(delta_x);
Surface_PIV = imSurfa.surfacePIVImg;    % [1 x w] pixel row coords

fprintf('PIV grid: %d x %d   Image: %d x %d\n', Nz_piv, Nx_piv, height, width);

u_raw = delta_x .* compVel.Mask .* DX/DT;   % m/s, lab frame
w_raw = delta_z .* compVel.Mask .* DX/DT;

figure;

set(gcf, 'DefaultAxesColor', [0.9 0.9 0.9])  % all subplots get gray NaN
subplot(1,2,1)
h = imagesc(x, z, delta_x .* compVel.Mask);
set(h, 'AlphaData', ~isnan(u_raw))
colorbar; colormap(gca, brewermap([],'Spectral'))
xlabel('x (m)'); ylabel('z (m)')
title('u\_raw (pix/dt)')
daspect([1,1,1])
clim([-2,20])

subplot(1,2,2)
h = imagesc(x, z, delta_z .* compVel.Mask);
set(h, 'AlphaData', ~isnan(w_raw))
colorbar; colormap(gca, brewermap([],'Spectral'))
xlabel('x (m)'); ylabel('z (m)')
title('w\_raw (pix/dt)')
daspect([1,1,1])
clim([-8,8])
% hold on
% quiver(x(1:skip:end), z(1:skip:end), ...
%        u(1:skip:end, 1:skip:end), ...
%        w(1:skip:end, 1:skip:end), ...
%        2, 'k')
% hold off

set(gcf,'color','white')

%% =========================================================================
%% STEP 2 — COORDINATE TRANSFORM + METHOD A ORBITALS
%% =========================================================================
% Pixel offset: maps surface-image coords to PIV-image coords
% (experiment-specific; from Main_LC.m)
SU_OFFSET = 0;%-1716 + 287;

pivRes.zPIV = compVel.zPIV;
pivRes.xPIV = compVel.xPIV;
pivRes.GS   = compVel.GS;
pivRes.mask  = compVel.Mask;

% no drift (default):
transfo = generateTransfo_LC_noLFV_2023(compVel, Surface_PIV, pivRes);

% % with wind drift (e.g. 0.05 m/s):
% opt.U_drift = 0.05;
% transfo = generateTransfo_LC_noLFV_2023(compVel, Surface_PIV, pivRes, opt);

SU   = transfo.SU(2:end,:);        % remove zeta=0 row (the surface itself)
SU   = SU + SU_OFFSET;             % surface-image → PIV-image pixel coords
ORBX = transfo.ORBX(2:end,:);      % horizontal orbital velocity
ORBZ = transfo.ORBZ(2:end,:);      % vertical orbital velocity

pivRes.pf_surf = SU(1,:);

fprintf('SU size: %d x %d   ORBX size: %d x %d\n', size(SU), size(ORBX));

% --- Plot SU grid lines on Cartesian velocity ---
figure;
imagesc(x,z,u_raw)
colorbar; colormap gray; hold on

% Plot every 20th constant-zeta line
n_su_lines = size(SU, 1);
line_skip = max(1, round(n_su_lines / 15));
for iz = 1:line_skip:n_su_lines
    plot((1:size(SU,2))* compVel.GS*  DX, SU(iz,:) * DX, '-r', 'LineWidth', 0.5)
end
% Surface line
plot((1:length(Surface_PIV))*DX, Surface_PIV*DX, '-r', 'LineWidth', 0.5)

xlabel('x (m)'); ylabel('z (m)')
daspect([1,1,1])
title(sprintf('Constant-\\zeta lines on u (m/s) frame %s', pair_str))
drawnow

%% =========================================================================
%% STEP 3 — TRANSFORM VELOCITY TO WAVE-FOLLOWING
%% =========================================================================
intrp_u_raw = transformVelField_decay_forFab(u_raw, pivRes, SU);
intrp_w_raw = transformVelField_decay_forFab(w_raw, pivRes, SU);

% figure out where nan is at the bottom 
% deepest row where the full width is valid
last_valid_row = find(all(~isnan(intrp_u_raw), 2), 1, 'last');

% rows of intrp_u correspond to altitude(2:end) = compVel.zPIV - GS/2
z_ax = (compVel.zPIV - compVel.GS/2) * DX;  % m

figure;
subplot(1,2,1)
imagesc(x,z_ax,intrp_u_raw)
colorbar; colormap(gca, brewermap([],'Spectral'))
xlabel('x (m)'); ylabel('\zeta (levels)')
title('u\_raw wave-following (m/s)')
clim([-0.01,0.12])
axis equal;axis tight
ylim([z_ax(1), z_ax(last_valid_row)])
xlim([compVel.xPIV(1), compVel.xPIV(end)] * DX)

subplot(1,2,2)
imagesc(x,z_ax,intrp_w_raw)
colorbar; colormap(gca, brewermap([],'Spectral'))
xlabel('x (m)'); ylabel('\zeta (levels)')
title('w\_raw wave-following (m/s)')
clim([-0.04,0.04])
sgtitle(sprintf('Wave-following velocity frame %s', pair_str))
axis equal;axis tight
ylim([z_ax(1), z_ax(last_valid_row)])
xlim([compVel.xPIV(1), compVel.xPIV(end)] * DX)
drawnow

%% =========================================================================
%% STEP 4 — METHOD A ORBITALS
%% =========================================================================
 % Convert pixel/frame -> m/s
    ORBX_ms = ORBX * DX / DT;
    ORBZ_ms = ORBZ * DX / DT;

figure;
subplot(1,2,1)
imagesc(x,z_ax,ORBX * DX/DT)
colorbar; colormap(gca, brewermap([],'Spectral'))
xlabel('x (m)'); ylabel('\zeta (m)')
title('u orbital (m/s)')
axis equal;axis tight
clim([-0.01,0.12])
ylim([z_ax(1), z_ax(last_valid_row)])
xlim([compVel.xPIV(1), compVel.xPIV(end)] * DX)

subplot(1,2,2)
imagesc(x,z_ax,ORBZ * DX/DT)
colorbar; colormap(gca, brewermap([],'Spectral'))
xlabel('x (m)'); ylabel('\zeta (m)')
title('w orbital (m/s)')
axis equal;axis tight
clim([-0.04,0.04])
ylim([z_ax(1), z_ax(last_valid_row)])
xlim([compVel.xPIV(1), compVel.xPIV(end)] * DX)

sgtitle(sprintf('linear wave theory orbital velocities frame %s', pair_str))
drawnow

%% =========================================================================
%% STEP 5 — SUBTRACT ORBITALS aka RESIDUAL
%% =========================================================================
u_res = intrp_u_raw - ORBX_ms;   % mean+turb residual, wave-following
w_res = intrp_w_raw - ORBZ_ms;

figure(6); clf;
subplot(3,2,1)
imagesc(x,z_ax,intrp_u_raw);
colorbar; colormap(gca, brewermap([],'Spectral'))
xlabel('x (m)'); ylabel('\zeta (m)')
title('u\_raw wave-following (m/s)')
axis equal;axis tight
clim([-0.01,0.12])
ylim([z_ax(1), z_ax(last_valid_row)])
xlim([compVel.xPIV(1), compVel.xPIV(end)] * DX)

subplot(3,2,2)
imagesc(x,z_ax,intrp_w_raw);
colorbar; colormap(gca, brewermap([],'Spectral'))
xlabel('x (m)'); ylabel('\zeta (m)')
title('w\_raw wave-following (m/s)')
axis equal;axis tight
clim([-0.04,0.04])
ylim([z_ax(1), z_ax(last_valid_row)])
xlim([compVel.xPIV(1), compVel.xPIV(end)] * DX)

subplot(3,2,3)
imagesc(x,z_ax,ORBX_ms);
colorbar; colormap(gca, brewermap([],'Spectral'))
xlabel('x (m)'); ylabel('\zeta (m)')
title('u orbital (m/s)')
axis equal;axis tight
clim([-0.04,0.04])
ylim([z_ax(1), z_ax(last_valid_row)])
xlim([compVel.xPIV(1), compVel.xPIV(end)] * DX)

subplot(3,2,4)
imagesc(x,z_ax,ORBZ_ms);
colorbar; colormap(gca, brewermap([],'Spectral'))
xlabel('x (m)'); ylabel('\zeta (m)')
title('w orbital (m/s)')
axis equal;axis tight
clim([-0.04,0.04])
ylim([z_ax(1), z_ax(last_valid_row)])
xlim([compVel.xPIV(1), compVel.xPIV(end)] * DX)

subplot(3,2,5)
imagesc(x,z_ax,u_res);
colorbar; colormap(gca, brewermap([],'Spectral'))
xlabel('x (m)'); ylabel('\zeta (m)')
title('u (mean + turb) (m/s)')
axis equal;axis tight
clim([-0.01,0.12])
ylim([z_ax(1), z_ax(last_valid_row)])
xlim([compVel.xPIV(1), compVel.xPIV(end)] * DX)

subplot(3,2,6)
imagesc(x,z_ax,w_res);
colorbar; colormap(gca, brewermap([],'Spectral'))
xlabel('x (m)'); ylabel('\zeta (m)')
title('w (mean + turb) (m/s)')
axis equal;axis tight
clim([-0.04,0.04])
ylim([z_ax(1), z_ax(last_valid_row)])
xlim([compVel.xPIV(1), compVel.xPIV(end)] * DX)

set(gcf,'Color','white')

%% Transform wave following frame back to lab frame
u_res_lab  = reverseTransformVelField_decay_forFab(u_res,       pivRes, SU);
w_res_lab  = reverseTransformVelField_decay_forFab(w_res,       pivRes, SU);
u_orb_lab  = reverseTransformVelField_decay_forFab(ORBX_ms,     pivRes, SU);
w_orb_lab  = reverseTransformVelField_decay_forFab(ORBZ_ms,     pivRes, SU);
u_raw_lab  = reverseTransformVelField_decay_forFab(intrp_u_raw, pivRes, SU);
w_raw_lab  = reverseTransformVelField_decay_forFab(intrp_w_raw, pivRes, SU);

figure(7); clf;
subplot(3,2,1)
imagesc(x,z,u_raw_lab);
colorbar; colormap(gca, brewermap([],'Spectral'))
xlabel('x (m)'); ylabel('z (m)')
title('u\_raw lab (m/s)')
axis equal;axis tight
clim([-0.01,0.12])

subplot(3,2,2)
imagesc(x,z,w_raw_lab);
colorbar; colormap(gca, brewermap([],'Spectral'))
xlabel('x (m)'); ylabel('z (m)')
title('w\_raw lab (m/s)')
axis equal;axis tight
clim([-0.04,0.04])

subplot(3,2,3)
imagesc(x,z,u_orb_lab);
colorbar; colormap(gca, brewermap([],'Spectral'))
xlabel('x (m)'); ylabel('z (m)')
title('u\_orb lab (m/s)')
axis equal;axis tight
clim([-0.04,0.04])

subplot(3,2,4)
imagesc(x,z,w_orb_lab);
colorbar; colormap(gca, brewermap([],'Spectral'))
xlabel('x (m)'); ylabel('z (m)')
title('w\_orb lab (m/s)')
axis equal;axis tight
clim([-0.04,0.04])

subplot(3,2,5)
imagesc(x,z,u_res_lab);
colorbar; colormap(gca, brewermap([],'Spectral'))
xlabel('x (m)'); ylabel('z (m)')
title('u\_res lab — mean+turb (m/s)')
axis equal;axis tight
clim([-0.01,0.12])

subplot(3,2,6)
imagesc(x,z,w_res_lab);
colorbar; colormap(gca, brewermap([],'Spectral'))
xlabel('x (m)'); ylabel('z (m)')
title('w\_res lab — mean+turb (m/s)')
axis equal;axis tight
clim([-0.04,0.04])

set(gcf,'Position',[700,100,600,600])
set(gcf,'Color','white')