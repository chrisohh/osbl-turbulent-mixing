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
image_pair_number = 216;    % sanity-check frame

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

% --- Raw PIV image ---
load([load_path '/PIVRaw/PIV/' exp_name '_Piv_' pair_str '_a.mat']);
IM_a = imgPiv;
load([load_path '/PIVRaw/PIV/' exp_name '_Piv_' pair_str '_b.mat']);
IM_b = imgPiv;
[height, width] = size(IM_a);

% --- Surface detection ---
imSurfa = FindSurfaceCapillary( ...
    [load_path '/PIVRaw/PIVSURF/' exp_name '_Pivsurf_' pair_str '_a.mat'], ...
    findMask = true);
imSurfb = FindSurfaceCapillary( ...
    [load_path '/PIVRaw/PIVSURF/' exp_name '_Pivsurf_' pair_str '_b.mat'], ...
    findMask = true);

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
IntrWndw = [256 128 64 32 16 8];
GrdSpc   = [128 64 32 16 8 4];

chris_path = [load_path '/Chris_recompute/'];
piv_save   = [chris_path 'PIVMat/'];
turb_save  = [chris_path 'PIVMat_TURB/'];

compvel_file = [piv_save exp_name '_compVel_' pair_str '.mat'];
if exist(compvel_file, 'file')
    fprintf('Loading compVel from %s\n', compvel_file);
    load(compvel_file);
else
    fprintf('compVel not found — computing from raw images...\n');
    compVel = ComputeVelocities_Quick_NoFilt_Deform_Water( ...
        IM_a, IM_b, imSurfa.mask, imSurfb.mask, IntrWndw, GrdSpc);
    compVel.DX = DX;
    compVel.DT = DT;
end

%% Extract fields
delta_x = compVel.delta_x .* compVel.Mask;   % in pixel per frame
delta_z = compVel.delta_z .* compVel.Mask;
x = compVel.xPIV.*DX; %in m
z = compVel.zPIV.*DX; %in m
[Nz_piv, Nx_piv] = size(delta_x);
Surface_PIV = imSurfa.surfacePIVImg;    % [1 x w] pixel row coords

fprintf('PIV grid: %d x %d   Image: %d x %d\n', Nz_piv, Nx_piv, height, width);

% remove outlier
delta_x_remove = removeOutliers(delta_x, compVel.dcor);
delta_x_smooth = smoothn(delta_x_remove,0.4,'robust');
delta_z_remove = removeOutliers(delta_z, compVel.dcor);
delta_z_smooth = smoothn(delta_z_remove,0.4,'robust');

u=delta_x_smooth.* compVel.Mask.*DX/DT;
w=delta_z_smooth.* compVel.Mask.*DX/DT;

figure;
subplot(1,2,1)
imagesc(x, z, u)
colorbar; colormap(gca, brewermap([],'Spectral'))
xlabel('x (m)'); ylabel('z (m)')
title('u (m/s)')
daspect([1,1,1])
clim([-0.01,0.12])

subplot(1,2,2)
imagesc(x, z, w)
colorbar; colormap(gca, brewermap([],'Spectral'))
xlabel('x (m)'); ylabel('z (m)')
title('w (m/s)')
daspect([1,1,1])
clim([-0.04,0.04])


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
imagesc(x,z,u)
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
intrp_u = transformVelField_decay_forFab(u, pivRes, SU);
intrp_w = transformVelField_decay_forFab(w, pivRes, SU);

% figure out where nan is at the bottom 
% deepest row where the full width is valid
last_valid_row = find(all(~isnan(intrp_u), 2), 1, 'last');

% rows of intrp_u correspond to altitude(2:end) = compVel.zPIV - GS/2
z_ax = (compVel.zPIV - compVel.GS/2) * DX;  % m

figure;
subplot(1,2,1)
imagesc(x,z_ax,intrp_u)
colorbar; colormap(gca, brewermap([],'Spectral'))
xlabel('x (m)'); ylabel('\zeta (levels)')
title('u in wave-following (m/s)')
clim([-0.01,0.12])
axis equal;axis tight
ylim([z_ax(1), z_ax(last_valid_row)])
xlim([compVel.xPIV(1), compVel.xPIV(end)] * DX)

subplot(1,2,2)
imagesc(x,z_ax,intrp_w)
colorbar; colormap(gca, brewermap([],'Spectral'))
xlabel('x (m)'); ylabel('\zeta (levels)')
title('w in wave-following (m/s)')
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
u_res = intrp_u - ORBX* DX/DT;
w_res = intrp_w - ORBZ* DX/DT;

figure(6); clf;
subplot(3,2,1)
imagesc(x,z_ax,intrp_u);
colorbar; colormap(gca, brewermap([],'Spectral'))
xlabel('x (m)'); ylabel('\zeta (m)')
title('u mean (m/s)')
axis equal;axis tight
clim([-0.01,0.12])
ylim([z_ax(1), z_ax(last_valid_row)])
xlim([compVel.xPIV(1), compVel.xPIV(end)] * DX)

subplot(3,2,2)
imagesc(x,z_ax,intrp_w);
colorbar; colormap(gca, brewermap([],'Spectral'))
xlabel('x (m)'); ylabel('\zeta (m)')
title('w mean (m/s)')
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
u_res_wf = reverseTransformVelField_decay_forFab( u_res, pivRes, SU );
w_res_wf = reverseTransformVelField_decay_forFab( w_res, pivRes, SU );
u_orb_wf = reverseTransformVelField_decay_forFab( ORBX_ms, pivRes, SU );
w_orb_wf = reverseTransformVelField_decay_forFab( ORBZ_ms, pivRes, SU );
u_mean_wf = reverseTransformVelField_decay_forFab( intrp_u, pivRes, SU );
w_mean_wf = reverseTransformVelField_decay_forFab( intrp_w, pivRes, SU );

figure(7); clf;
subplot(3,2,1)
imagesc(x,z,u_mean_wf);
colorbar; colormap(gca, brewermap([],'Spectral'))
xlabel('x (m)'); ylabel('z (m)')
title('u mean (m/s)')
axis equal;axis tight
clim([-0.01,0.12])

subplot(3,2,2)
imagesc(x,z,w_mean_wf);
colorbar; colormap(gca, brewermap([],'Spectral'))
xlabel('x (m)'); ylabel('z (m)')
title('w mean (m/s)')
axis equal;axis tight
clim([-0.04,0.04])


subplot(3,2,3)
imagesc(x,z,u_orb_wf);
colorbar; colormap(gca, brewermap([],'Spectral'))
xlabel('x (m)'); ylabel('z (m)')
title('u orbital (m/s)')
axis equal;axis tight
clim([-0.04,0.04])


subplot(3,2,4)
imagesc(x,z,w_orb_wf);
colorbar; colormap(gca, brewermap([],'Spectral'))
xlabel('x (m)'); ylabel('z (m)')
title('w orbital (m/s)')
axis equal;axis tight
clim([-0.04,0.04])


subplot(3,2,5)
imagesc(x,z,u_res_wf);
colorbar; colormap(gca, brewermap([],'Spectral'))
xlabel('x (m)'); ylabel('z (m)')
title('u (mean + turb) (m/s)')
axis equal;axis tight
clim([-0.01,0.12])


subplot(3,2,6)
imagesc(x,z,w_res_wf);
colorbar; colormap(gca, brewermap([],'Spectral'))
xlabel('x (m)'); ylabel('z (m)')
title('w (mean + turb) (m/s)')
axis equal;axis tight
clim([-0.04,0.04])

set(gcf,'Position',[700,100,600,600])
set(gcf,'Color','white')