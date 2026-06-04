% run_decomposition_loop.m
% Batch loop: compute Method A (LWT) wave-turbulence decomposition for all
% frames. Saves per-frame snapshots and accumulates ensemble average.
%
% Run run_decomposition_linearWaveTheory.m first to sanity-check a single
% frame before running this.

clear; clc;

%% =========================================================================
%% PARAMETERS  (keep in sync with run_decomposition_linearWaveTheory.m)
%% =========================================================================
LONG = 'D:\DelawareDataBackup\Longitudinal\PIV\';

rootpath = 'C:\Users\airsealab\Documents\GitHub';
% rootpath = 'D:\Scripps';
addpath(strcat(rootpath,'\GC-Wave-Gen\M-Files_FabMarcNovDec2014\'));
addpath(strcat(rootpath,'\GC-Wave-Gen\M-Files_FabMarcNovDec2014\FabriceScripts\'));
addpath(strcat(rootpath,'\GC-Wave-Gen\M-Files_FabMarcNovDec2014\CrapperOptimizedFindSurface\'));

ii            = 4;
num_of_digits = 3;

DX = 1/17697.69;   % m per pixel
DT = 10e-3;        % sec per image pair

SU_OFFSET = 0;     % surface-image → PIV-image pixel offset

% recompute_piv = true : rerun ComputeVelocities even if PIVMat file exists
%                        (reuses saved imSurfa/imSurfb, skips FindSurfaceCapillary)
%               false : use cached compVel if available
recompute_piv = true;

% inter-pass UOD settings (same for all frames)
iuod.enabled  = true;
iuod.remove   = 2.0;
iuod.reinsert = 3.0;
iuod.minvec   = 5;

% validatePIV settings
val_opts.do_uod    = true;
val_opts.do_groups = true;
val_opts.fill_gaps = true;

%% =========================================================================
%% PATHS
%% =========================================================================
DIRS = dir(LONG);
DIRS = DIRS(3:end);
exp_name  = DIRS(ii).name;
load_path = [LONG exp_name];

chris_path = [load_path '/Chris_recompute/'];
piv_save   = [chris_path 'PIVMat/'];
turb_save  = [chris_path 'PIVMat_TURB/'];
if ~exist(piv_save,  'dir'), mkdir(piv_save);  end
if ~exist(turb_save, 'dir'), mkdir(turb_save); end

%% =========================================================================
%% LOOP
%% =========================================================================
raw_files = dir([load_path '/PIVRaw/PIV/' exp_name '_Piv_*_a.mat']);
N_frames  = length(raw_files);
fprintf('Experiment: %s   Frames: %d\n', exp_name, N_frames);

ensembleSum_u   = [];   % initialised from first frame
ensembleSum_w   = [];
ensembleStats_u = [];
ensembleStats_w = [];
N_processed = 0;

for ff = 1:N_frames
    tok = regexp(raw_files(ff).name, '_Piv_(\d+)_a\.mat$', 'tokens');
    pair_num = str2double(tok{1}{1});
    ps = sprintf(['%0' num2str(num_of_digits) 'd'], pair_num);

    % --- load or compute compVel + surfaces ---
    cf_orig = [piv_save exp_name '_compVel_' ps '.mat'];
    cache_hit = exist(cf_orig, 'file');
    if cache_hit
        tmp     = load(cf_orig, 'compVel', 'imSurfa', 'imSurfb');
        compVel = tmp.compVel;
        imSurfa = tmp.imSurfa;
        imSurfb = tmp.imSurfb;
    end
    if ~cache_hit || recompute_piv
        if ~cache_hit
            fprintf('  computing compVel for pair %s...\n', ps);
            imSurfa = FindSurfaceCapillary( ...
                [load_path '/PIVRaw/PIVSURF/' exp_name '_Pivsurf_' ps '_a.mat'], findMask=true);
            imSurfb = FindSurfaceCapillary( ...
                [load_path '/PIVRaw/PIVSURF/' exp_name '_Pivsurf_' ps '_b.mat'], findMask=true);
        else
            fprintf('  recomputing PIV for pair %s (reusing saved surfaces)...\n', ps);
        end
        load([load_path '/PIVRaw/PIV/' exp_name '_Piv_' ps '_a.mat']); IMa = imgPiv;
        load([load_path '/PIVRaw/PIV/' exp_name '_Piv_' ps '_b.mat']); IMb = imgPiv;

        % frame-dependent pyramid
        if pair_num <= 210
            IntrWndw = [128 32;
                         64 32;
                         32 16];
        else
            IntrWndw = [256 64;
                        128 32;
                         64 32;
                         32 16];
        end
        GrdSpc = [IntrWndw(:,1)/2, IntrWndw(:,2)/2];
        GrdSpc(end,1) = IntrWndw(end,1)/4;

        compVel = ComputeVelocities_Quick_Filt_Deform_Water_dcorFilt( ...
            IMa, IMb, imSurfa.mask, imSurfb.mask, IntrWndw, GrdSpc, 0, iuod);
    end
    compVel.DX = DX; compVel.DT = DT;

    Surface_PIV = imSurfa.surfacePIVImg;

    % --- Cartesian velocity: reuse cached u,w if present, else compute ---
    if isfield(compVel, 'u_raw') && isfield(compVel, 'w_raw') && ~recompute_piv
        u_raw = compVel.u_raw; w_raw = compVel.w_raw;
    else
        compVel = validatePIV(compVel, val_opts);
        compVel.delta_x = smoothn(compVel.delta_x, 0.1, 'robust');
        compVel.delta_z = smoothn(compVel.delta_z, 0.1, 'robust');
        u_raw = compVel.delta_x .* compVel.Mask .* DX/DT;   % m/s, lab frame
        w_raw = compVel.delta_z .* compVel.Mask .* DX/DT;
        compVel.u_raw = u_raw; compVel.w_raw = w_raw;
    end

    % --- save PIVMat only if newly computed ---
    if ~cache_hit
        save([piv_save exp_name '_compVel_' ps '.mat'], 'compVel', 'imSurfa', 'imSurfb');
    end

    % --- wave-following transform ---
    pivRes.zPIV = compVel.zPIV; pivRes.xPIV = compVel.xPIV;
    pivRes.GS   = compVel.GS;   pivRes.mask  = compVel.Mask;
    transfo     = generateTransfo_LC_noLFV_2023(compVel, Surface_PIV, pivRes);
    SU          = transfo.SU(2:end,:)   + SU_OFFSET;
    ORBX        = transfo.ORBX(2:end,:);    % pixels/frame
    ORBZ        = transfo.ORBZ(2:end,:);
    ORBX_ms     = ORBX * DX/DT;             % m/s
    ORBZ_ms     = ORBZ * DX/DT;
    pivRes.pf_surf = SU(1,:);

    % transform raw velocity to wave-following frame
    intrp_u_raw = transformVelField_decay_forFab(u_raw, pivRes, SU);   % m/s, wave-following
    intrp_w_raw = transformVelField_decay_forFab(w_raw, pivRes, SU);

    % subtract orbitals → mean+turb residual in wave-following frame
    intrp_u_res = intrp_u_raw - ORBX_ms;
    intrp_w_res = intrp_w_raw - ORBZ_ms;

    % reverse-transform to lab frame
    u_res_lab  = reverseTransformVelField_decay_forFab(intrp_u_res, pivRes, SU);
    w_res_lab  = reverseTransformVelField_decay_forFab(intrp_w_res, pivRes, SU);
    u_orb_lab  = reverseTransformVelField_decay_forFab(ORBX_ms,     pivRes, SU);
    w_orb_lab  = reverseTransformVelField_decay_forFab(ORBZ_ms,     pivRes, SU);

    % --- build and save decomposedVel to PIVMat_TURB ---
    decomposedVel.compVel.u_raw        = single(u_raw);        % measured, lab frame
    decomposedVel.compVel.w_raw        = single(w_raw);
    decomposedVel.compVel.intrp_u_raw  = single(intrp_u_raw);  % measured, wave-following
    decomposedVel.compVel.intrp_w_raw  = single(intrp_w_raw);
    decomposedVel.compVel.ORBX_ms      = single(ORBX_ms);      % orbital, wave-following
    decomposedVel.compVel.ORBZ_ms      = single(ORBZ_ms);
    decomposedVel.compVel.intrp_u_res  = single(intrp_u_res);  % mean+turb, wave-following
    decomposedVel.compVel.intrp_w_res  = single(intrp_w_res);
    decomposedVel.compVel.u_res_lab    = single(u_res_lab);    % mean+turb, lab frame
    decomposedVel.compVel.w_res_lab    = single(w_res_lab);
    decomposedVel.compVel.u_orb_lab    = single(u_orb_lab);    % orbital, lab frame
    decomposedVel.compVel.w_orb_lab    = single(w_orb_lab);
    decomposedVel.compVel.SU           = single(SU);
    decomposedVel.compVel.pf_surf      = Surface_PIV;

    save([turb_save exp_name '_compVel_' ps '.mat'], 'decomposedVel', 'pivRes');

%     % --- accumulate 2D ensemble average ---
%     if isempty(ensembleSum_u)
%         ensembleSum_u   = zeros(size(intrp_u_res));
%         ensembleSum_w   = zeros(size(intrp_w_res));
%         ensembleStats_u = zeros(size(intrp_u_res));
%         ensembleStats_w = zeros(size(intrp_w_res));
%     end
%     valid_u = isfinite(intrp_u_res);
%     valid_w = isfinite(intrp_w_res);
%     ensembleSum_u(valid_u)   = ensembleSum_u(valid_u)   + intrp_u_res(valid_u);
%     ensembleSum_w(valid_w)   = ensembleSum_w(valid_w)   + intrp_w_res(valid_w);
%     ensembleStats_u          = ensembleStats_u + valid_u;
%     ensembleStats_w          = ensembleStats_w + valid_w;

    N_processed = N_processed + 1;

    if mod(ff, 50) == 0
        fprintf('  processed %d / %d frames\n', ff, N_frames);
    end
end

% %% =========================================================================
% %% ENSEMBLE MEAN
% %% =========================================================================
% u_mean_wf      = ensembleSum_u ./ ensembleStats_u;   % [Nz x Nx] time-mean, wave-following, m/s
% w_mean_wf      = ensembleSum_w ./ ensembleStats_w;
% u_mean_profile = nanmean(u_mean_wf, 2);               % [Nz x 1] depth profile
% w_mean_profile = nanmean(w_mean_wf, 2);
% 
% fprintf('Ensemble average complete over %d frames.\n', N_processed);
% 
% % derive plot axes from last processed frame
% x    = compVel.xPIV .* DX * 1e3;                      % mm
% z_ax = (compVel.zPIV - compVel.GS/2) .* DX * 1e3;    % mm
% last_valid_row = find(all(~isnan(u_mean_wf), 2), 1, 'last');
% 
% figure;
% subplot(2,1,1)
% imagesc(x, z_ax, u_mean_wf)
% colorbar; colormap(gca, brewermap([], 'Spectral'))
% xlabel('x (mm)'); ylabel('depth (mm)')
% title(sprintf('u mean — %d frames (m/s)', N_processed))
% ylim([z_ax(1), z_ax(last_valid_row)])
% xlim([x(1), x(end)])
% clim([-0.01, 0.12])
% 
% subplot(2,1,2)
% imagesc(x, z_ax, w_mean_wf)
% colorbar; colormap(gca, brewermap([], 'Spectral'))
% xlabel('x (mm)'); ylabel('depth (mm)')
% title(sprintf('w mean — %d frames (m/s)', N_processed))
% ylim([z_ax(1), z_ax(last_valid_row)])
% xlim([x(1), x(end)])
% clim([-0.04, 0.04])
% drawnow
% 
% fprintf('\nDone.\n  PIVMat      → %s\n  PIVMat_TURB → %s\n', piv_save, turb_save);
