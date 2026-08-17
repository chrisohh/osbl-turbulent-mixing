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
rootpath = 'C:\Users\airsealab\Documents\GitHub';
% rootpath = 'D:\Scripps';
addpath(strcat(rootpath,'\GC-Wave-Gen\M-Files_FabMarcNovDec2014\'));
addpath(strcat(rootpath,'\GC-Wave-Gen\M-Files_FabMarcNovDec2014\FabriceScripts\'));
addpath(strcat(rootpath,'\GC-Wave-Gen\M-Files_FabMarcNovDec2014\CrapperOptimizedFindSurface\'));

p = get_piv_params('ExpLCL_2_02');

exp_name      = p.exp_name;
DX            = p.DX;       % 1/17697.69 m/pixel
DT            = p.DT;       % 10e-3 s, inter-frame (A→B), for velocity scaling
num_of_digits = 3;

SU_OFFSET = 0;     % surface-image → PIV-image pixel offset

% recompute_piv = true : rerun ComputeVelocities even if PIVMat file exists
%                        (reuses saved imSurfa/imSurfb, skips FindSurfaceCapillary)
%               false : use cached compVel if available
recompute_piv = false;

% do_decomposition = false : produce only PIVMat/ (raw velocity + dcor) — enough
%                            for make_raw_qc_video.m. Skips wave removal entirely.
%                  true  : also run the setup-specific decomposition -> PIVMat_TURB/
do_decomposition = true;

GLINT_BUFFER = 20;  %290 for transverse and 20 for longitudinal % px below detected surface to exclude (glint/foam band)
fprintf('Glint buffer: %d px = %.2f mm\n', GLINT_BUFFER, GLINT_BUFFER * DX * 1e3);

IntrWndw = [128 128;
             64 64;    %   shallow z to stay below surface
             32 32;
             16 16];
GrdSpc = [IntrWndw(:,1)/2, IntrWndw(:,2)/2];

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
load_path      = p.load_path;
piv_save       = p.piv_save;
turb_save      = p.turb_save;
piv_folder     = p.piv_folder;     % 'PIV' (LCL) or 'PIVCC' (transverse)
pivsurf_folder = p.pivsurf_folder; % 'PIVSURF' (LCL) or 'PIVSURFCC' (transverse)
if ~exist(piv_save,  'dir'), mkdir(piv_save);  end
if ~exist(turb_save, 'dir'), mkdir(turb_save); end

%% =========================================================================
%% LOOP
%% =========================================================================
raw_files = dir([load_path '/PIVRaw/' piv_folder '/' exp_name '_Piv_*_a.mat']);
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
    cf_orig = fullfile(piv_save, [exp_name '_compVel_' ps '.mat']);
    cache_hit = exist(cf_orig, 'file');
    if cache_hit
        tmp     = load(cf_orig, 'compVel', 'imSurfa', 'imSurfb');
        compVel = tmp.compVel;
        imSurfa = tmp.imSurfa;
        imSurfb = tmp.imSurfb;
    end
    if ~cache_hit || recompute_piv
        is_transverse = ismember(p.setup, {'LCTA', 'LCTB'});
        if ~cache_hit
            fprintf('  computing compVel for pair %s...\n', ps);
            if is_transverse
                % Transverse calibration ported from Main_LC_TRAN2.m -- unverified
                % against a real frame, see FindSurfaceCapillaryTransverse.m.
                imSurfa = FindSurfaceCapillaryTransverse( ...
                    [load_path '/PIVRaw/' pivsurf_folder '/' exp_name '_Pivsurf_' ps '_a.mat'], findMask=true);
                imSurfb = FindSurfaceCapillaryTransverse( ...
                    [load_path '/PIVRaw/' pivsurf_folder '/' exp_name '_Pivsurf_' ps '_b.mat'], findMask=true);
            else
                imSurfa = FindSurfaceCapillary( ...
                    [load_path '/PIVRaw/' pivsurf_folder '/' exp_name '_Pivsurf_' ps '_a.mat'], findMask=true);
                imSurfb = FindSurfaceCapillary( ...
                    [load_path '/PIVRaw/' pivsurf_folder '/' exp_name '_Pivsurf_' ps '_b.mat'], findMask=true);
            end
        else
            fprintf('  recomputing PIV for pair %s (reusing saved surfaces)...\n', ps);
        end
        load([load_path '/PIVRaw/' piv_folder '/' exp_name '_Piv_' ps '_a.mat']); IMa = imgPiv;
        load([load_path '/PIVRaw/' piv_folder '/' exp_name '_Piv_' ps '_b.mat']); IMb = imgPiv;
        if is_transverse
            % Rotate+crop to match the pixel grid FindSurfaceCapillaryTransverse's mask assumes.
            IMa = RectifyPIVFrameTransverse(IMa);
            IMb = RectifyPIVFrameTransverse(IMb);
        end

%         % frame-dependent pyramid
%         if pair_num <= 210
%             IntrWndw = [128 32;
%                          64 32;
%                          32 16];
%         else
%             IntrWndw = [256 64;
%                         128 32;
%                          64 32;
%                          32 16];
%         end
%         GrdSpc = [IntrWndw(:,1)/2, IntrWndw(:,2)/2];
%         GrdSpc(end,1) = IntrWndw(end,1)/4;
        Mask_a = apply_glint_buffer(imSurfa, GLINT_BUFFER);
        Mask_b = apply_glint_buffer(imSurfb, GLINT_BUFFER);
        compVel = ComputeVelocities_Quick_Filt_Deform_Water_dcorFilt( ...
            IMa, IMb, Mask_a, Mask_b, IntrWndw, GrdSpc, 0, iuod);
    end
    compVel.DX = DX; compVel.DT = DT;

    Surface_PIV = imSurfa.surfacePIVImg;

    % --- Cartesian velocity: reuse cached u,w if present, else compute ---
    if isfield(compVel, 'u_raw') && isfield(compVel, 'w_raw') && ...
       isfield(compVel, 'u_raw_nan') && isfield(compVel, 'w_raw_nan') && ~recompute_piv
        u_raw     = compVel.u_raw;     w_raw     = compVel.w_raw;
        u_raw_nan = compVel.u_raw_nan; w_raw_nan = compVel.w_raw_nan;
    else
        % Snapshot before gap-fill/smoothn: no post-PIV interpolation or smoothing.
        % NaN preserved at air (Mask) and PIV-gate-failed locations (delta_x already NaN).
        u_raw_nan = compVel.delta_x .* compVel.Mask .* DX/DT;
        w_raw_nan = compVel.delta_z .* compVel.Mask .* DX/DT;
        compVel.u_raw_nan = u_raw_nan;
        compVel.w_raw_nan = w_raw_nan;

        compVel = validatePIV(compVel, val_opts);
        compVel.delta_x = smoothn(compVel.delta_x, 0.1, 'robust');
        compVel.delta_z = smoothn(compVel.delta_z, 0.1, 'robust');
        u_raw = compVel.delta_x .* compVel.Mask .* DX/DT;   % m/s, lab frame
        w_raw = compVel.delta_z .* compVel.Mask .* DX/DT;
        compVel.u_raw = u_raw; compVel.w_raw = w_raw;
    end

    % --- save PIVMat (new or recomputed) — slim imSurfa/imSurfb to mask+surface only ---
    if ~cache_hit || recompute_piv
        imSurfa = struct('mask', imSurfa.mask, 'surfacePIVImg', imSurfa.surfacePIVImg, ...
                         'ImgScaledToPIVSmallCrop', imSurfa.ImgScaledToPIVSmallCrop);
        imSurfb = struct('mask', imSurfb.mask, 'surfacePIVImg', imSurfb.surfacePIVImg, ...
                         'ImgScaledToPIVSmallCrop', imSurfb.ImgScaledToPIVSmallCrop);
        save(fullfile(piv_save, [exp_name '_compVel_' ps '.mat']), 'compVel', 'imSurfa', 'imSurfb', '-v7');
    end

    if do_decomposition
    % --- wave-turbulence decomposition (setup-specific) ---
    switch p.setup
        case 'LCL'    % longitudinal: image x-axis is wave-propagation direction
            [D, pivRes] = decompose_longitudinal(u_raw, w_raw, compVel, Surface_PIV, SU_OFFSET);
        case {'LCTA', 'LCTB'}   % transverse: cross-wave plane, different physics
            [D, pivRes] = decompose_transverse(u_raw, w_raw, compVel, Surface_PIV, SU_OFFSET);
        otherwise
            error('Unknown setup for decomposition: %s', p.setup);
    end

    % --- build and save decomposedVel to PIVMat_TURB ---
    decomposedVel.compVel.u_raw        = single(D.u_raw);        % measured, lab frame (gap-filled + smoothed)
    decomposedVel.compVel.w_raw        = single(D.w_raw);
    decomposedVel.compVel.u_raw_nan    = single(u_raw_nan);      % raw PIV output, NaN at bad-dcor/air (no gap-fill, no smoothn)
    decomposedVel.compVel.w_raw_nan    = single(w_raw_nan);
    decomposedVel.compVel.intrp_u_raw  = single(D.intrp_u_raw);  % measured, wave-following
    decomposedVel.compVel.intrp_w_raw  = single(D.intrp_w_raw);
    decomposedVel.compVel.ORBX_ms      = single(D.ORBX_ms);      % orbital, wave-following
    decomposedVel.compVel.ORBZ_ms      = single(D.ORBZ_ms);
    decomposedVel.compVel.intrp_u_res  = single(D.intrp_u_res);  % mean+turb, wave-following
    decomposedVel.compVel.intrp_w_res  = single(D.intrp_w_res);
    decomposedVel.compVel.u_res_lab    = single(D.u_res_lab);    % mean+turb, lab frame
    decomposedVel.compVel.w_res_lab    = single(D.w_res_lab);
    decomposedVel.compVel.u_orb_lab    = single(D.u_orb_lab);    % orbital, lab frame
    decomposedVel.compVel.w_orb_lab    = single(D.w_orb_lab);
    decomposedVel.compVel.SU           = single(D.SU);
    decomposedVel.compVel.pf_surf      = Surface_PIV;
    decomposedVel.compVel.dcor         = single(compVel.dcor);  % correlation quality (NaN=air; use dcor<0.4 to mask before Reynolds stresses)

    save(fullfile(turb_save, [exp_name '_compVel_' ps '.mat']), 'decomposedVel', 'pivRes');
    end   % do_decomposition

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
