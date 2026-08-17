% run_decomposition_transverse_singleframe.m
% Single-frame sanity check of the TRANSVERSE wave removal (decompose_transverse).
% Mirrors run_decomposition_linearWaveTheory.m but for the cross-wind (y,z) plane:
% there is no wave-following transform, ORBX/ORBZ are uniform in y, and the
% residual is the turbulent fluctuation (y-mean method) or mean+turb (eta method).
%
% Implements OSBL-notes/wave_removal.tex. Run on one frame before the batch loop.

clear; clc;

%% =========================================================================
%% PARAMETERS
%% =========================================================================
rootpath = 'C:\Users\airsealab\Documents\GitHub';
addpath(strcat(rootpath,'\GC-Wave-Gen\M-Files_FabMarcNovDec2014\'));
addpath(strcat(rootpath,'\GC-Wave-Gen\M-Files_FabMarcNovDec2014\FabriceScripts\'));
addpath(strcat(rootpath,'\GC-Wave-Gen\M-Files_FabMarcNovDec2014\CrapperOptimizedFindSurface\'));

p = get_piv_params('ExpLCTB_2_01');     % a transverse run (LCTA/LCTB)

exp_name          = p.exp_name;
DX                = p.DX;
DT                = p.DT;
image_pair_number = 216;
num_of_digits     = 3;

method = 'ymean';   % 'ymean' (self-contained) or 'eta' (needs longitudinal eta)
rerunsurf = 1;
% NOTE: 290 was tuned against the (wrong) longitudinal calibration as a
% vertical-registration fudge, not a real glint/foam band width -- see
% chat history. Now that FindSurfaceCapillaryTransverse uses transverse-
% specific calibration, re-check this against a real frame; it should
% likely shrink toward something like the longitudinal ~20-40 px.
GLINT_BUFFER = 30;
IntrWndw = [64 64;
            32 32;
            16 16];
GrdSpc = [IntrWndw(:,1)/2, IntrWndw(:,2)/2];
iuod.enabled = true; iuod.remove = 2.0; iuod.reinsert = 3.0; iuod.minvec = 5;
val_opts.do_uod = true; val_opts.do_groups = true; val_opts.fill_gaps = true;

%% =========================================================================
%% LOAD OR COMPUTE compVel  (same pipeline as the loop)
%% =========================================================================
load_path      = p.load_path;
piv_save       = p.piv_save;
piv_folder     = p.piv_folder;
pivsurf_folder = p.pivsurf_folder;
ps = sprintf(['%0' num2str(num_of_digits) 'd'], image_pair_number);
fprintf('Experiment: %s   Pair: %s   method: %s\n', exp_name, ps, method);

cf = fullfile(piv_save, [exp_name '_compVel_' ps '.mat']);
if exist(cf, 'file') && rerunsurf==0
    fprintf('Loading cached compVel from %s\n', cf);
    load(cf, 'compVel', 'imSurfa');
else
    fprintf('Computing compVel from raw images...\n');
    % Transverse calibration (dewarp, resize, crop, vertical offset) ported from
    % Main_LC_TRAN2.m -- see FindSurfaceCapillaryTransverse.m for provenance/caveats.
    % Unlike the longitudinal FindSurfaceCapillary, these numbers are unverified;
    % use measure_registration_offset_transverse.m to check/refine them.
    imSurfa = FindSurfaceCapillaryTransverse( ...
        [load_path '/PIVRaw/' pivsurf_folder '/' exp_name '_Pivsurf_' ps '_a.mat'], findMask=true);
    imSurfb = FindSurfaceCapillaryTransverse( ...
        [load_path '/PIVRaw/' pivsurf_folder '/' exp_name '_Pivsurf_' ps '_b.mat'], findMask=true);
    load([load_path '/PIVRaw/' piv_folder '/' exp_name '_Piv_' ps '_a.mat']); IMa = imgPiv;
    load([load_path '/PIVRaw/' piv_folder '/' exp_name '_Piv_' ps '_b.mat']); IMb = imgPiv;
    % Rotate+crop to match the pixel grid FindSurfaceCapillaryTransverse's mask assumes.
    IMa = RectifyPIVFrameTransverse(IMa);
    IMb = RectifyPIVFrameTransverse(IMb);

    %% --- Surface-detection sanity plots ---
    % Raw surface camera
    figure;
    imagesc(imSurfa.ImgScaledToPIVSmallCrop,[0,300]); hold on;
    plot(imSurfa.surface_raw, '-r')
    colormap bone
    axis tight; axis equal
    title(sprintf('Surface Camera: %s frame %s', exp_name, ps), ...
        'Interpreter', 'none')
    set(gcf,'Color','white')

    % Raw PIV camera
    figure;
    subplot(1,2,1)
    imagesc(IMa,[0,500]);
    hold on; plot(imSurfa.surfacePIVImg, '-r', 'LineWidth', 1)
    colormap gray
    axis tight; axis equal
    title(sprintf('PIV %s frame %s - A', exp_name, ps), ...
        'Interpreter', 'none')

    subplot(1,2,2)
    imagesc(IMb,[0,500]);
    hold on; plot(imSurfb.surfacePIVImg, '-r', 'LineWidth', 1)
    colormap gray
    axis tight; axis equal
    title(sprintf('PIV %s frame %s - B', exp_name, ps), ...
        'Interpreter', 'none')

    set(gcf,'Color','white')
    drawnow

    Mask_a = apply_glint_buffer(imSurfa, GLINT_BUFFER);
    Mask_b = apply_glint_buffer(imSurfb, GLINT_BUFFER);

    %% --- Mask overlay sanity plot ---
    figure;
    imagesc(IMa, [0,500]);
    colormap(gca, gray)
    hold on;
    red_overlay = cat(3, ones(size(Mask_a)), zeros(size(Mask_a)), zeros(size(Mask_a)));
    h = imagesc(red_overlay);
    h.AlphaData = 0.3 * isnan(Mask_a);
    axis tight; axis equal
    title(sprintf('PIV %s frame %s - A (mask overlay)', exp_name, ps), ...
        'Interpreter', 'none')
    set(gcf,'Color','white')
    drawnow

    compVel = ComputeVelocities_Quick_Filt_Deform_Water_dcorFilt( ...
        IMa, IMb, Mask_a, Mask_b, IntrWndw, GrdSpc, 0, iuod);
    compVel = validatePIV(compVel, val_opts);
    compVel.delta_x = smoothn(compVel.delta_x, 0.1, 'robust');
    compVel.delta_z = smoothn(compVel.delta_z, 0.1, 'robust');
end
compVel.DX = DX; compVel.DT = DT;

Surface_PIV = imSurfa.surfacePIVImg;
u_raw = compVel.delta_x .* compVel.Mask .* DX/DT;   % m/s, lab frame (cross-wind v)
w_raw = compVel.delta_z .* compVel.Mask .* DX/DT;   % m/s, lab frame (vertical w)

%% =========================================================================
%% TRANSVERSE DECOMPOSITION
%% =========================================================================
opt = struct('method', method);
% For method='eta' also supply the matched longitudinal surface line, e.g.:
%   opt.eta    = eta_x_at_this_time;   % metres, sampled along x
%   opt.eta_dx = eta_dx;               % metres
[D, pivRes] = decompose_transverse(u_raw, w_raw, compVel, Surface_PIV, 0, opt);


%% =========================================================================
%% PLOTS
%% =========================================================================
x = compVel.xPIV .* DX * 1e3;    % cross-wind y (mm)
z = compVel.zPIV .* DX * 1e3;    % depth (mm)
surf_z = interp1(1:numel(imSurfa.surfacePIVImg), imSurfa.surfacePIVImg, compVel.xPIV, 'linear', 'extrap') * DX* 1e3;

figure('Color','w','Position',[100 100 900 700]);
set(gcf,'DefaultAxesColor',[0.9 0.9 0.9]);
ttl = {'u_{raw} (y,z)','w_{raw} (y,z)', ...
       'ORBX (u wave/mean, uniform in y)','ORBZ (w wave, uniform in y)', ...
       'u_{res} = u_{raw} - ORBX','w_{res} = w_{raw} - ORBZ'};
flds = {u_raw, w_raw, D.ORBX_ms, D.ORBZ_ms, D.intrp_u_res, D.intrp_w_res};
for k = 1:6
    subplot(3,2,k)
    h = imagesc(x, z, flds{k}); set(h,'AlphaData',~isnan(flds{k}));
    hold on; plot(x, surf_z, 'k-', 'LineWidth', 1.5);
    colorbar; colormap(gca, brewermap([],'Spectral'));
    xlabel('y (mm)'); ylabel('depth (mm)'); title(ttl{k});
    axis tight; daspect([1 1 1]);
    clim([-0.04,0.04])
end
sgtitle(sprintf('Transverse decomposition (%s) — %s frame %s', method, exp_name, ps), ...
        'Interpreter','none');
