% generate_pivlab_masks.m
% Generates per-frame binary mask .mat files for PIVlab "Import pixel mask".
% One file per image (A and B) per frame.
%
% Convention: mask = 0 (black) → valid water region
%             mask = 255 (white) → excluded air region above surface
%
% In PIVlab: Import pixel mask, select all files for the frame range,
% starting from the first loaded image pair.
%
% Output files saved to: <chris_path>/PIVlab_masks/
%   <exp_name>_mask_<NNN>_a.mat
%   <exp_name>_mask_<NNN>_b.mat

%% =========================================================================
%% PARAMETERS — edit these
%% =========================================================================
LONG          = 'D:\DelawareDataBackup\Longitudinal\PIV\';
ii            = 4;          % experiment index
frame_start   = 216;        % first frame number
frame_end     = 216;        % last frame number
num_of_digits = 3;          % zero-padding width in filenames
glint_buffer_px = 40;       % rows below surface to exclude (covers glint/foam band)

%% =========================================================================
%% SETUP
%% =========================================================================
DIRS      = dir(LONG);
DIRS      = DIRS(3:end);
exp_name  = DIRS(ii).name;
load_path = [LONG exp_name];

chris_path = [load_path '/Chris_recompute/'];
piv_save   = [chris_path 'PIVMat/'];
surf_path  = [load_path '/PIVRaw/PIVSURF/'];

out_dir = [chris_path 'PIVlab_masks/'];
if ~exist(out_dir, 'dir'), mkdir(out_dir); end

fprintf('Experiment : %s\n', exp_name);
fprintf('Frames     : %d – %d\n', frame_start, frame_end);
fprintf('Output     : %s\n\n', out_dir);

%% =========================================================================
%% LOOP
%% =========================================================================
for frame = frame_start:frame_end
    pair_str = sprintf(['%0' num2str(num_of_digits) 'd'], frame);

    %% --- Load surface detection ---
    compvel_file = [piv_save exp_name '_compVel_' pair_str '.mat'];
    if exist(compvel_file, 'file')
        data    = load(compvel_file, 'imSurfa', 'imSurfb');
        imSurfa = data.imSurfa;
        imSurfb = data.imSurfb;
    else
        surf_a = [surf_path exp_name '_Pivsurf_' pair_str '_a.mat'];
        surf_b = [surf_path exp_name '_Pivsurf_' pair_str '_b.mat'];
        if ~exist(surf_a, 'file') || ~exist(surf_b, 'file')
            fprintf('  [SKIP] frame %s — surface files not found\n', pair_str);
            continue
        end
        imSurfa = FindSurfaceCapillary(surf_a, findMask=true);
        imSurfb = FindSurfaceCapillary(surf_b, findMask=true);
    end

    %% --- Binary exclusion masks (255 = excluded/air, 0 = valid water) ---
    mask_a = uint8(~(imSurfa.mask == 1)) * 255;
    mask_b = uint8(~(imSurfb.mask == 1)) * 255;

    % Extend exclusion zone downward to cover surface glint/foam band
    [h, w] = size(mask_a);
    rows = (1:h)';
    surf_a = round(imSurfa.surfacePIVImg(:)');
    mask_a(  (rows > surf_a) & (rows <= min(surf_a + glint_buffer_px, h))  ) = 255;
    surf_b = round(imSurfb.surfacePIVImg(:)');
    mask_b(  (rows > surf_b) & (rows <= min(surf_b + glint_buffer_px, h))  ) = 255;

    %% --- Save ---
    out_a = [out_dir exp_name '_mask_' pair_str '_a.mat'];
    out_b = [out_dir exp_name '_mask_' pair_str '_b.mat'];
    mask = mask_a; save(out_a, 'mask'); %#ok<NASGU>
    mask = mask_b; save(out_b, 'mask'); %#ok<NASGU>

    fprintf('  frame %s — saved\n', pair_str);
end

fprintf('\nDone.\n');
fprintf('In PIVlab: navigate to frame %d, then Import pixel mask\n', frame_start);
fprintf('and select all _a and _b files in order (a then b per pair).\n');
