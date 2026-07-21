% test_window_frame157.m
% One-frame A/B test of interrogation-window pyramids, to check the spurious
% ~35 px/frame near-surface displacement seen with the 256-start pyramid.
%
% NON-DESTRUCTIVE: computes compVel in memory for each case, plots dx side by
% side. Saves nothing, overwrites no precomputed compVel files.
%
% Mirrors the data-loading of run_decomposition_linearWaveTheory.m.

clear; clc;

%% ---------------- Params ----------------
LONG     = 'D:\DelawareDataBackup\Longitudinal\PIV\';
rootpath = 'D:\Scripps';
addpath(strcat(rootpath,'\GC-Wave-Gen\M-Files_FabMarcNovDec2014\'));
addpath(strcat(rootpath,'\GC-Wave-Gen\M-Files_FabMarcNovDec2014\FabriceScripts\'));
addpath(strcat(rootpath,'\GC-Wave-Gen\M-Files_FabMarcNovDec2014\CrapperOptimizedFindSurface\'));

ii = 4;                    % experiment index
image_pair_number = 157;   % frame to test
DX = 1/17697.69;           % m per pixel
DT = 10e-3;                % s per pair

% Pyramids to compare:  { label, IntrWndw, GrdSpc }
% (both square — the current ComputeVelocities only accepts scalar windows)
cases = {
    '256-start (current)', [256 128 64 32 16 8], [128 64 32 16 8 4]
    '64-start (test)',     [64 32 16 8],         [32 16 8 4]
    };

% Common color limits so both panels are directly comparable (px/frame)
CLIM = [-5 40];

%% ---------------- Load frame ----------------
DIRS = dir(LONG); DIRS = DIRS(3:end);
exp_name  = DIRS(ii).name;
load_path = [LONG exp_name];
pair_str  = sprintf('%03d', image_pair_number);
fprintf('Experiment: %s   Frame: %s\n\n', exp_name, pair_str);

A = load([load_path '/PIVRaw/PIV/' exp_name '_Piv_' pair_str '_a.mat']); IM_a = A.imgPiv;
B = load([load_path '/PIVRaw/PIV/' exp_name '_Piv_' pair_str '_b.mat']); IM_b = B.imgPiv;

imSurfa = FindSurfaceCapillary( ...
    [load_path '/PIVRaw/PIVSURF/' exp_name '_Pivsurf_' pair_str '_a.mat'], findMask=true);
imSurfb = FindSurfaceCapillary( ...
    [load_path '/PIVRaw/PIVSURF/' exp_name '_Pivsurf_' pair_str '_b.mat'], findMask=true);

surf_x = (1:numel(imSurfa.surfacePIVImg)) * DX;
surf_y = imSurfa.surfacePIVImg * DX;

%% ---------------- Compute + plot each case ----------------
nC = size(cases,1);
figure('Color','w','Position',[60 80 700*nC 640]);

for k = 1:nC
    label = cases{k,1};
    IW    = cases{k,2};
    GS    = cases{k,3};

    fprintf('Computing %-22s  IW=[%s]\n', label, num2str(IW));
    cv = ComputeVelocities_Quick_NoFilt_Deform_Water( ...
        IM_a, IM_b, imSurfa.mask, imSurfb.mask, IW, GS);

    dx = cv.delta_x .* cv.Mask;        % raw masked displacement (px/frame)
    x  = cv.xPIV * DX;
    z  = cv.zPIV * DX;

    v = dx(isfinite(dx));
    fprintf('   grid %dx%d   dx (px/frame): median %.1f   95th pct %.1f   max %.1f\n\n', ...
        size(dx,1), size(dx,2), median(v), prctile(v,95), max(v));

    subplot(1,nC,k)
    imagesc(x, z, dx); axis image; colorbar
    try
        colormap(gca, brewermap([],'Spectral'));
    catch
        colormap(gca, parula);
    end
    clim(CLIM)
    hold on
    plot(surf_x, surf_y, '-k', 'LineWidth', 1.2)
    xlabel('x (m)'); ylabel('z (m)')
    title(sprintf('%s\ndx raw (px/frame)', label), 'Interpreter','none')
end

sgtitle(sprintf('Window-size A/B — %s frame %s  (shared clim [%g %g])', ...
    exp_name, pair_str, CLIM(1), CLIM(2)), 'Interpreter','none')
drawnow
