% make_decomposition_video.m
% Render a 3x2 lab-frame decomposition video from saved PIVMat_TURB files.
% Mirrors figure(7) of run_decomposition_linearWaveTheory.m.
%
% Run run_decomposition_loop.m first to populate PIVMat_TURB/.

clear; clc;

%% =========================================================================
%% PARAMETERS
%% =========================================================================
LONG     = 'D:\DelawareDataBackup\Longitudinal\PIV\';
ii       = 4;            % experiment index
fps      = 7.2;
DX       = 1/17697.69;   % m/pixel
DT       = 10e-3;        % s/frame

%% =========================================================================
%% PATHS
%% =========================================================================
DIRS      = dir(LONG); DIRS = DIRS(3:end);
exp_name  = DIRS(ii).name;
turb_save = [LONG exp_name '/Chris_recompute/PIVMat_TURB/'];

files = dir([turb_save exp_name '_compVel_*.mat']);
nums  = cellfun(@(s) sscanf(s, [exp_name '_compVel_%d.mat']), {files.name});
[~,o] = sort(nums); files = files(o);
N     = length(files);
fprintf('Experiment: %s   Frames: %d\n', exp_name, N);

out_name = strcat('D:\DelawareDataResult\',sprintf('decomposition_%s.mp4', exp_name));
vw = VideoWriter(out_name, 'MPEG-4');
vw.FrameRate = fps;
open(vw);

%% =========================================================================
%% LOOP
%% =========================================================================
fig = figure('Position',[700,100,800,800],'Color','white');

for ff = 1:N
    S      = load([turb_save files(ff).name], 'decomposedVel', 'pivRes');
    cv     = S.decomposedVel.compVel;
    pivRes = S.pivRes;

    x = pivRes.xPIV * DX;   % m
    z = pivRes.zPIV * DX;   % m

    clf(fig);
    panels = { ...
       cv.u_mean_lab, 'u mean (m/s)',        [-0.01, 0.12]; ...
       cv.w_mean_lab, 'w mean (m/s)',        [-0.04, 0.04]; ...
       cv.u_orb_lab,  'u orbital (m/s)',     [-0.04, 0.04]; ...
       cv.w_orb_lab,  'w orbital (m/s)',     [-0.04, 0.04]; ...
       cv.u_res,      'u (mean+turb) (m/s)', [-0.01, 0.12]; ...
       cv.w_res,      'w (mean+turb) (m/s)', [-0.04, 0.04]};
    for k = 1:6
        subplot(3,2,k);
        imagesc(x, z, panels{k,1});
        colorbar; colormap(gca, brewermap([],'Spectral'));
        xlabel('x (m)'); ylabel('z (m)');
        title(panels{k,2});
        axis equal; axis tight;
        clim(panels{k,3});
    end
    sgtitle(sprintf('%s  frame %d/%d', exp_name, ff, N), 'Interpreter','none');
    drawnow;
    writeVideo(vw, getframe(fig));

    if mod(ff,50) == 0
        fprintf('  rendered %d / %d frames\n', ff, N);
    end
end
close(vw);
fprintf('Wrote %s (%d frames)\n', out_name, N);

