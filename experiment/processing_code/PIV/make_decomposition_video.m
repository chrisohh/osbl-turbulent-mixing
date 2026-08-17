% make_decomposition_video.m
% Render a 3x2 lab-frame decomposition video from saved PIVMat_TURB files.
% Mirrors figure(7) of run_decomposition_linearWaveTheory.m.
%
% Run run_decomposition_loop.m first to populate PIVMat_TURB/.

clear; clc;

%% =========================================================================
%% PARAMETERS
%% =========================================================================
p = get_piv_params('ExpLCTB_2_01');

exp_name  = p.exp_name;
turb_save = p.turb_save;
DX        = p.DX;           % 1/17697.69 m/pixel
DT        = p.DT;           % 10e-3 s, inter-frame (A→B), for velocity scaling
Fs_cam    = p.Fs_cam;       % 14.4 Hz — camera acquisition rate (images)
Fs_PIV    = p.Fs_PIV;       % 7.2 Hz  — velocity map rate (one pair = two images)
fps       = Fs_PIV;         % video frame rate matches acquisition

makeVideo = 1;          % rerun video
frame     = 157;        % frame to plot

%% =========================================================================
%% PATHS
%% =========================================================================

files = dir(fullfile(turb_save, [exp_name '_compVel_*.mat']));
nums  = cellfun(@(s) sscanf(s, [exp_name '_compVel_%d.mat']), {files.name});
[~,o] = sort(nums); files = files(o);
N     = length(files);
fprintf('Experiment: %s   Frames: %d\n', exp_name, N);

if makeVideo == 1
    clearvars vw
    out_name = strcat('D:\DelawareDataResult\',sprintf('decomposition_%s.mp4', exp_name));
    vw = VideoWriter(out_name, 'MPEG-4');
    vw.FrameRate = fps;
    open(vw);
end
%% =========================================================================
%% LOOP
%% =========================================================================
fig = figure('Position',[700,100,800,800],'Color','white');

set(gcf, 'DefaultAxesColor', [0.9 0.9 0.9])  % all subplots get gray NaN

if makeVideo == 1
    frame_list = 1:N;
else
    frame_list = frame + 1;   % +1: title shows ff-1, so file index = frame+1
end

for ff = frame_list

    S      = load(fullfile(turb_save, files(ff).name), 'decomposedVel', 'pivRes');
    cv     = S.decomposedVel.compVel;
    pivRes = S.pivRes;
    surf   = cv.pf_surf;

    x = pivRes.xPIV * DX;   % m
    z = pivRes.zPIV * DX;   % m

    % True surface line in mm (interpolated to PIV x-grid)
    surf_z = interp1(1:numel(surf), surf, pivRes.xPIV, 'linear', 'extrap') * DX;


    clf(fig);
    panels = { ...
       cv.u_raw,  'u_{raw} (m/s)',            [-0.01, 0.12]; ...
       cv.w_raw,  'w_{raw} (m/s)',            [-0.04, 0.04]; ...
       cv.u_orb_lab,  'u_{orb} (m/s)',            [-0.04, 0.04]; ...
       cv.w_orb_lab,  'w_{orb} (m/s)',            [-0.04, 0.04]; ...
       cv.u_res_lab,  'u_{raw}-u_{orb} (m/s)',    [-0.01, 0.12]; ...
       cv.w_res_lab,  'w_{raw}-w_{orb} (m/s)',    [-0.04, 0.04]};
    for k = 1:6
        subplot(3,2,k);
        h=imagesc(x, z, panels{k,1});
        set(h, 'AlphaData', ~isnan(panels{k,1}))
        colorbar; colormap(gca, brewermap([],'Spectral'));
        hold on; plot(x, surf_z, 'k-', 'LineWidth', 1.5);
        xlabel('x (m)'); ylabel('z (m)');
        title(panels{k,2});
        axis equal; axis tight;
        clim(panels{k,3});
    end
    sgtitle(sprintf('%s  frame %d/%d', exp_name, ff-1, N-1), 'Interpreter','none');
    drawnow;

    if makeVideo == 1
        writeVideo(vw, getframe(fig));
        if mod(ff,50) == 0
            fprintf('  rendered %d / %d frames\n', ff, N);
        end
    end
end

if makeVideo == 1
    close(vw);
    fprintf('Wrote %s (%d frames)\n', out_name, N);
end

