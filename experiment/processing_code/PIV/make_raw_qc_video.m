% make_raw_qc_video.m
% Lightweight QC video: u_raw, w_raw, dcor — loaded from PIVMat/.
% Run after run_decomposition_loop.m (or even mid-run once files appear).

clear; clc;

%% =========================================================================
%% PARAMETERS
%% =========================================================================
p = get_piv_params('ExpLCL_2_01');

exp_name = p.exp_name;
piv_save = p.piv_save;
DX       = p.DX;            % 1/17697.69 m/pixel
DT       = p.DT;            % 10e-3 s, inter-frame (A→B), for velocity scaling
Fs_cam   = p.Fs_cam;        % 14.4 Hz — camera acquisition rate (images)
Fs_PIV   = p.Fs_PIV;        % 7.2 Hz  — velocity map rate (one pair = two images)
fps      = Fs_PIV;          % video frame rate matches acquisition

clim_u = [-0.01  0.15];   % m/s  — adjust to your flow
clim_w = [-0.05  0.05];   % m/s

%% =========================================================================
%% PATHS
%% =========================================================================

files = dir([piv_save exp_name '_compVel_*.mat']);
nums  = cellfun(@(s) sscanf(s, [exp_name '_compVel_%d.mat']), {files.name});
[~,o] = sort(nums); files = files(o);
N     = length(files);
fprintf('Experiment: %s   Frames: %d\n', exp_name, N);

out_name = strcat('D:\DelawareDataResult\', sprintf('raw_qc_%s.mp4', exp_name));
vw = VideoWriter(out_name, 'MPEG-4');
vw.FrameRate = fps;
open(vw);

%% bwr colormap (blue=negative, white=zero, red=positive)
n   = 256;
r   = [linspace(0,1,n/2), ones(1,n/2)];
g   = [linspace(0,1,n/2), linspace(1,0,n/2)];
b   = [ones(1,n/2), linspace(1,0,n/2)];
bwr = [r; g; b]';

%% =========================================================================
%% LOOP
%% =========================================================================
fig = figure('Position',[100,100,900,700],'Color','white');
set(fig, 'DefaultAxesColor', [0.85 0.85 0.85]);   % gray for NaN/masked

for ff = 1:N
    S      = load([piv_save files(ff).name], 'compVel', 'imSurfa');
    cv     = S.compVel;
    surf   = S.imSurfa;

    x = cv.xPIV * DX * 1e3;   % mm
    z = cv.zPIV * DX * 1e3;   % mm

    u = cv.u_raw;
    w = cv.w_raw;
    d = cv.dcor;

    % True surface line in mm (interpolated to PIV x-grid)
    surf_z = interp1(1:numel(surf.surfacePIVImg), surf.surfacePIVImg, cv.xPIV, 'linear', 'extrap') * DX * 1e3;

    clf(fig);

    % --- u_raw ---
    ax1 = subplot(2,2,1);
    h = imagesc(x, z, u);
    set(h, 'AlphaData', ~isnan(u));
    colormap(ax1, brewermap([], 'Spectral')); clim(clim_u); colorbar;
    hold on; plot(x, surf_z, 'k-', 'LineWidth', 1.5);
    xlabel('x (mm)'); ylabel('z (mm)');
    title('u_{raw} (m/s)'); axis tight;
    clim([-0.01, 0.12])
axis equal;axis tight
    % --- w_raw ---
    ax2 = subplot(2,2,2);
    h = imagesc(x, z, w);
    set(h, 'AlphaData', ~isnan(w));
    colormap(ax2, brewermap([], 'Spectral')); clim(clim_w); colorbar;
    hold on; plot(x, surf_z, 'k-', 'LineWidth', 1.5);
    xlabel('x (mm)'); ylabel('z (mm)');
    title('w_{raw} (m/s)'); axis tight;
clim([-0.04, 0.04])
axis equal;axis tight
    % --- dcor (full bottom row) ---
    ax3 = subplot(2,2,[3 4]);
    h = imagesc(x, z, d);
    set(h, 'AlphaData', ~isnan(d));
    colormap(ax3, brewermap([], 'Spectral')); clim([0 1]); colorbar;
    hold on; plot(x, surf_z, 'k-', 'LineWidth', 1.5);
    xlabel('x (mm)'); ylabel('z (mm)');
    title('dcor  (0=no peak, 1=sharp peak)'); axis tight;
axis equal;axis tight
    sgtitle(sprintf('%s  frame %d/%d', exp_name, ff-1, N-1), 'Interpreter','none');
    drawnow;
    writeVideo(vw, getframe(fig));

    if mod(ff, 50) == 0
        fprintf('  rendered %d / %d frames\n', ff, N);
    end
end

close(vw);
fprintf('Wrote %s (%d frames)\n', out_name, N);
