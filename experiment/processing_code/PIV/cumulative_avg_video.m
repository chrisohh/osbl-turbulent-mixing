% cumulative_avg_video.m
% Cumulative-average video of u_raw, w_raw — loaded from PIVMat/.
% Frame ff shows the running mean over frames 1..ff
% (e.g. frame 5 = mean of 1-5, frame 10 = mean of 1-10).
% Run after run_decomposition_loop.m (or even mid-run once files appear).

clear; clc;

%% =========================================================================
%% PARAMETERS
%% =========================================================================
p = get_piv_params('ExpLCL_2_01');   % match the run processed by run_decomposition_loop.m

exp_name = p.exp_name;
piv_save = p.piv_save;
DX       = p.DX;            % 1/17697.69 m/pixel
DT       = p.DT;            % 10e-3 s, inter-frame (A→B), for velocity scaling
Fs_cam   = p.Fs_cam;        % 14.4 Hz — camera acquisition rate (images)
Fs_PIV   = p.Fs_PIV;        % 7.2 Hz  — velocity map rate (one pair = two images)
fps      = Fs_PIV;          % video frame rate matches acquisition

clim_u = [-0.01  0.12];   % m/s  — adjust to your flow
clim_w = [-0.04  0.04];   % m/s

%% =========================================================================
%% PATHS
%% =========================================================================

files = dir(fullfile(piv_save, [exp_name '_compVel_*.mat']));
nums  = cellfun(@(s) sscanf(s, [exp_name '_compVel_%d.mat']), {files.name});
[~,o] = sort(nums); files = files(o);
N     = length(files);
fprintf('Experiment: %s   Frames: %d\n', exp_name, N);

out_name = strcat('D:\DelawareDataResult\', sprintf('cumulative_avg_%s.mp4', exp_name));
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
fig = figure('Position',[100,100,900,400],'Color','white');
set(fig, 'DefaultAxesColor', [0.85 0.85 0.85]);   % gray for NaN/masked

% Running accumulators (NaN-aware): sum of values and count of valid samples.
u_sum = []; u_cnt = [];
w_sum = []; w_cnt = [];
surf_sum = []; surf_cnt = [];   % cumulative mean surface line (mm)

ff_start = 1;   % first frame to start cumulating from

for ff = ff_start:N
    S      = load(fullfile(piv_save, files(ff).name), 'compVel', 'imSurfa');
    cv     = S.compVel;
%     surf   = S.imSurfa;

    x = cv.xPIV * DX * 1e3;   % mm
    z = cv.zPIV * DX * 1e3;   % mm

    % --- accumulate cumulative sums / counts (ignoring NaNs) ---
    u_this = cv.u_raw;
    w_this = cv.w_raw;
    if ff == ff_start
        u_sum = zeros(size(u_this)); u_cnt = zeros(size(u_this));
        w_sum = zeros(size(w_this)); w_cnt = zeros(size(w_this));
    end
    u_valid = ~isnan(u_this);
    w_valid = ~isnan(w_this);
    u_sum(u_valid) = u_sum(u_valid) + u_this(u_valid);
    w_sum(w_valid) = w_sum(w_valid) + w_this(w_valid);
    u_cnt = u_cnt + u_valid;
    w_cnt = w_cnt + w_valid;

    % cumulative mean over frames 1..ff (NaN where no valid samples yet)
    u = u_sum ./ u_cnt;   u(u_cnt == 0) = NaN;
    w = w_sum ./ w_cnt;   w(w_cnt == 0) = NaN;

%     % Instantaneous surface line in mm (interpolated to PIV x-grid)
%     surf_z_this = interp1(1:numel(surf.surfacePIVImg), surf.surfacePIVImg, cv.xPIV, 'linear', 'extrap') * DX * 1e3;
% 
%     % Cumulative-mean surface line (matches the time-averaged velocity fields)
%     if ff == 1
%         surf_sum = zeros(size(surf_z_this)); surf_cnt = zeros(size(surf_z_this));
%     end
%     surf_valid = ~isnan(surf_z_this);
%     surf_sum(surf_valid) = surf_sum(surf_valid) + surf_z_this(surf_valid);
%     surf_cnt = surf_cnt + surf_valid;
%     surf_z = surf_sum ./ surf_cnt;   surf_z(surf_cnt == 0) = NaN;

    clf(fig);

    % --- u_raw ---
    ax1 = subplot(1,2,1);
    h = imagesc(x, z, u);
    set(h, 'AlphaData', ~isnan(u));
    colormap(ax1, brewermap([], 'Spectral')); clim(clim_u); colorbar;
%     hold on; plot(x, surf_z, 'k-', 'LineWidth', 1.5);
    xlabel('x (mm)'); ylabel('z (mm)');
    title('cumulative mean u_{raw} (m/s)'); axis tight;
%     clim([-0.01, 0.12])
axis equal;axis tight
    % --- w_raw ---
    ax2 = subplot(1,2,2);
    h = imagesc(x, z, w);
    set(h, 'AlphaData', ~isnan(w));
    colormap(ax2, brewermap([], 'Spectral')); clim(clim_w); colorbar;
%     hold on; plot(x, surf_z, 'k-', 'LineWidth', 1.5);
    xlabel('x (mm)'); ylabel('z (mm)');
    title('cumulative mean w_{raw} (m/s)'); axis tight;
% clim([-0.04, 0.04])
axis equal;axis tight
  
    sgtitle(sprintf('%s  cumulative avg of frames %d-%d/%d', exp_name, ff_start, ff, N), 'Interpreter','none');
    drawnow;
    writeVideo(vw, getframe(fig));

    if mod(ff, 50) == 0
        fprintf('  rendered %d / %d frames\n', ff, N);
    end
end

close(vw);
fprintf('Wrote %s (%d frames)\n', out_name, N);
