%% Longitudinal view - cam 1 optical-depth video
% Writes an MP4 of camera 1 (Flare 12M125 MCL C0100800) in optical-depth
% mode, -log10(I/I_ref), with fixed colour limits.  No spatial calibration.
%
% Frames are colour-mapped directly to RGB and written to the video (no
% figure rendering), so this stays fast over thousands of frames.
clear; close all;
addpath(fullfile(fileparts(mfilename('fullpath')), 'util'));

%% ---------------- User settings ----------------
coreview_num = 89;
raw_root     = 'D:\HLAB_2026\Core';
cam          = 'Flare 12M125 MCL C0100800';

ref_frame    = 1;           % dye-free background frame
frame_start  = 1;
frame_step   = 5;           % every Nth acquired frame
frame_end    = [];          % [] = all frames found in the folder

fs           = 50;          % Hz acquisition rate
clim_od      = [0 0.05];    % colour limits on -log10(I/I_ref)
cmap_name    = 'bone';
downsample   = 2;           % spatial decimation (2 -> 2048 x 1536)
flip_vert    = false;       % true = flip so row 1 ends up at the bottom

out_file     = fullfile(raw_root, sprintf('CoreView_%d_cam1_od.mp4', coreview_num));
video_fps    = 10;
video_quality = 95;

%% ---------------- Set up ----------------
cam_dir = fullfile(raw_root, sprintf('CoreView_%d', coreview_num), cam);
fname = @(n) fullfile(cam_dir, sprintf('CoreView_%d_%s_%04d.raw', ...
                                       coreview_num, cam, n));

listing = dir(fullfile(cam_dir, sprintf('CoreView_%d_%s_*.raw', coreview_num, cam)));
nAvail = numel(listing);
if nAvail == 0
    error('No raw frames found in %s', cam_dir);
end
if isempty(frame_end), frame_end = nAvail; end
frames = frame_start:frame_step:frame_end;
nFrames = numel(frames);

fprintf('Camera dir : %s\n', cam_dir);
fprintf('Frames     : %d to %d step %d  (%d frames of %d available)\n', ...
        frame_start, frame_end, frame_step, nFrames, nAvail);
fprintf('clim        : [%g %g]   colormap: %s   downsample: %d\n', ...
        clim_od(1), clim_od(2), cmap_name, downsample);
fprintf('Output      : %s\n\n', out_file);

Iref = lif_load_raw(fname(ref_frame));
Iref = Iref(1:downsample:end, 1:downsample:end);
Iref = max(Iref, 1);

cmap = feval(cmap_name, 256);

v = VideoWriter(out_file, 'MPEG-4');
v.FrameRate = video_fps;
v.Quality   = video_quality;
open(v);
cleanupObj = onCleanup(@() close(v));

%% ---------------- Write frames ----------------
tic;
for i = 1:nFrames
    n = frames(i);
    f = fname(n);
    if ~exist(f, 'file')
        warning('Missing frame %d, skipping.', n);
        continue;
    end

    I = lif_load_raw(f);
    I = I(1:downsample:end, 1:downsample:end);

    OD = -log10(max(I, 1) ./ Iref);
    if flip_vert, OD = flipud(OD); end

    % Map to the colormap with fixed limits
    idx = round((OD - clim_od(1)) / (clim_od(2) - clim_od(1)) * 255) + 1;
    idx = min(max(idx, 1), 256);
    rgb = ind2rgb(uint8(idx - 1), cmap);

    writeVideo(v, im2uint8(rgb));

    if mod(i, 25) == 0 || i == nFrames
        el = toc;
        fprintf('  %d/%d (%.1f%%)  t = %.2f s  --  %.2f fps  ETA %.0f s\n', ...
                i, nFrames, 100*i/nFrames, (n-1)/fs, i/el, (nFrames-i)/(i/el));
    end
end

clear cleanupObj;   % closes the writer
fprintf('\nWrote %s  (%d frames, %.1f s)\n', out_file, nFrames, toc);
