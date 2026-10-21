%% CISG slope video: raw | Sx | Sy | eta, one frame per video frame
% Video counterpart of the cisg_main.m preview (same setup via util/cisg_setup_run,
% same figure via util/cisg_plot_frame, slopes from cisg_calculate_slopes, eta from
% cisg_slope_to_eta).  Per-frame slopes are cached in the folder / format of the
% cisg_main cache (<out_root>\CISG_slopes_CoreView_<N>\frame_NNNN.mat, ROI-cropped
% single Sx, Sy), so re-rendering is cheap and the cache is shared.  bin > 1 uses a
% separate folder (..._bin<B>).  Colour limits are fixed over the video (given, or
% +-3 rms of the first rendered frame).
clear; close all;
worker_gb = 1.5;                    % estimated peak RAM per parallel worker at bin 2 (GB); adjust to your measurement
here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, 'util'), fullfile(here, 'calibration'), fullfile(fileparts(here), 'IR Camera', 'calibration'));

%% ---------------- User settings ----------------
cfg.coreview_num = 95;
cfg.frames       = 3000:1:3009;    % start:stride:end
cfg.roi_mode     = 'pick';         % 'pick' (first time; then the saved pick is reused), 'default' (central 60%) or a roi struct
cfg.roi_repick   = false;          % true = pick the ROI again
cfg.out_root     = 'D:\HLAB_2026\SlopeGauge';
cfg.out_path     = '';             % '' = <out_root>/CISG_slopes_CoreView_<N>_frames_<a>-<b>.mp4
cfg.out_fps      = 10;
cfg.quality      = 90;
cfg.cache_slopes = true;           % save / reuse per-frame slopes
cfg.lim_s        = [];             % slope colour limit (+-); [] = 3 * rms of the first frame, then held
cfg.lim_eta      = [];             % eta limit (mm, +-);       [] = likewise
cfg.bin          = 2;              % block-average factor before the slopes (1 = full resolution; 2 matched bin 1 within 0.5% rms on frame 3000)
cfg.fig_size     = [1280 800];     % px, fixed for the whole video
cfg.use_parfor   = true;           % compute missing slopes in parallel (needs cfg.cache_slopes)
cfg.n_workers    = 'auto';         % 'auto' = min(cores - 4, 0.6 * free RAM / worker_gb, 32), or a number
cfg.show_live    = true;           % show the frame being computed / written (falls back to false without a desktop, e.g. matlab -batch)

%% ---------------- Setup ----------------
R = cisg_setup_run(struct('coreview_num', cfg.coreview_num, 'roi_mode', cfg.roi_mode, ...
                          'roi_repick', cfg.roi_repick));
[ref_b, setup_b] = cisg_bin_setup(R.ref_image, R.setup, cfg.bin);
frames = cfg.frames;
live = cfg.show_live && usejava('desktop');
if cfg.show_live && ~live, disp('No desktop: show_live switched off'); end
cache_dir = '';
if cfg.cache_slopes
    sub = sprintf('CISG_slopes_CoreView_%d', cfg.coreview_num);
    if cfg.bin > 1, sub = sprintf('%s_bin%d', sub, cfg.bin); end
    % Cached slopes are ROI-cropped, so a non-default ROI gets its own folder
    % (the default central-60% cache keeps its old name)
    r = R.setup.roi;
    if ~R.setup.roi_is_default
        sub = sprintf('%s_roi_r%d-%d_c%d-%d', sub, r.row_min, r.row_max, r.col_min, r.col_max);
    end
    cache_dir = fullfile(cfg.out_root, sub);
    if ~exist(cache_dir, 'dir'), mkdir(cache_dir); end
end
if isempty(cfg.out_path)
    cfg.out_path = fullfile(cfg.out_root, sprintf('CISG_slopes_CoreView_%d_frames_%d-%d.mp4', ...
                            cfg.coreview_num, frames(1), frames(end)));
end
if ~exist(fileparts(cfg.out_path), 'dir'), mkdir(fileparts(cfg.out_path)); end
fprintf('CoreView %d: %d frames (%d to %d), ROI rows %d-%d cols %d-%d, bin %d\n', cfg.coreview_num, numel(frames), ...
        frames(1), frames(end), R.setup.roi.row_min, R.setup.roi.row_max, R.setup.roi.col_min, R.setup.roi.col_max, cfg.bin);
%% ---------------- Slopes (parallel pre-pass, optional) ----------------
if cfg.use_parfor
    if ~cfg.cache_slopes, error('cisg_slope_video:parfor', 'use_parfor needs cache_slopes = true.'); end
    have_pct = license('test', 'Distrib_Computing_Toolbox') && ~isempty(ver('parallel'));
    todo = frames(~arrayfun(@(n) isfile(fullfile(cache_dir, sprintf('frame_%04d.mat', n))), frames));
    if have_pct && ~isempty(todo)
        nw = cfg.n_workers;
        [~, sv] = memory;  ram_gb = sv.PhysicalMemory.Available / 1e9;
        if ischar(nw)                                        % 'auto': leave 4 cores and 40% of the free RAM
            nw = max(1, min([feature('numcores') - 4, floor(0.6 * ram_gb / worker_gb), 32]));
        end
        nw = min(nw, numel(todo));
        fprintf('Parallel: %d workers (cores %d, free RAM %.0f GB, %.1f GB per worker assumed), %d frames to compute\n', ...
                nw, feature('numcores'), ram_gb, worker_gb, numel(todo));
        pool = gcp('nocreate');
        if ~isempty(pool) && pool.NumWorkers ~= nw, delete(pool); pool = []; end
        t_pool = tic;
        if isempty(pool), parpool('local', nw); end
        paths = {fullfile(here, 'util'), fullfile(here, 'calibration'), fullfile(fileparts(here), 'IR Camera', 'calibration')};
        pctRunOnAll(sprintf('addpath(''%s'', ''%s'', ''%s'')', paths{:}));
        fprintf('pool ready: %.1f s\n', toc(t_pool));
        t0 = tic;
        ff = R.frame_file;  st = setup_b;  rb = ref_b;  bn = cfg.bin;  dtt = R.dt;  cd_ = cache_dir;
        dq = [];  fp = [];
        if live                                              % workers send back a thumbnail of each finished frame
            dq = parallel.pool.DataQueue;
            fp = figure('Name', 'computing slopes', 'Position', [50 50 700 560]);
            axes(fp);  fp.UserData = 0;                      %#ok<LAXES>  (frame counter)
            afterEach(dq, @(d) show_progress(fp, d, numel(todo), t0));
        end
        parfor i = 1:numel(todo)
            [sx, ~] = cisg_cache_frame(todo(i), rb, st, ff, cd_, bn, dtt);
            if ~isempty(dq), send(dq, struct('n', todo(i), 'Sx', sx(1:8:end, 1:8:end))); end
        end
        if live, close(fp); end
        fprintf('parfor: %d frames in %.1f s (%.1f s/frame wall)\n', numel(todo), toc(t0), toc(t0) / numel(todo));
    elseif ~have_pct
        fprintf('Parallel Computing Toolbox not available: computing serially in the render loop\n');
    end
end
%% ---------------- Single frame: just the figure, no video ----------------
if isscalar(frames)
    n = frames;
    [Sx, Sy] = cisg_cache_frame(n, ref_b, setup_b, R.frame_file, cache_dir, cfg.bin, R.dt);
    obs = cisg_load_coreview(R.frame_file(n));
    if cfg.bin > 1, obs = uint8(imresize(obs, 1 / cfg.bin, 'box')); end
    fig = figure('Name', sprintf('CISG CoreView %d frame %d', cfg.coreview_num, n), ...
                 'Position', [50 50 cfg.fig_size]);
    opts = struct('bin', cfg.bin, 'lim_s', cfg.lim_s, 'lim_eta', cfg.lim_eta, ...
                  'head', sprintf('CoreView %d, frame %d, t = %.2f s', cfg.coreview_num, n, (n - 1) * R.dt));
    info = cisg_plot_frame(fig, obs, double(Sx), double(Sy), setup_b, R.cal, opts);
    fprintf('frame %d: rms Sx %.3f Sy %.3f eta %.2f mm, masked %.1f%% (single frame: no video written)\n', ...
            n, info.rms_sx, info.rms_sy, info.rms_eta, info.masked);
    return
end

%% ---------------- Render ----------------
fig = figure('Visible', matlab.lang.OnOffSwitchState(live), 'Position', [50 50 cfg.fig_size]);   % visible: the frame being written
vw = VideoWriter(cfg.out_path, 'MPEG-4');
vw.FrameRate = cfg.out_fps;  vw.Quality = cfg.quality;
open(vw);
opts = struct('bin', cfg.bin, 'lim_s', cfg.lim_s, 'lim_eta', cfg.lim_eta);
tic;
for i = 1:numel(frames)
    n = frames(i);
    t1 = tic;
    [Sx, Sy] = cisg_cache_frame(n, ref_b, setup_b, R.frame_file, cache_dir, cfg.bin, R.dt);
    obs = cisg_load_coreview(R.frame_file(n));
    if cfg.bin > 1, obs = uint8(imresize(obs, 1 / cfg.bin, 'box')); end
    opts.head = sprintf('CoreView %d, frame %d, t = %.2f s', cfg.coreview_num, n, (n - 1) * R.dt);
    info = cisg_plot_frame(fig, obs, double(Sx), double(Sy), setup_b, R.cal, opts);
    if i == 1                                                % hold the colour limits from here on
        if isempty(opts.lim_s),   opts.lim_s   = info.lim_s;   fprintf('slope colour limit +-%.3f (first frame, held)\n', opts.lim_s); end
        if isempty(opts.lim_eta), opts.lim_eta = info.lim_eta; fprintf('eta colour limit +-%.2f mm (first frame, held)\n', opts.lim_eta); end
    end
    fig.Position(3:4) = cfg.fig_size;
    if live, drawnow; end
    writeVideo(vw, print(fig, '-RGBImage', '-r96'));            % getframe fails on an invisible figure
    el = toc;
    fprintf('  frame %d (%d/%d): %.1f s, rms Sx %.3f Sy %.3f eta %.2f mm, masked %.1f%% | elapsed %.0f s, ETA %.0f s\n', ...
            n, i, numel(frames), toc(t1), info.rms_sx, info.rms_sy, info.rms_eta, info.masked, el, el / i * (numel(frames) - i));
end
close(vw);  close(fig);
d = dir(cfg.out_path);
fprintf('Wrote %s (%.1f MB), %.1f s total\n', cfg.out_path, d.bytes / 1e6, toc);

function show_progress(fp, d, total, t0)
    % DataQueue callback (client): thumbnail of the slopes a worker just finished
    fp.UserData = fp.UserData + 1;
    k = fp.UserData;
    ax = fp.CurrentAxes;
    imagesc(ax, d.Sx);  axis(ax, 'image');  colormap(ax, 'parula');
    el = toc(t0);
    title(ax, sprintf('computing slopes: frame %d (%d/%d), ETA %.0f s', d.n, k, total, el / k * (total - k)));
    drawnow limitrate;
end
