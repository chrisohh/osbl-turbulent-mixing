%% CISG slope video: raw | Sx | Sy | eta, one frame per video frame
% Video counterpart of the cisg_main.m preview (same setup via util/cisg_setup_run,
% same figure via util/cisg_plot_frame, slopes from cisg_calculate_slopes, eta from
% cisg_slope_to_eta).  Per-frame slopes are cached in the folder / format of the
% cisg_main cache (<out_root>\CISG_slopes_CoreView_<N>\frame_NNNN.mat, ROI-cropped
% single Sx, Sy), so re-rendering is cheap and the cache is shared.  bin > 1 uses a
% separate folder (..._bin<B>).  Colour limits are fixed over the video (given, or
% +-3 rms of the first rendered frame).
clearvars -except cfg_override; close all;   % cfg_override: set by cisg_pipeline.m, replaces the user settings below
worker_gb = 1.5;                    % estimated peak RAM per parallel worker at bin 2 (GB); adjust to your measurement
here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, 'util'));
cisg_fix_path(here);                 % drop stale copies under .claude/worktrees, put Slope_Gauge first
w_ = which('cisg_load_coreview');
if ~startsWith(w_, here), error('cisg_slope_video:path', 'cisg_load_coreview resolves to %s, not the copy next to this script', w_); end

%% ---------------- User settings ----------------
cfg.coreview_num = 101;
cfg.frames       = 1000:10:3000;    % start:stride:end
cfg.roi_mode     = 'auto';         % 'pick' (first time; then the saved pick is reused), 'default' (central 60%), 'auto' (ROI from the reference: avoids wall shadow / flat colour map) or a roi struct
cfg.auto_margin  = 0;              % 'auto': extra inset (px) so shifted pixels still land on usable colours (100 is safe)
cfg.dye_mask     = false;          % mask cross-stream R6G dye stripes (row-wise blue/brightness test); changes the slopes, so it gets its own cache folder
cfg.dye_model    = '';             % '' = by run (cisg_run_defaults: 'joint' for the R6G runs 95, 96, else 'none'), or 'none' | 'joint'; joint also saves the dye map c4 in the cache
cfg.roi_repick   = false;          % true = pick the ROI again
cfg.out_root     = 'D:\HLAB_2026\SlopeGauge';
cfg.out_path     = '';             % '' = <out_root>/CISG_slopes_CoreView_<N>_frames_<a>-<b>.mp4
cfg.out_fps      = 5;
cfg.quality      = 90;
cfg.cache_slopes = true;           % save / reuse per-frame slopes
cfg.lim_s        = [];             % slope colour limit (+-); [] = 3 * rms of the first frame, then held
cfg.lim_eta      = [];             % eta limit (mm, +-);       [] = likewise
cfg.bin          = 2;              % block-average factor before the slopes (1 = full resolution; 2 matched bin 1 within 0.5% rms on frame 3000)
cfg.fig_size     = [1280 800];     % px, fixed for the whole video
cfg.use_parfor   = true;           % compute missing slopes in parallel (needs cfg.cache_slopes)
cfg.n_workers    = 'auto';         % 'auto' = min(cores - 4, 0.6 * free RAM / worker_gb, 32), or a number
cfg.show_live    = true;           % show the frame being computed / written (falls back to false without a desktop, e.g. matlab -batch)
cfg.bad_frame    = 'blank';        % unreadable raw file: 'blank' (empty panels, title says MISSING), 'hold' (repeat the previous good frame) or 'skip' (no video frame)
cfg.warm_start   = false;          % true: start the colour solver from the previous frame's shift too (and from 0; per pixel the better one); own cache folder (_warm).
                                   % Off: on run 101 frame 2000 it cut the masked area only 35.6% -> 33.9% for ~3x the compute (cisg_pipeline 5c)
cfg.show_series  = true;           % wind (hot-wire, or fan command) and wave-gauge time series under the maps, red line = frame (cisg_wind_wave_series)
cfg.series_xlim  = [];             % time axis of those panels (s since wind start); [] = whole record
cfg.warm_lead    = 3;             % frames solved (not saved) before a frame whose previous frame is not available (stride > 1, start of a parallel block)
cfg.ref_frames   = 2:4:58;         % still-water reference = median of these frames (1 = the old single frame 1, which is ~0.8 counts off on run 101)
cfg.bias_frames  = 60:10:200;      % flat-water frames whose mean slope (fixed bias) is subtracted from every frame; [] = off (cisg_flat_bias)
cfg.drift_corr   = 'none';         % per-frame brightness correction before the inversion: 'none' | 'offset' | 'gain' (cisg_colour_displacement)
cfg.colour_channels = 'RG';        % 'RG' | 'RGB' | 'hybrid' (R, G; R, G, B only where R, G fails or folds). Blue: dye-free runs only (cisg_pipeline 5f / 5h / 5i)
if exist('cfg_override', 'var')
    for f_ = fieldnames(cfg_override).', cfg.(f_{1}) = cfg_override.(f_{1}); end
    clear cfg_override f_
end

%% ---------------- Setup ----------------
R = cisg_setup_run(struct('coreview_num', cfg.coreview_num, 'roi_mode', cfg.roi_mode, ...
                          'roi_repick', cfg.roi_repick, 'auto_margin', cfg.auto_margin, 'dye_model', cfg.dye_model, 'setup', struct('dye_mask', cfg.dye_mask, 'drift_corr', cfg.drift_corr, 'colour_channels', cfg.colour_channels), ...
                          'ref_frames', cfg.ref_frames, 'bin', cfg.bin));
[ref_b, setup_b] = cisg_bin_setup(R.ref_image, R.setup, cfg.bin);
setup_b = cisg_flat_bias(R, ref_b, setup_b, cfg.bin, cfg.bias_frames);   % fixed flat-water bias, subtracted in cisg_cache_frame
% frame time in the titles: since wind (fan) start when the run has a wind clock, plus since camera start
if isfinite(R.cam_delay_s)
    tstr = @(n) sprintf('t = %.2f s after wind start (camera %.2f s)', R.t_wind(n), (n - 1) * R.dt);
else
    tstr = @(n) sprintf('t = %.2f s after camera start (no wind clock)', (n - 1) * R.dt);
end
ts = [];
if cfg.show_series && isfinite(R.cam_delay_s)
    ts = cisg_wind_wave_series(cfg.coreview_num);
    if isempty(ts.wind) && isempty(ts.wave), ts = []; end
    if ~isempty(ts), cfg.fig_size = max(cfg.fig_size, [1600 1000]); end   % raw | Sx | Sy across, series below
end
frames = cfg.frames;
bad_list = [];  bad_msgs = {};                       % frames whose raw file could not be read
live = cfg.show_live && usejava('desktop');
if cfg.show_live && ~live, disp('No desktop: show_live switched off'); end
cache_dir = '';
if cfg.cache_slopes
    sub = sprintf('CISG_slopes_CoreView_%d', cfg.coreview_num);
    if cfg.bin > 1, sub = sprintf('%s_bin%d', sub, cfg.bin); end
    if strcmp(R.setup.dye_model, 'joint'), sub = [sub '_joint2']; end       % 2-channel and joint slopes never share a cache
    sub = [sub '_sigma'];                                                   % v3 masks: physical precision (sigma_S), despeckle, dark mask
    sub = [sub '_yup'];                                                     % Sy = d eta/dy with y up (sign fixed 2026-10-25): old caches not reused
    if cfg.warm_start, sub = [sub '_warm']; end                             % previous-frame warm start changes the slopes
    if ~isequal(cfg.ref_frames, 1)                                          % median reference (frame 1 alone keeps the old name)
        sub = sprintf('%s_ref%d-%d', sub, min(cfg.ref_frames), max(cfg.ref_frames));
    end
    if ~strcmp(cfg.drift_corr, 'none'), sub = [sub '_' cfg.drift_corr]; end  % per-frame brightness correction
    if strcmpi(cfg.colour_channels, 'RGB'), sub = [sub '_rgb']; end         % blue in the colour solve
    if strcmpi(cfg.colour_channels, 'hybrid'), sub = [sub '_hyb']; end      % blue only where R, G fails / folds
    if ~isempty(cfg.bias_frames), sub = [sub '_flat']; end                  % flat-water bias subtracted
    if cfg.dye_mask, sub = [sub '_dye']; end
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
        % same path fix on the workers (a pool started earlier may hold old function versions or a shadowing copy)
        pctRunOnAll(sprintf('addpath(''%s''); cisg_fix_path(''%s'');', fullfile(here, 'util'), here));
        wf = parfevalOnAll(@() which('cisg_load_coreview'), 1);  w_all = cellstr(fetchOutputs(wf, 'UniformOutput', false));
        if ~all(startsWith(w_all, here)), error('cisg_slope_video:path', 'A worker resolves cisg_load_coreview to %s', w_all{find(~startsWith(w_all, here), 1)}); end
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
        fs_ = R.img_size;
        % contiguous blocks of frames, solved in order inside a block so each frame can start from the one before
        % (warm start); only a block's first frame needs the lead-in.  2 blocks per worker for load balance.
        nb = min(numel(todo), 2 * nw);  edges = round(linspace(0, numel(todo), nb + 1));
        blocks = arrayfun(@(k) todo(edges(k) + 1:edges(k + 1)), 1:nb, 'UniformOutput', false);
        warm0 = [];  if cfg.warm_start, warm0 = struct('prev', [], 'lead', cfg.warm_lead); end
        bf = cell(nb, 1);  bm = cell(nb, 1);  bi = cell(nb, 1);
        parfor k = 1:nb
            fr = blocks{k};  wk = warm0;  f_ = false(numel(fr), 1);  m_ = cell(numel(fr), 1);  i_ = cell(numel(fr), 1);
            for j = 1:numel(fr)
                try
                    [sx, sy] = cisg_cache_frame(fr(j), rb, st, ff, cd_, bn, dtt, fs_, wk);
                    if ~isempty(wk), wk.prev = struct('n', fr(j), 'Sx', sx, 'Sy', sy); end
                    if ~isempty(dq), send(dq, struct('n', fr(j), 'Sx', sx(1:8:end, 1:8:end))); end
                catch ME                                     % record the failure (identifier + message), judged below
                    f_(j) = true;  m_{j} = ME.message;  i_{j} = ME.identifier;
                    if ~isempty(wk), wk.prev = []; end
                end
            end
            bf{k} = f_;  bm{k} = m_;  bi{k} = i_;
        end
        bad_f = vertcat(bf{:});  bad_m = vertcat(bm{:});  bad_id = vertcat(bi{:});
        if live, close(fp); end
        % only 'cisg:badFrame' (empty / truncated file) is skippable; anything else is a bug: fail loudly
        other = bad_f & ~strcmp(bad_id, 'cisg:badFrame');
        if any(other)
            [um, iu] = unique(bad_m(other));  fo = todo(other);  io = find(other);
            fprintf(2, '%d frame(s) failed with an unexpected error. First distinct message(s):\n', nnz(other));
            for k = 1:min(3, numel(um)), fprintf(2, '  [%s] frame %d: %s\n', bad_id{io(iu(k))}, fo(iu(k)), um{k}); end
            error('cisg_slope_video:workerError', 'Aborting: unexpected error in %d of %d frames (%s)', nnz(other), numel(todo), um{1});
        end
        bad_list = todo(bad_f);
        bad_msgs = arrayfun(@(k) bad_reason(ff(todo(k)), bad_m{k}), find(bad_f)', 'UniformOutput', false);
        fprintf('parfor: %d frames in %.1f s (%.1f s/frame wall)\n', numel(todo), toc(t0), toc(t0) / numel(todo));
        report_bad(bad_list, bad_msgs, cache_dir);
    elseif ~have_pct
        fprintf('Parallel Computing Toolbox not available: computing serially in the render loop\n');
    end
end
%% ---------------- Single frame: just the figure, no video ----------------
if isscalar(frames)
    n = frames;
    warm = [];  if cfg.warm_start, warm = struct('prev', [], 'lead', cfg.warm_lead); end
    [Sx, Sy] = cisg_cache_frame(n, ref_b, setup_b, R.frame_file, cache_dir, cfg.bin, R.dt, R.img_size, warm);
    obs = cisg_load_coreview(R.frame_file(n), R.img_size);
    if cfg.bin > 1, obs = uint8(imresize(obs, 1 / cfg.bin, 'box')); end
    fig = figure('Name', sprintf('CISG CoreView %d frame %d', cfg.coreview_num, n), ...
                 'Position', [50 50 cfg.fig_size]);
    opts = struct('bin', cfg.bin, 'lim_s', cfg.lim_s, 'lim_eta', cfg.lim_eta, ...
                  'head', sprintf('CoreView %d, frame %d, %s', cfg.coreview_num, n, tstr(n)), ...
                  'ts', ts, 't_now', R.t_wind(n), 'ts_xlim', cfg.series_xlim);
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
opts = struct('bin', cfg.bin, 'lim_s', cfg.lim_s, 'lim_eta', cfg.lim_eta, 'ts', ts, 't_now', [], 'ts_xlim', cfg.series_xlim, ...
              'eta_hist', struct('t', [], 'v', []));       % CISG eta at (0, 0) of the frames so far, over the wave gauge
have_good = false;
warm = [];  if cfg.warm_start, warm = struct('prev', [], 'lead', cfg.warm_lead); end
tic;
for i = 1:numel(frames)
    n = frames(i);
    t1 = tic;
    blank = false;  missing = '';
    try
        [Sx_n, Sy_n] = cisg_cache_frame(n, ref_b, setup_b, R.frame_file, cache_dir, cfg.bin, R.dt, R.img_size, warm);
        if ~isempty(warm), warm.prev = struct('n', n, 'Sx', Sx_n, 'Sy', Sy_n); end
        obs_n = cisg_load_coreview(R.frame_file(n), R.img_size);
        if cfg.bin > 1, obs_n = uint8(imresize(obs_n, 1 / cfg.bin, 'box')); end
        Sx = Sx_n;  Sy = Sy_n;  obs = obs_n;  have_good = true;
    catch ME
        if ~strcmp(ME.identifier, 'cisg:badFrame'), rethrow(ME); end   % only a bad file is skippable
        reason = bad_reason(R.frame_file(n), ME.message);
        if ~ismember(n, bad_list), bad_list(end+1) = n;  bad_msgs{end+1} = reason; end %#ok<SAGROW>
        warning('cisg_slope_video:badFrame', 'Frame %d unreadable (%s)', n, reason);
        switch cfg.bad_frame
            case 'skip',  continue;
            case 'hold'
                if ~have_good, continue; end
                missing = ' (missing)';                      % repeat the previous good frame
            otherwise                                        % 'blank': same layout, empty panels
                blank = true;
                nr = setup_b.roi.row_max - setup_b.roi.row_min + 1;  nc = setup_b.roi.col_max - setup_b.roi.col_min + 1;
                Sx = nan(nr, nc);  Sy = Sx;  obs = uint8(128 * ones(size(ref_b)));
        end
    end
    if blank
        opts.head = sprintf('CoreView %d, frame %d (%s) - MISSING (%s)', cfg.coreview_num, n, tstr(n), reason);
    else
        opts.head = sprintf('CoreView %d, frame %d%s, %s', cfg.coreview_num, n, missing, tstr(n));
    end
    opts.blank = blank;  opts.t_now = R.t_wind(n);
    info = cisg_plot_frame(fig, obs, double(Sx), double(Sy), setup_b, R.cal, opts);
    opts.eta_hist.t(end+1) = opts.t_now;  opts.eta_hist.v(end+1) = info.eta_centre_mm;   % NaN for a missing frame
    if blank                                                % big label on the raw panel
        ax1 = info.ax_raw;
        text(ax1, mean(xlim(ax1)), mean(ylim(ax1)), sprintf('frame %d missing', n), 'Color', 'w', 'FontSize', 28, ...
             'FontWeight', 'bold', 'HorizontalAlignment', 'center');
    else                                                     % hold the colour limits from the first good frame
        if isempty(opts.lim_s),   opts.lim_s   = info.lim_s;   fprintf('slope colour limit +-%.3f (held)\n', opts.lim_s); end
        if isempty(opts.lim_eta), opts.lim_eta = info.lim_eta; fprintf('eta colour limit +-%.2f mm (held)\n', opts.lim_eta); end
    end
    fig.Position(3:4) = cfg.fig_size;
    if live, drawnow; end
    writeVideo(vw, print(fig, '-RGBImage', '-r96'));            % getframe fails on an invisible figure
    el = toc;
    fprintf('  frame %d (%d/%d): %.1f s, rms Sx %.3f Sy %.3f eta %.2f mm, masked %.1f%% | elapsed %.0f s, ETA %.0f s\n', ...
            n, i, numel(frames), toc(t1), info.rms_sx, info.rms_sy, info.rms_eta, info.masked, el, el / i * (numel(frames) - i));
end
close(vw);  close(fig);
report_bad(bad_list, bad_msgs, cache_dir);
d = dir(cfg.out_path);
if isempty(d)
    fprintf('No video written (no frame was rendered): %s\n', cfg.out_path);
else
    fprintf('Wrote %s (%.1f MB), %.1f s total\n', cfg.out_path, d.bytes / 1e6, toc);
end

function show_progress(fp, d, total, t0)
    % DataQueue callback (client): thumbnail of the slopes a worker just finished
    fp.UserData = fp.UserData + 1;
    k = fp.UserData;
    ax = fp.CurrentAxes;
    imagesc(ax, d.Sx, 'AlphaData', ~isnan(d.Sx));  axis(ax, 'image');  colormap(ax, 'parula');
    set(ax, 'Color', [0.8 0.8 0.8]);                             % NaN (masked) shows grey
    el = toc(t0);
    title(ax, sprintf('computing slopes: frame %d (%d/%d), ETA %.0f s', d.n, k, total, el / k * (total - k)));
    drawnow limitrate;
end

function report_bad(bad_list, bad_msgs, cache_dir)
    % print the unreadable frames and save "frame  reason" lines to <cache_dir>/bad_frames.txt (merged with an earlier list)
    if isempty(bad_list), return; end
    fprintf('%d unreadable frame(s): %s\n', numel(bad_list), mat2str(sort(bad_list)));
    if isempty(cache_dir), return; end
    f = fullfile(cache_dir, 'bad_frames.txt');
    lines = containers.Map('KeyType', 'double', 'ValueType', 'char');
    if isfile(f)
        for L = strsplit(strtrim(fileread(f)), newline)
            k = sscanf(L{1}, '%d', 1);
            if ~isempty(k), lines(k) = strtrim(L{1}(find(L{1} == ' ', 1):end)); end
        end
    end
    for k = 1:numel(bad_list), lines(bad_list(k)) = bad_msgs{k}; end
    fid = fopen(f, 'w');
    for k = sort(cell2mat(keys(lines)))
        fprintf(fid, '%d  %s\n', k, lines(k));
    end
    fclose(fid);
    fprintf('saved %s\n', f);
end

function s = bad_reason(file, msg)
    % short reason for a bad frame: 'empty file' (header only), 'truncated file (N bytes)' or the error text
    d = dir(file);
    if isempty(d), s = 'file not found';
    elseif d.bytes <= 20, s = 'empty file';
    elseif contains(msg, 'Bad CoreView frame'), s = sprintf('bad or truncated file (%d bytes)', d.bytes);
    else, s = strtok(msg, newline);
    end
end
