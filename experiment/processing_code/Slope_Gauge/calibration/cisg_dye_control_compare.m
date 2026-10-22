%% cisg_dye_control_compare.m
% Dye artefact or real convergence lines?  Streamwise stripes in Sx / Sy for
%   run 95 (full rhodamine 6G), run 96 (less R6G, same calibration/height 48 cm),
%   run 101 (no R6G, height 50.1 cm, plate frame CoreView 106) at matched times (same rpm ramp).
% Per frame (bin 2, 2-channel colour inversion, 'auto' ROI per run; joint dye solve also for 95 and 96):
% rms slopes, stripe strength = var_y(streamwise row-mean) / total variance, mean-slope and variance
% profiles vs y.  Slopes are computed here but nothing is cached on D:.  Writes calibration/cisg_dye_control_compare.png.
% Uses a parallel pool of cfg.n_workers (the user's own MATLAB must be idle).
clear; clc;
here = fileparts(mfilename('fullpath'));
sg = fileparts(here);
addpath(here, fullfile(sg, 'util'), sg, fullfile(fileparts(sg), 'IR Camera', 'calibration'));

cfg.runs      = [95 96 101];
cfg.n_workers = 12;
cfg.dense     = {1460:1559, 2400:2449, 3000:3049};      % every frame (2 s, 1 s, 1 s) at matched times
cfg.stride    = 10;                                     % sparse sampling over the whole run
cfg.results   = fullfile(getenv('TEMP'), 'cisg_dye_control_results.mat');

%% ---- per-run setup ----
nr = numel(cfg.runs);  RR = cell(1, nr);
for i = 1:nr
    run_n = cfg.runs(i);
    S = cisg_scan_run(run_n);
    R = cisg_setup_run(struct('coreview_num', run_n, 'roi_mode', 'auto', 'auto_margin', 100));
    [ref_b, st] = cisg_bin_setup(R.ref_image, R.setup, 2);
    r = st.roi;  rows = r.row_min:r.row_max;
    uv = [repmat(mean([r.col_min r.col_max]), numel(rows), 1), rows(:)];
    xy = surf_px2world(R.cal, [(uv(:,1) - 0.5) * 2 + 0.5, (uv(:,2) - 0.5) * 2 + 0.5]);
    RR{i} = struct('run', run_n, 'scan', S, 'ref_b', ref_b, 'st', st, 'R', R, 'ycm', xy(:,2) / 10);
    fprintf('run %d: L_eff %.1f cm, H %.1f cm, r %.3f, h_est %.1f cm, ROI rows %d-%d cols %d-%d (full res), cm/px %.5f, y range %.1f..%.1f cm\n', ...
            run_n, R.setup.L_eff_cm, R.depth.H_mm/10, R.depth.r, R.depth.h_est_mm/10, R.setup.roi.row_min, R.setup.roi.row_max, ...
            R.setup.roi.col_min, R.setup.roi.col_max, R.setup.cm_per_pixel_x, min(xy(:,2))/10, max(xy(:,2))/10);
end
ylo = max(cellfun(@(q) min(q.ycm), RR));  yhi = min(cellfun(@(q) max(q.ycm), RR));
yq = (ceil(ylo * 4) / 4 : 0.25 : floor(yhi * 4) / 4)';
nmax = min(cellfun(@(q) numel(q.scan.frames), RR));
frames = unique([1:cfg.stride:nmax, cfg.dense{:}]);  frames = frames(frames <= nmax);
fprintf('\ncommon y grid %.2f..%.2f cm (%d points); %d frames per run (%d sparse + dense windows)\n', yq(1), yq(end), numel(yq), numel(frames), numel(1:cfg.stride:nmax));

%% ---- compute ----
pool = gcp('nocreate');
if isempty(pool) || pool.NumWorkers ~= cfg.n_workers, delete(pool); parpool('local', cfg.n_workers); end
pctRunOnAll(sprintf('addpath(''%s''); cisg_fix_path(''%s'');', fullfile(sg, 'util'), sg));
M = cell(1, nr);
for i = 1:nr
    q = RR{i};  t0 = tic;
    ff = q.R.frame_file;  rb = q.ref_b;  st = q.st;  fs_ = q.R.img_size;  yc = q.ycm;  dj = q.run ~= 101;
    out = cell(1, numel(frames));
    parfor k = 1:numel(frames)
        out{k} = cisg_frame_metrics(frames(k), rb, st, ff, 2, fs_, yc, yq, dj);
    end
    M{i} = out;
    ok = cellfun(@(m) m.ok, out);
    fprintf('run %d: %d frames in %.0f s (%d unreadable: %s)\n', q.run, numel(frames), toc(t0), nnz(~ok), mat2str(frames(~ok)));
end
save(cfg.results, 'M', 'RR', 'frames', 'yq', 'cfg', '-v7.3');

%% ---- tables ----
seg = [1 1000; 1001 2000; 2001 3000; 1460 1559];  segn = {'frames 1-1000', '1001-2000', '2001-3000', 'dense 1460-1559'};
fprintf('\nmetrics (time-mean over the frames in each segment; 2-channel inversion; str = stripe strength var_y(row-mean)/var_total)\n');
fprintf('%-18s %-6s | rms Sx  rms Sy | strSx   strSy  | joint: strSx  strSy | mean c   corr(c,rmSx) corr(c,rmSy)\n', 'segment', 'run');
for s = 1:size(seg, 1)
    for i = 1:nr
        out = M{i};  ok = cellfun(@(m) m.ok, out) & frames >= seg(s, 1) & frames <= seg(s, 2);
        g = @(f) mean(cellfun(@(m) m.(f), out(ok)), 'omitnan');
        if RR{i}.run ~= 101
            fprintf('%-18s %-6d | %.3f  %.3f  | %.4f %.4f | %.4f %.4f | %.4f  %+.2f  %+.2f   (c: p99 %.3f)\n', segn{s}, RR{i}.run, g('rmsSx'), g('rmsSy'), g('strSx'), g('strSy'), g('strSxj'), g('strSyj'), ...
                    g('c_mean'), g('corr_c_Sx'), g('corr_c_Sy'), g('c_p99'));
        else
            fprintf('%-18s %-6d | %.3f  %.3f  | %.4f %.4f |\n', segn{s}, RR{i}.run, g('rmsSx'), g('rmsSy'), g('strSx'), g('strSy'));
        end
    end
end

%% ---- figure ----
col = {[0.85 0.2 0.2], [0.9 0.55 0.1], [0.1 0.4 0.8]};  nm = {'95 (full R6G)', '96 (less R6G)', '101 (no R6G)'};
fig = figure('Visible', 'off', 'Position', [50 50 1700 950]);
ts = frames / 50;
sm = @(v) movmean(v, 5, 'omitnan');
panels = {'rmsSx', 'rms S_x'; 'rmsSy', 'rms S_y'; 'strSx', 'stripe strength S_x'; 'strSy', 'stripe strength S_y'};
for p = 1:4
    subplot(2, 3, p); hold on;
    for i = 1:nr
        ok = cellfun(@(m) m.ok, M{i});  v = nan(1, numel(frames));  v(ok) = cellfun(@(m) m.(panels{p, 1}), M{i}(ok));
        sel = frames <= 3000 & mod(frames, cfg.stride) == 1;                % sparse samples only (a clean time series)
        plot(ts(sel), sm(v(sel)), '-', 'Color', col{i}, 'LineWidth', 1.3);
        if p > 2 && RR{i}.run ~= 101
            vj = nan(1, numel(frames));  vj(ok) = cellfun(@(m) m.([panels{p, 1} 'j']), M{i}(ok));
            plot(ts(sel), sm(vj(sel)), '--', 'Color', col{i});
        end
    end
    xlabel('t (s)'); title(panels{p, 2}); grid on;
    if p == 3, legend([nm, {'(dashed: joint dye solve)'}], 'Location', 'northwest'); end
end
dense = frames >= 1460 & frames <= 1559;
subplot(2, 3, 5); hold on;
for i = 1:nr
    out = M{i}(dense & cellfun(@(m) m.ok, M{i}));
    plot(yq, mean(cell2mat(cellfun(@(m) m.rmSx(:)', out(:), 'UniformOutput', false)), 1, 'omitnan'), 'Color', col{i}, 'LineWidth', 1.3);
    plot(yq, mean(cell2mat(cellfun(@(m) m.rmSy(:)', out(:), 'UniformOutput', false)), 1, 'omitnan'), ':', 'Color', col{i}, 'LineWidth', 1.3);
end
xlabel('y (cm)'); title('streamwise-mean S_x (solid) and S_y (dotted) vs y, frames 1460-1559'); grid on;
subplot(2, 3, 6); hold on;
for i = 1:nr
    out = M{i}(dense & cellfun(@(m) m.ok, M{i}));
    plot(yq, mean(cell2mat(cellfun(@(m) m.rvSx(:)', out(:), 'UniformOutput', false)), 1, 'omitnan'), 'Color', col{i}, 'LineWidth', 1.3);
    plot(yq, mean(cell2mat(cellfun(@(m) m.rvSy(:)', out(:), 'UniformOutput', false)), 1, 'omitnan'), ':', 'Color', col{i}, 'LineWidth', 1.3);
end
xlabel('y (cm)'); title('variance along x of S_x (solid), S_y (dotted) vs y'); grid on;
exportgraphics(fig, fullfile(here, 'cisg_dye_control_compare.png'), 'Resolution', 90); close(fig);
fprintf('\nSaved calibration/cisg_dye_control_compare.png\n');
