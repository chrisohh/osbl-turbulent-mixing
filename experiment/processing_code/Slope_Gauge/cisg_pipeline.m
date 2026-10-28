%% CISG pipeline: every step from raw files to slope cache / video, for one run
% Run section by section (Ctrl+Enter) in order; each section can be re-run on
% its own after section 0.  Per-run files (plate frame, IR snapshot, IR flat,
% plate parity, dye model) come from util/cisg_run_defaults.m: add a new run
% there, not here.  Everything printed is also appended to
% calibration/cisg_pipeline_CoreView_<N>.log (MATLAB diary), so the numbers of
% every pass (rms, laser shift, L_eff, ...) are kept.
%
% Steps that need clicks: 2 (IR, 4 corners, first time only: seeds are saved),
% 3 (laser-sheet line).  Steps that overwrite files back the old file up first
% (<name>_<yyyymmdd_HHMM>.mat).
%
%   0  settings, paths, log
%   1  scan raw files (bad / headerless frames)
%   2  IR surface calibration       (skipped if the IR cal exists and redo_ir = false)
%   3  CISG plate calibration       (laser line -> x = 0; skipped if done and redo_plate = false)
%   4  checks: IR <-> CISG colocation, apparent depth (L_eff, h_est)
%   5  preview one frame (nothing saved); 5b why pixels are masked; 5c warm start / outlier test vs cold start
%   6  slope cache + video (cisg_slope_video with the settings below)
%   7  optional, skipped by default: space-time outlier pass over the cache

%% 0. Settings
clear P;
P.run_num     = 101;             % CoreView run
P.redo_ir     = false;           % true = redo the IR surface calibration (needs rd.ir_flat_path)
P.redo_plate  = false;           % true = redo the CISG plate calibration even if the laser line is set
P.laser_x_px  = [];              % [] = click the laser line; a number = reuse that image column (no click)

P.video = struct('frames', [], 'roi_mode', 'full', 'bin', 2, 'dye_model', '', ...
                 'use_parfor', true, 'n_workers', 'auto', 'show_live', true, 'bad_frame', 'blank', ...
                 'cache_slopes', true, 'out_root', 'D:\HLAB_2026\SlopeGauge', ...
                 'warm_start', false, 'warm_lead', 3, ...  % previous-frame warm start (off: 35.6% -> 33.9% masked on 101/2000 for ~3x compute)
                 'ref_frames', 2:4:58, ...      % reference = median of these flat frames (frame 1 is ~0.8 counts brighter: section 5d)
                 'bias_frames', 60:10:200, ...  % flat frames whose mean slope is subtracted from every frame ([] = off)
                 'drift_corr', 'none', ...      % per-frame brightness correction: 'none' | 'offset' | 'gain' (section 5e)
                 'colour_channels', 'hybrid');  % 'RG' | 'RGB' | 'hybrid' (R, G; R, G, B only where R, G fails or folds: 5f, 5h, 5i);
                                                % blue only for dye-free runs: 'RG' for rhodamine (95, 96)
P.outlier = struct('q', 8, 'w', 5, 'thr', 3, 'eps0', 0.01);   % space-time median test (section 7, cisg_outlier_mask)

% Run Section (Ctrl+Enter) executes a temp copy, so mfilename points to a temp folder: take the first candidate
% folder that holds util/cisg_run_defaults.m
c_ = {fileparts(mfilename('fullpath')), fileparts(which('cisg_pipeline')), pwd};
try c_{end+1} = fileparts(matlab.desktop.editor.getActiveFilename); catch, end
c_ = c_(cellfun(@(d) ~isempty(d) && isfile(fullfile(d, 'util', 'cisg_run_defaults.m')), c_));
if isempty(c_), error('cisg_pipeline:here', 'Open cisg_pipeline.m in the editor or cd to Slope_Gauge first'); end
P.here = c_{1};  clear c_
addpath(fullfile(P.here, 'util'), fullfile(P.here, 'calibration'), fullfile(P.here, '..', 'IR Camera', 'calibration'));
cisg_fix_path(P.here);
P.rd = cisg_run_defaults(P.run_num);
P.log = fullfile(P.here, 'calibration', sprintf('cisg_pipeline_CoreView_%d.log', P.run_num));
diary(P.log);  diary on;
fprintf('\n========== cisg_pipeline CoreView %d, %s ==========\n', P.run_num, datestr(now, 'yyyy-mm-dd HH:MM'));
fprintf('plate CoreView_%d, IR Snap-%06d, ablated parity %d, dye model %s\n', ...
        P.rd.plate, P.rd.ir_snap, P.rd.ablated_parity, P.rd.dye_model);
P.backup = @(f) copyfile(f, regexprep(f, '\.mat$', ['_' datestr(now, 'yyyymmdd_HHMM') '.mat']));

%% 1. Scan raw files
fprintf('\n--- 1. raw file scan ---\n');
P.scan = cisg_scan_run(P.run_num);

%% 2. IR surface calibration (4 corner clicks the first time: TL, TR, BR, BL inner corners)
fprintf('\n--- 2. IR surface calibration (%s) ---\n', P.rd.ir_snap_path);
if isfile(P.rd.ir_cal_path) && ~P.redo_ir
    fprintf('kept %s (redo_ir = false)\n', P.rd.ir_cal_path);
elseif isempty(P.rd.ir_flat_path)
    error('cisg_pipeline:flat', 'No board-free IR recording for Snap-%06d: set ir_flat_path in cisg_run_defaults', P.rd.ir_snap);
else
    if isfile(P.rd.ir_cal_path), P.backup(P.rd.ir_cal_path); end
    vr = FlirMovieReader(P.rd.ir_flat_path);  vr.unit = 'temperatureFactory';
    acc = 0;  for i = 1:50, acc = acc + double(step(vr)); end
    flat = acc / 50;  clear vr acc i
    fprintf('flat: mean of 50 frames of %s\n', P.rd.ir_flat_path);
    cal_ir = calibrate_surface_single(P.rd.ir_snap_path, struct('flat', flat, 'seed_clicks', true));
    fprintf('IR cal: rms %.3f px -> %s\n', cal_ir.rms_px, P.rd.ir_cal_path);
end

%% 3. CISG plate calibration (click the laser-sheet line: x = 0 there)
fprintf('\n--- 3. CISG plate calibration (%s) ---\n', P.rd.plate_frame);
done_ = false;
if isfile(P.rd.cal_path)
    S_ = load(P.rd.cal_path, 'cal');
    done_ = isfield(S_.cal, 'laser_x_px') && ~isnan(S_.cal.laser_x_px);
end
if done_ && ~P.redo_plate
    fprintf('kept %s (laser line at column %.0f; redo_plate = false)\n', P.rd.cal_path, S_.cal.laser_x_px);
else
    if isfile(P.rd.cal_path), P.backup(P.rd.cal_path); end
    cal = cisg_plate_calibrate(struct('frame_path', P.rd.plate_frame, 'ir_cal_path', P.rd.ir_cal_path, ...
                                      'out_path', P.rd.cal_path, 'ablated_parity', P.rd.ablated_parity, ...
                                      'use_ir', isfile(P.rd.ir_cal_path), 'laser_x_px', P.laser_x_px));
    fprintf('plate cal: rms %.3f px, theta %.2f deg, %.5f mm/px at centre\n', cal.rms_px, cal.theta_deg, cal.mm_per_px_center);
    fprintf('laser line: clicked column %.0f -> snapped %.0f, x origin moved %.0f mm\n', ...
            cal.laser_x_px, cal.laser_snapped_px, cal.laser_x_shift_mm);
end
clear S_ done_

%% 4. Checks: colocation and apparent depth
fprintf('\n--- 4. checks ---\n');
run_num = P.run_num;
run(fullfile(P.here, 'calibration', 'cisg_check_colocation.m'));     % keeps run_num and P; writes the colocation png
S_ = load(P.rd.cal_path, 'cal');
ref_ = cisg_load_coreview(fullfile('D:\HLAB_2026\Core', sprintf('CoreView_%d', P.run_num), 'Flare 12M125 CCL C1112A00', ...
                          sprintf('CoreView_%d_Flare 12M125 CCL C1112A00_0001.raw', P.run_num)), [3072 4100]);
P.depth = cisg_apparent_depth(ref_, S_.cal, struct('fig_path', ...
          fullfile(P.here, 'calibration', sprintf('cisg_apparent_depth_CoreView_%d.png', P.run_num))));
fprintf('L_eff %.1f mm, h_est %.1f mm, H %.0f mm\n', P.depth.L_eff_mm, P.depth.h_est_mm, P.depth.H_mm);
clear S_ ref_ run_num

%% 5. Preview one frame (nothing saved)
P.preview_frame = 10;
fprintf('\n--- 5. preview frame %d ---\n', P.preview_frame);
R_ = cisg_setup_run(struct('coreview_num', P.run_num, 'roi_mode', P.video.roi_mode, 'dye_model', P.video.dye_model, ...
                           'bin', P.video.bin, 'ref_frames', P.video.ref_frames, 'setup', struct('drift_corr', P.video.drift_corr, 'colour_channels', P.video.colour_channels)));
[ref_b_, setup_b_] = cisg_bin_setup(R_.ref_image, R_.setup, P.video.bin);
setup_b_ = cisg_flat_bias(R_, ref_b_, setup_b_, P.video.bin, P.video.bias_frames);   % same flat-water bias as section 6
[Sx_, Sy_] = cisg_cache_frame(P.preview_frame, ref_b_, setup_b_, R_.frame_file, '', P.video.bin, R_.dt, R_.img_size);
obs_ = cisg_load_coreview(R_.frame_file(P.preview_frame), R_.img_size);
if P.video.bin > 1, obs_ = uint8(imresize(obs_, 1 / P.video.bin, 'box')); end
ts_ = cisg_wind_wave_series(P.run_num);              % wind + wave gauge panels (as in the video)
info_ = cisg_plot_frame(figure('Name', sprintf('CISG preview CoreView %d', P.run_num), 'Position', [50 50 1600 1000]), ...
                        obs_, double(Sx_), double(Sy_), setup_b_, R_.cal, struct('bin', P.video.bin, ...
                        'ts', ts_, 't_now', R_.t_wind(P.preview_frame), ...
                        'head', sprintf('CoreView %d, frame %d, t = %.2f s after wind start (camera %.2f s)', P.run_num, ...
                                P.preview_frame, R_.t_wind(P.preview_frame), (P.preview_frame - 1) * R_.dt)));
fprintf('frame %d: rms Sx %.3f Sy %.3f eta %.2f mm, masked %.1f%%\n', ...
        P.preview_frame, info_.rms_sx, info_.rms_sy, info_.rms_eta, info_.masked);
clear Sx_ Sy_ info_

%% 5b. Why pixels are masked (preview frame; run after section 5)
% Each masked pixel gets the first reason in this order; the bars count every reason a pixel has.
[~, ~, cinfo_, D_] = cisg_colour_displacement(ref_b_, obs_, setup_b_);
names_ = {'dark in reference', 'shift out of image', 'shift clamped', 'ill-conditioned J', 'not converged', ...
          'colour residual', 'precision sigma_S > max'};
vals_  = [64 1 2 4 8 16 128];
first_ = zeros(size(D_.bits));
for k_ = numel(vals_):-1:1, first_(bitand(D_.bits, vals_(k_)) > 0) = k_; end
frac_ = arrayfun(@(v) 100 * mean(bitand(D_.bits(:), v) > 0), vals_);
fprintf('masked %.1f%% of the ROI; by reason (overlapping):\n', 100 * mean(D_.bits(:) > 0));
for k_ = 1:numel(vals_), fprintf('  %-26s %5.1f%%\n', names_{k_}, frac_(k_)); end
figure('Name', sprintf('CISG mask reasons CoreView %d frame %d', P.run_num, P.preview_frame), 'Position', [80 80 1500 520]);
subplot(1, 3, 1);
imagesc(first_, [-0.5 numel(vals_) + 0.5]);  axis image;
colormap(gca, [1 1 1; lines(numel(vals_))]);
colorbar('Ticks', 0:numel(vals_), 'TickLabels', [{'kept'} names_]);
title('first mask reason');
subplot(1, 3, 2);
if cinfo_.sigma_mask
    imagesc(log10(max(cinfo_.sig_sx, cinfo_.sig_sy)), log10([0.001 0.2]));  axis image;  colorbar;
    title(sprintf('log_{10} \\sigma_S (mask above %.3f)', setup_b_.max_slope_sigma));
else
    title('precision mask OFF (setup lacks noise / scale fields)');
end
subplot(1, 3, 3);
r_ = setup_b_.roi;
imagesc(mean(double(ref_b_(r_.row_min:r_.row_max, r_.col_min:r_.col_max, :)), 3) / 255, [0 0.3]);  axis image;  colorbar;
title('reference brightness (mean RGB / 255)');
clear cinfo_ D_ names_ vals_ first_ frac_ k_ r_ R_ ref_b_ setup_b_ obs_

%% 5c. Test: warm start and outlier test vs cold start (frames 1995-2005, nothing saved, ~5-10 min)
P.test_frames = 1995:2005;                         % consecutive; the middle one is plotted
if ~isfield(P.video, 'warm_lead'), P.video.warm_lead = 3; end   % P from before section 0 had warm start
fprintf('\n--- 5c. warm-start test, frames %d-%d ---\n', P.test_frames(1), P.test_frames(end));
R_ = cisg_setup_run(struct('coreview_num', P.run_num, 'roi_mode', P.video.roi_mode, 'dye_model', P.video.dye_model, ...
                           'bin', P.video.bin, 'ref_frames', P.video.ref_frames, 'setup', struct('drift_corr', P.video.drift_corr, 'colour_channels', P.video.colour_channels)));
[ref_b_, setup_b_] = cisg_bin_setup(R_.ref_image, R_.setup, P.video.bin);
setup_b_ = cisg_flat_bias(R_, ref_b_, setup_b_, P.video.bin, P.video.bias_frames);
cf_ = @(n, w) cisg_cache_frame(n, ref_b_, setup_b_, R_.frame_file, '', P.video.bin, R_.dt, R_.img_size, w);
nf_ = numel(P.test_frames);  cold_ = cell(nf_, 2);  warm_ = cell(nf_, 2);  wk_ = struct('prev', [], 'lead', 0);
fprintf('frame   masked cold   masked warm   from previous\n');
for i_ = 1:nf_
    n_ = P.test_frames(i_);
    [cold_{i_, 1}, cold_{i_, 2}] = cf_(n_, []);
    [warm_{i_, 1}, warm_{i_, 2}, wi_] = cf_(n_, wk_);
    wk_.prev = struct('n', n_, 'Sx', warm_{i_, 1}, 'Sy', warm_{i_, 2});
    fprintf('%5d   %8.1f%%   %10.1f%%   %10.1f%%\n', n_, 100 * mean(isnan(cold_{i_, 1}), 'all'), ...
            100 * mean(isnan(warm_{i_, 1}), 'all'), 100 * wi_.warm_fraction);
end
im_ = ceil(nf_ / 2);  nm_ = P.test_frames(im_);
[ldx_, ~] = cf_(nm_, struct('prev', [], 'lead', P.video.warm_lead));   % stride > 1 case: lead-in only, no solved n-1
fprintf('frame %d with only a %d-frame lead-in (stride-10 case): masked %.1f%%\n', nm_, P.video.warm_lead, 100 * mean(isnan(ldx_), 'all'));
[out_, oi_] = cisg_outlier_mask(warm_{im_, 1}, warm_{im_, 2}, struct('Sx', warm_{im_-1, 1}, 'Sy', warm_{im_-1, 2}), ...
                                struct('Sx', warm_{im_+1, 1}, 'Sy', warm_{im_+1, 2}), P.outlier);
fprintf('outlier test on frame %d (warm): %.2f%% of the ROI flagged\n', nm_, 100 * oi_.pixel_fraction);
both_ = isfinite(cold_{im_, 1}) & isfinite(warm_{im_, 1});
fprintf('where cold and warm both kept a pixel: rms difference Sx %.4f (rms Sx %.3f)\n', ...
        rms(cold_{im_, 1}(both_) - warm_{im_, 1}(both_)), rms(warm_{im_, 1}(both_)));
sw_ = warm_{im_, 1};  sw_(out_) = NaN;
lim_ = 3 * rms(warm_{im_, 1}, 'all', 'omitnan');
figure('Name', sprintf('Warm start test CoreView %d frame %d', P.run_num, nm_), 'Position', [60 60 1600 480]);
pan_ = {cold_{im_, 1}, warm_{im_, 1}, sw_, double(out_) + 2 * double(isnan(cold_{im_, 1}) & isfinite(warm_{im_, 1}))};
ttl_ = {'S_x cold start', 'S_x warm start', 'S_x warm + outlier mask', 'recovered by warm (2) / outliers (1)'};
for k_ = 1:4
    subplot(1, 4, k_);
    if k_ < 4
        imagesc(pan_{k_}, 'AlphaData', isfinite(pan_{k_}), [-lim_ lim_]);  set(gca, 'Color', [0.8 0.8 0.8]);
    else
        imagesc(pan_{k_}, [0 3]);  colormap(gca, [1 1 1; 1 0 0; 0 0.6 0; 0 0 0]);
    end
    axis image;  colorbar;  title(ttl_{k_});
end
clear R_ ref_b_ setup_b_ cf_ nf_ cold_ warm_ wk_ wi_ i_ n_ im_ nm_ ldx_ out_ oi_ both_ sw_ lim_ pan_ ttl_ k_

%% 5d. Diagnostic: flat-water bias against the reference (lighting drift? bad reference?), ~2 min, nothing saved
% Early frames are flat, so their slopes vs the reference should be ~0.  If the frame-mean Sx, Sy follow the
% brightness change of R / G, the bias is lighting drift.  The crops around a point (P.diag_xy, tank cm as
% plotted) show which frame differs: a feature in |f2 - f1| and |f10 - f1| but not in |f10 - f2| is in frame 1 only.
% P.diag_ref = 1 repeats the original test (frame 1 alone); the default checks the median reference (no bias
% subtraction here, so the table shows what the reference alone leaves).
P.diag_frames = [2 10 50 100 200];
P.diag_xy     = [-10 7];                           % the persistent patch (x, y in cm, as on the maps)
P.diag_ref    = P.video.ref_frames;
fprintf('\n--- 5d. flat-water bias, reference = median of frames %s ---\n', mat2str(P.diag_ref));
R_ = cisg_setup_run(struct('coreview_num', P.run_num, 'roi_mode', P.video.roi_mode, 'dye_model', P.video.dye_model, ...
                           'bin', P.video.bin, 'ref_frames', P.diag_ref, 'setup', struct('drift_corr', P.video.drift_corr, 'colour_channels', P.video.colour_channels)));
[ref_b_, setup_b_] = cisg_bin_setup(R_.ref_image, R_.setup, P.video.bin);
rd_ = setup_b_.dark_mask;  if isempty(rd_), rd_ = false(size(ref_b_, 1), size(ref_b_, 2)); end
ok_ = ~rd_;                                         % on-target pixels for the brightness means
ld_ = @(n) uint8(imresize(cisg_load_coreview(R_.frame_file(n), R_.img_size), 1 / P.video.bin, 'box'));
fprintf('frame  t_wind   mean Sx    mean Sy    rms Sx   rms Sy  |  dR      dG      dB   (counts vs reference, on target)\n');
for n_ = P.diag_frames
    [sx_, sy_] = cisg_cache_frame(n_, ref_b_, setup_b_, R_.frame_file, '', P.video.bin, R_.dt, R_.img_size);
    ob_ = double(ld_(n_));  dc_ = zeros(1, 3);
    for c_ = 1:3, d_ = ob_(:, :, c_) - double(ref_b_(:, :, c_));  dc_(c_) = mean(d_(ok_)); end
    fprintf('%5d  %6.2f  %+8.4f  %+8.4f  %7.4f  %7.4f  |  %+6.2f  %+6.2f  %+6.2f\n', n_, R_.t_wind(n_), ...
            mean(sx_, 'all', 'omitnan'), mean(sy_, 'all', 'omitnan'), rms(sx_, 'all', 'omitnan'), rms(sy_, 'all', 'omitnan'), dc_);
end
% crops around the patch: tank cm (as plotted, y relative to the channel centre) -> binned px
y0_ = 0;  if isfield(R_.setup, 'y0_cm'), y0_ = R_.setup.y0_cm; end
uv_ = surf_world2px(R_.cal, 10 * [P.diag_xy(1), P.diag_xy(2) + y0_]) / P.video.bin;
hw_ = 120;  rr_ = max(1, round(uv_(2)) - hw_):min(size(ref_b_, 1), round(uv_(2)) + hw_);
cc_ = max(1, round(uv_(1)) - hw_):min(size(ref_b_, 2), round(uv_(1)) + hw_);
f1_ = double(ld_(1));  f1_ = f1_(rr_, cc_, :);  f2_ = double(ld_(2));  f2_ = f2_(rr_, cc_, :);  f10_ = double(ld_(10));  f10_ = f10_(rr_, cc_, :);
dm_ = @(a, b) sum(abs(a - b), 3);
figure('Name', sprintf('Flat-water bias CoreView %d, around (%g, %g) cm', P.run_num, P.diag_xy), 'Position', [60 60 1500 520]);
pan_ = {uint8(f1_), uint8(f10_), dm_(f2_, f1_), dm_(f10_, f1_), dm_(f10_, f2_)};
ttl_ = {'frame 1', 'frame 10', '|f2 - f1| (R+G+B)', '|f10 - f1|', '|f10 - f2|'};
cl_ = [0, max(5, prctile(pan_{4}(:), 99.5))];
for k_ = 1:5
    subplot(1, 5, k_);
    if k_ <= 2, image(pan_{k_}); else, imagesc(pan_{k_}, cl_); colorbar; end
    axis image;  title(ttl_{k_});
    hold on;  plot(uv_(1) - cc_(1) + 1, uv_(2) - rr_(1) + 1, 'c+', 'MarkerSize', 12);
end
clear R_ ref_b_ setup_b_ rd_ ok_ ld_ n_ sx_ sy_ ob_ dc_ c_ d_ y0_ uv_ hw_ rr_ cc_ f1_ f2_ f10_ dm_ pan_ ttl_ cl_ k_

%% 5e. Test: per-frame brightness correction 'none' | 'offset' | 'gain' (~3 min, nothing saved, no bias subtraction)
% Flat frames still flicker by ~0.15 counts (5d), which reads as slope wherever the colour map is flat or has a
% sharp feature (the patch at (-10, 7) and the target edges): such an artefact changes sign with the flicker, so
% the fixed bias map cannot remove it.  Per frame, 'offset' removes median(obs - ref) per channel.
P.drift_frames = [10 100 200 2000];
P.drift_modes  = {'none', 'offset', 'gain'};
fprintf('\n--- 5e. brightness correction test ---\n');
R_ = cisg_setup_run(struct('coreview_num', P.run_num, 'roi_mode', P.video.roi_mode, 'dye_model', P.video.dye_model, ...
                           'bin', P.video.bin, 'ref_frames', P.video.ref_frames));
[ref_b_, setup_b_] = cisg_bin_setup(R_.ref_image, R_.setup, P.video.bin);
sx10_ = cell(1, numel(P.drift_modes));
fprintf('frame   mode     mean Sx   mean Sy   rms Sx   rms Sy   masked\n');
for n_ = P.drift_frames
    for m_ = 1:numel(P.drift_modes)
        st_ = setup_b_;  st_.drift_corr = P.drift_modes{m_};
        [sx_, sy_] = cisg_cache_frame(n_, ref_b_, st_, R_.frame_file, '', P.video.bin, R_.dt, R_.img_size);
        fprintf('%5d   %-6s  %+8.4f  %+8.4f  %7.4f  %7.4f  %6.1f%%\n', n_, P.drift_modes{m_}, mean(sx_, 'all', 'omitnan'), ...
                mean(sy_, 'all', 'omitnan'), rms(sx_, 'all', 'omitnan'), rms(sy_, 'all', 'omitnan'), 100 * mean(isnan(sx_), 'all'));
        if n_ == P.drift_frames(1), sx10_{m_} = sx_; end
    end
end
figure('Name', sprintf('Brightness correction, CoreView %d frame %d', P.run_num, P.drift_frames(1)), 'Position', [60 60 1500 450]);
lim_ = 3 * rms(sx10_{1}, 'all', 'omitnan');
for m_ = 1:numel(P.drift_modes)
    subplot(1, numel(P.drift_modes), m_);
    imagesc(sx10_{m_}, 'AlphaData', isfinite(sx10_{m_}), [-lim_ lim_]);  axis image;  colorbar;  set(gca, 'Color', [0.8 0.8 0.8]);
    title(sprintf('S_x frame %d, drift\\_corr = %s (rms %.4f)', P.drift_frames(1), P.drift_modes{m_}, rms(sx10_{m_}, 'all', 'omitnan')));
end
colormap(interp1([0 0.5 1], [0.2 0.3 0.8; 1 1 1; 0.8 0.2 0.2], linspace(0, 1, 256)));
clear R_ ref_b_ setup_b_ sx10_ n_ m_ st_ sx_ sy_ lim_

%% 5f. Diagnostic: why the top-left does not converge on flat water (~1 min, nothing saved)
% On a flat frame the shift is ~0, so a pixel only fails to converge where the R/G colour map folds (|det J| ~ 0).
% Two candidate causes: (1) a channel clipped at 255 (no gradient), (2) the edge light falloff times the colour
% ramp, which makes a channel peak a few cm inside the edge.  Also shown: how much the blue channel would add
% (smallest singular value of the 3x2 RGB Jacobian vs the 2x2 RG one).
P.diag_frame = 10;
P.sat_level  = 250;                                % counts (8-bit): at or above = clipped
fprintf('\n--- 5f. saturation / colour-fold diagnostic, frame %d ---\n', P.diag_frame);
R_ = cisg_setup_run(struct('coreview_num', P.run_num, 'roi_mode', P.video.roi_mode, 'dye_model', P.video.dye_model, ...
                           'bin', P.video.bin, 'ref_frames', P.video.ref_frames, 'setup', struct('drift_corr', P.video.drift_corr, 'colour_channels', P.video.colour_channels)));
[ref_b_, setup_b_] = cisg_bin_setup(R_.ref_image, R_.setup, P.video.bin);
obs_ = cisg_load_coreview(R_.frame_file(P.diag_frame), R_.img_size);
if P.video.bin > 1, obs_ = uint8(imresize(obs_, 1 / P.video.bin, 'box')); end
[~, ~, ~, D_] = cisg_colour_displacement(ref_b_, obs_, setup_b_);
r_ = setup_b_.roi;  rr_ = r_.row_min:r_.row_max;  cc_ = r_.col_min:r_.col_max;
nc_ = bitand(D_.bits, 8) > 0;  ill_ = bitand(D_.bits, 4) > 0;  dark_ = bitand(D_.bits, 64) > 0;

% (1) saturation, from the full-resolution reference (binning would average it away)
sat3_ = R_.ref_image >= P.sat_level;
satb_ = imresize(double(any(sat3_, 3)), [size(ref_b_, 1) size(ref_b_, 2)], 'box') > 0;   % any clipped px in the bin
satb_ = satb_(rr_, cc_);
mx_ = max(ref_b_(rr_, cc_, :), [], 3);
on_ = ~dark_;
fprintf('reference at >= %d counts (whole image, full res): R %.2f%%, G %.2f%%, B %.2f%%\n', P.sat_level, ...
        arrayfun(@(k) 100 * mean(reshape(sat3_(:, :, k), [], 1)), 1:3));
near_ = imdilate(satb_, strel('disk', 3));
fprintf('not-converged pixels: %.2f%% of on-target; %.0f%% of them within 3 px of a clipped pixel\n', ...
        100 * nnz(nc_ & on_) / max(nnz(on_), 1), 100 * nnz(nc_ & near_) / max(nnz(nc_), 1));

% (2) colour-map folds: smoothed reference as the solver sees it, and its gradients
C_ = imgaussfilt(single(ref_b_) / 255, setup_b_.colour_sigma);
gx_ = cell(1, 3);  gy_ = cell(1, 3);
for k_ = 1:3, [gx_{k_}, gy_{k_}] = gradient(C_(:, :, k_));  gx_{k_} = gx_{k_}(rr_, cc_);  gy_{k_} = gy_{k_}(rr_, cc_); end
% smallest singular value of J (rows = channels, columns = d/dx, d/dy): sqrt of the smaller eigenvalue of J'J
smin_ = @(ch) sqrt(max(0, (sum(cat(3, gx_{ch}) .^ 2, 3) + sum(cat(3, gy_{ch}) .^ 2, 3)) / 2 ...
              - sqrt(((sum(cat(3, gx_{ch}) .^ 2, 3) - sum(cat(3, gy_{ch}) .^ 2, 3)) / 2) .^ 2 + sum(cat(3, gx_{ch}) .* cat(3, gy_{ch}), 3) .^ 2)));
sRG_ = smin_([1 2]);  sRGB_ = smin_([1 2 3]);
med_ = median(sRG_(on_), 'omitnan');
fprintf('smallest singular value of J (/ on-target median of RG): not-converged px median RG %.3f, RGB %.3f (blue gains x%.1f)\n', ...
        median(sRG_(nc_)) / med_, median(sRGB_(nc_)) / med_, median(sRGB_(nc_) ./ max(sRG_(nc_), eps)));
% which gradient vanishes there: share of not-converged px where each component is below 10% of its on-target median
nm_ = {'dR/dx', 'dR/dy', 'dG/dx', 'dG/dy'};  gg_ = {gx_{1}, gy_{1}, gx_{2}, gy_{2}};
for k_ = 1:4
    fprintf('  %-6s small (< 10%% of median) at %3.0f%% of not-converged px (%3.0f%% of all on-target px)\n', nm_{k_}, ...
            100 * mean(abs(gg_{k_}(nc_)) < 0.1 * median(abs(gg_{k_}(on_)))), 100 * mean(abs(gg_{k_}(on_)) < 0.1 * median(abs(gg_{k_}(on_)))));
end

% profiles across the two arms of the L: the column (left 40%) and the row (top 40%) with the most failures
w_ = numel(cc_);  h_ = numel(rr_);
[~, ic_] = max(sum(nc_(:, 1:round(0.4 * w_)), 1));    % vertical arm: column with most failures
[~, ir_] = max(sum(nc_(1:round(0.4 * h_), :), 2));    % horizontal arm: row with most failures
prow_ = round(median(find(nc_(:, ic_))));             % horizontal profile through the vertical arm
pcol_ = round(median(find(nc_(ir_, :))));             % vertical profile through the horizontal arm
if ~isfinite(prow_), prow_ = round(h_ / 2); end                 % no failures: profiles through the middle
if ~isfinite(pcol_), pcol_ = round(w_ / 2); end

figure('Name', sprintf('Colour fold diagnostic, CoreView %d frame %d', P.run_num, P.diag_frame), 'Position', [30 30 1800 900]);
tl_ = tiledlayout(2, 4, 'TileSpacing', 'compact');
[yy_, xx_] = find(nc_);
nexttile(tl_);  image(ref_b_(rr_, cc_, :));  axis image;  hold on;
plot(xx_, yy_, 'r.', 'MarkerSize', 2);  contour(double(satb_), [0.5 0.5], 'w', 'LineWidth', 1);
plot([1 w_], [prow_ prow_], 'c-', [pcol_ pcol_], [1 h_], 'm-');
title('reference; red = not converged, white = clipped, cyan / magenta = profiles');
nexttile(tl_);  imagesc(mx_, [150 255]);  axis image;  colorbar;  hold on;  contour(double(satb_), [0.5 0.5], 'w');
title(sprintf('brightest channel (counts); white contour = >= %d', P.sat_level));
nexttile(tl_);  imagesc(log10(sRG_ / med_), [-2 0.5]);  axis image;  colorbar;  hold on;  plot(xx_, yy_, 'r.', 'MarkerSize', 1);
title('log_{10} \sigma_{min}(J_{RG}) / median: folds are dark');
nexttile(tl_);  imagesc(log10(sRGB_ / med_), [-2 0.5]);  axis image;  colorbar;  hold on;  plot(xx_, yy_, 'r.', 'MarkerSize', 1);
title('log_{10} \sigma_{min}(J_{RGB}) / median: with blue');
for k_ = 1:4
    nexttile(tl_);  imagesc(sign(gg_{k_}) .* min(abs(gg_{k_}) / median(abs(gg_{k_}(on_))), 1), [-1 1]);  axis image;  hold on;
    plot(xx_, yy_, 'k.', 'MarkerSize', 1);
    title(sprintf('%s (sign, clipped at median): zero line = fold', nm_{k_}));
    colormap(gca, interp1([0 0.5 1], [0.2 0.3 0.8; 1 1 1; 0.8 0.2 0.2], linspace(0, 1, 256)));
end
figure('Name', sprintf('Colour profiles across the L, CoreView %d', P.run_num), 'Position', [60 60 1400 450]);
ch_ = {'r', [0 0.6 0], 'b'};
subplot(1, 2, 1);  hold on;
for k_ = 1:3, plot(cc_ - cc_(1) + 1, 255 * C_(rr_(prow_), cc_, k_), 'Color', ch_{k_}); end
xline(ic_, 'k--', 'most failures');  yline(P.sat_level, 'k:');
xlabel('column (binned, ROI)');  ylabel('smoothed counts');  title(sprintf('across the vertical arm (row %d)', prow_));  grid on;
subplot(1, 2, 2);  hold on;
for k_ = 1:3, plot(rr_ - rr_(1) + 1, 255 * C_(rr_, cc_(pcol_), k_), 'Color', ch_{k_}); end
xline(ir_, 'k--', 'most failures');  yline(P.sat_level, 'k:');
xlabel('row (binned, ROI)');  ylabel('smoothed counts');  title(sprintf('across the horizontal arm (column %d)', pcol_));  grid on;
legend({'R', 'G', 'B'});
clear R_ ref_b_ setup_b_ obs_ D_ r_ rr_ cc_ nc_ ill_ dark_ sat3_ satb_ mx_ on_ near_ C_ gx_ gy_ k_ smin_ sRG_ sRGB_ med_ ...
      nm_ gg_ w_ h_ ic_ ir_ prow_ pcol_ tl_ yy_ xx_ ch_

%% 5g. Compare runs: is the top-left non-convergence in run 101 only? (~2 min, nothing saved)
% The colour solve uses only the reference and the frame (the plate cal sets cm/px, L_eff and the origin, and enters
% only the sigma_S test), so the same mask on another run with the same target means it is the target / lighting.
% Dye solve off everywhere, so all runs use the same 2-channel solver.
P.cmp_runs  = [100 101 96];                        % run 96: low rhodamine
P.cmp_plate = [106 106 97];                        % plate cal per run (97 has no laser line; does not affect convergence)
P.cmp_frame = 10;                                  % flat (before waves) in each run?  check the printed frame times
names_ = {'dark in reference', 'shift out of image', 'shift clamped', 'ill-conditioned J', 'not converged', ...
          'colour residual', 'precision sigma_S > max'};
vals_  = [64 1 2 4 8 16 128];
fprintf('\n--- 5g. run comparison, frame %d ---\n', P.cmp_frame);
figure('Name', sprintf('Mask reasons by run, frame %d', P.cmp_frame), 'Position', [30 30 1800 800]);
tl_ = tiledlayout(2, numel(P.cmp_runs), 'TileSpacing', 'compact');
for k_ = 1:numel(P.cmp_runs)
    run_ = P.cmp_runs(k_);
    R_ = cisg_setup_run(struct('coreview_num', run_, 'roi_mode', P.video.roi_mode, 'dye_model', 'none', 'bin', P.video.bin, ...
                               'ref_frames', P.video.ref_frames, 'setup', struct('drift_corr', P.video.drift_corr, 'colour_channels', P.video.colour_channels), ...
                               'cal_path', fullfile(P.here, 'calibration', sprintf('cisg_cal_CoreView_%d.mat', P.cmp_plate(k_)))));
    [ref_b_, setup_b_] = cisg_bin_setup(R_.ref_image, R_.setup, P.video.bin);
    fprintf('run %d: frame %d is %.1f s after wind start (NaN = no wind clock)\n', run_, P.cmp_frame, R_.t_wind(P.cmp_frame));
    obs_ = cisg_load_coreview(R_.frame_file(P.cmp_frame), R_.img_size);
    if P.video.bin > 1, obs_ = uint8(imresize(obs_, 1 / P.video.bin, 'box')); end
    [~, ~, ~, D_] = cisg_colour_displacement(ref_b_, obs_, setup_b_);
    first_ = zeros(size(D_.bits));
    for j_ = numel(vals_):-1:1, first_(bitand(D_.bits, vals_(j_)) > 0) = j_; end
    fr_ = arrayfun(@(v) 100 * mean(bitand(D_.bits(:), v) > 0), vals_);
    fprintf('run %d (plate %d): masked %.1f%% |', run_, P.cmp_plate(k_), 100 * mean(D_.bits(:) > 0));
    fprintf(' %s %.1f', 'dark', fr_(1), 'oob', fr_(2), 'clamp', fr_(3), 'ill-c', fr_(4), 'noconv', fr_(5), 'resid', fr_(6), 'sigma', fr_(7));
    fprintf('\n');
    r_ = setup_b_.roi;
    nexttile(tl_, k_);  image(ref_b_(r_.row_min:r_.row_max, r_.col_min:r_.col_max, :));  axis image;
    title(sprintf('run %d reference (median of frames %d-%d)', run_, min(P.video.ref_frames), max(P.video.ref_frames)));
    nexttile(tl_, k_ + numel(P.cmp_runs));  imagesc(first_, [-0.5 numel(vals_) + 0.5]);  axis image;
    colormap(gca, [1 1 1; lines(numel(vals_))]);
    title(sprintf('run %d frame %d: not converged %.1f%%, ill-cond %.1f%%', run_, P.cmp_frame, fr_(5), fr_(4)));
end
colorbar(nexttile(tl_, 2 * numel(P.cmp_runs)), 'Ticks', 0:numel(vals_), 'TickLabels', [{'kept'} names_]);
clear names_ vals_ tl_ k_ run_ R_ ref_b_ setup_b_ obs_ D_ first_ j_ fr_ r_

%% 5h. Test: colour solve with R, G vs R, G, B (~2 min, nothing saved, no bias subtraction)
% 5f: where the R/G map folds (red peaks at the left, green at the top of the target) blue lifts sigma_min(J) ~4x.
% Frame 10 (flat): the L of not-converged pixels should go; frame 2000 (waves): masked area and agreement.
P.rgb_frames = [10 2000];
P.rgb_cases  = {'RG', 'RGB', 'hybrid'};
names_ = {'dark', 'oob', 'clamp', 'ill-c', 'noconv', 'resid', 'sigma'};
vals_  = [64 1 2 4 8 16 128];
fprintf('\n--- 5h. R,G vs R,G,B colour solve ---\n');
R_ = cisg_setup_run(struct('coreview_num', P.run_num, 'roi_mode', P.video.roi_mode, 'dye_model', P.video.dye_model, ...
                           'bin', P.video.bin, 'ref_frames', P.video.ref_frames, 'setup', struct('drift_corr', P.video.drift_corr)));
[ref_b_, setup_b_] = cisg_bin_setup(R_.ref_image, R_.setup, P.video.bin);
nf_ = numel(P.rgb_frames);  nc_ = numel(P.rgb_cases);  sxs_ = cell(nf_, nc_);  sys_ = sxs_;  fst_ = cell(nf_, nc_);
fprintf('frame  channels  masked |%s|  rms Sx  rms Sy     s\n', sprintf(' %6s', names_{:}));
for f_ = 1:nf_
    n_ = P.rgb_frames(f_);
    obs_ = cisg_load_coreview(R_.frame_file(n_), R_.img_size);
    if P.video.bin > 1, obs_ = uint8(imresize(obs_, 1 / P.video.bin, 'box')); end
    for c_ = 1:nc_
        st_ = setup_b_;  st_.colour_channels = P.rgb_cases{c_};
        [~, ~, ci_, D_] = cisg_colour_displacement(ref_b_, obs_, st_);
        [sx_, sy_] = cisg_cache_frame(n_, ref_b_, st_, R_.frame_file, '', P.video.bin, R_.dt, R_.img_size);
        sxs_{f_, c_} = sx_;  sys_{f_, c_} = sy_;
        first_ = zeros(size(D_.bits));
        for j_ = numel(vals_):-1:1, first_(bitand(D_.bits, vals_(j_)) > 0) = j_; end
        fst_{f_, c_} = first_;
        fr_ = arrayfun(@(v) 100 * mean(bitand(D_.bits(:), v) > 0), vals_);
        fprintf('%5d  %-8s  %5.1f%% |%s|  %6.4f  %6.4f  %4.0f', n_, P.rgb_cases{c_}, 100 * mean(isnan(sx_), 'all'), ...
                sprintf(' %5.1f%%', fr_), rms(sx_, 'all', 'omitnan'), rms(sy_, 'all', 'omitnan'), ci_.seconds);
        if strcmpi(P.rgb_cases{c_}, 'hybrid'), fprintf('   (%.1f%% of px from R,G,B)', 100 * ci_.hybrid_fraction); end
        fprintf('\n');
    end
    for c_ = 2:nc_
        k_ = isfinite(sxs_{f_, 1}) & isfinite(sxs_{f_, c_});
        fprintf('       pixels both kept: rms Sx difference %s - RG %.4f\n', P.rgb_cases{c_}, rms(sxs_{f_, c_}(k_) - sxs_{f_, 1}(k_)));
    end
end
% Physics check: a single-valued surface has curl(S) = dSx/dy - dSy/dx = 0 exactly (S = grad eta); measurement errors
% generally do not.  Score = rms(curl) / rms(div) on 4x4-px blocks (div = curvature, the real signal): 0 = perfectly
% integrable, ~1 = like random noise.  Lower = more physical.  Regions: each case's own kept pixels; pixels all cases
% kept (fair comparison); steep (|S| > 0.2 from RG); disputed (all kept, the cases' Sx differ by > 0.02).
% Sy is d eta / dy with y = image up (cisg_slope_to_eta), so d/dy = -d/d(row).
qb_ = 4;  mmx_ = 10 * setup_b_.cm_per_pixel_x;  mmy_ = 10 * setup_b_.cm_per_pixel_y;
fprintf('\ncurl / div (lower = more physical; 0 = integrable, ~1 = noise), %dx%d px blocks\n', qb_, qb_);
fprintf('frame  channels    own kept      all kept      steep         disputed\n');
for f_ = 1:nf_
    [nr_, nc2_] = size(sxs_{f_, 1});  nr_ = floor(nr_ / qb_);  nc2_ = floor(nc2_ / qb_);
    blk_ = @(A) squeeze(mean(reshape(double(A(1:nr_ * qb_, 1:nc2_ * qb_)), qb_, nr_, qb_, nc2_), [1 3]));   % NaN if any px is NaN
    com_ = true(size(sxs_{f_, 1}));  spr_ = zeros(size(com_));
    for c_ = 1:nc_, com_ = com_ & isfinite(sxs_{f_, c_}) & isfinite(sys_{f_, c_}); end
    for c_ = 2:nc_, spr_ = max(spr_, abs(sxs_{f_, c_} - sxs_{f_, 1})); end
    reg_ = {[], blk_(com_) > 0.5, blk_(com_ & hypot(sxs_{f_, 1}, sys_{f_, 1}) > 0.2) > 0.5, blk_(com_ & spr_ > 0.02) > 0.5};
    for c_ = 1:nc_
        sxm_ = sxs_{f_, c_};  sym_ = sys_{f_, c_};
        fprintf('%5d  %-8s', P.rgb_frames(f_), P.rgb_cases{c_});
        for g_ = 1:numel(reg_)
            if g_ > 1, sxm_(~com_) = NaN;  sym_(~com_) = NaN; end
            [ax_, ar_] = gradient(blk_(sxm_), qb_ * mmx_, qb_ * mmy_);
            [bx_, br_] = gradient(blk_(sym_), qb_ * mmx_, qb_ * mmy_);
            cu_ = -ar_ - bx_;  dv_ = ax_ - br_;                       % curl, div with y up
            m_ = isfinite(cu_) & isfinite(dv_);  if g_ > 1, m_ = m_ & reg_{g_}; end
            fprintf('  %5.3f (%5d)', rms(cu_(m_)) / rms(dv_(m_)), nnz(m_));
        end
        fprintf('\n');
    end
end
fprintf('(block counts in brackets)\n');
figure('Name', sprintf('R,G vs R,G,B, CoreView %d', P.run_num), 'Position', [30 30 1700 900]);
tl_ = tiledlayout(nf_, 2 * nc_, 'TileSpacing', 'compact');
for f_ = 1:nf_
    lim_ = 3 * rms(sxs_{f_, 1}, 'all', 'omitnan');
    for c_ = 1:nc_
        nexttile(tl_);
        imagesc(sxs_{f_, c_}, 'AlphaData', isfinite(sxs_{f_, c_}), [-lim_ lim_]);  axis image;  set(gca, 'Color', [0.8 0.8 0.8]);
        colormap(gca, interp1([0 0.5 1], [0.2 0.3 0.8; 1 1 1; 0.8 0.2 0.2], linspace(0, 1, 256)));  colorbar;
        title(sprintf('S_x frame %d, %s (rms %.4f)', P.rgb_frames(f_), P.rgb_cases{c_}, rms(sxs_{f_, c_}, 'all', 'omitnan')));
    end
    for c_ = 1:nc_
        nexttile(tl_);  imagesc(fst_{f_, c_}, [-0.5 numel(vals_) + 0.5]);  axis image;
        colormap(gca, [1 1 1; lines(numel(vals_))]);
        cb_ = colorbar;  cb_.Ticks = 0:numel(vals_);  cb_.TickLabels = [{'kept'} names_];
        title(sprintf('mask reason frame %d, %s', P.rgb_frames(f_), P.rgb_cases{c_}));
    end
end
clear names_ vals_ R_ ref_b_ setup_b_ nf_ nc_ sxs_ sys_ fst_ f_ n_ obs_ c_ st_ ci_ D_ sx_ sy_ first_ j_ fr_ k_ tl_ lim_ cb_ ...
      qb_ mmx_ mmy_ nr_ nc2_ blk_ com_ spr_ reg_ sxm_ sym_ g_ ax_ ar_ bx_ br_ cu_ dv_ m_

%% 5i. Diagnostic: why R,G,B disagrees with R,G under waves (frame 2000, ~1 min, nothing saved)
% 5h: on frame 2000 RGB doubled 'not converged' and differs from RG by ~0.035 rms in Sx where both kept the pixel.
% (a) Is blue consistent with the shift R and G agree on?  Residual of each channel at RG's converged solution
%     (blue is not in that fit, so its residual tests it) vs the camera noise.
% (b) Does blue move by a different amount (dispersion of water, lens chromatic aberration)?  If blue's shift is
%     (1 + kappa) x the R,G shift, blue's residual is kappa * grad(B) . delta: fit kappa.  Water dispersion ~ +1%.
% (c) Where do the extra not-converged pixels of RGB sit: steep wave faces (large |S|)?
P.res_frame = 2000;
fprintf('\n--- 5i. RGB vs RG under waves, frame %d ---\n', P.res_frame);
R_ = cisg_setup_run(struct('coreview_num', P.run_num, 'roi_mode', P.video.roi_mode, 'dye_model', P.video.dye_model, ...
                           'bin', P.video.bin, 'ref_frames', P.video.ref_frames, 'setup', struct('drift_corr', P.video.drift_corr)));
[ref_b_, setup_b_] = cisg_bin_setup(R_.ref_image, R_.setup, P.video.bin);
obs_ = cisg_load_coreview(R_.frame_file(P.res_frame), R_.img_size);
if P.video.bin > 1, obs_ = uint8(imresize(obs_, 1 / P.video.bin, 'box')); end
stg_ = setup_b_;  stg_.colour_channels = 'RG';
stb_ = setup_b_;  stb_.colour_channels = 'RGB';
[~, ~, ~, Dg_] = cisg_colour_displacement(ref_b_, obs_, stg_);
[~, ~, ~, Db_] = cisg_colour_displacement(ref_b_, obs_, stb_);
[sxg_, syg_] = cisg_cache_frame(P.res_frame, ref_b_, stg_, R_.frame_file, '', P.video.bin, R_.dt, R_.img_size);
[sxb_, ~]    = cisg_cache_frame(P.res_frame, ref_b_, stb_, R_.frame_file, '', P.video.bin, R_.dt, R_.img_size);
if ~isequal(size(sxg_), size(Dg_.bits)), error('5i: slope map %s and diagnostics %s differ in size', mat2str(size(sxg_)), mat2str(size(Dg_.bits))); end
r_ = setup_b_.roi;  rr_ = r_.row_min:r_.row_max;  cc_ = r_.col_min:r_.col_max;
cvg_ = Dg_.bits == 0;  cvb_ = Db_.bits == 0;          % converged and passed every test (before despeckle)

% (a) residuals at the converged pixels, counts (x 255), against the camera noise of C_obs - C_ref
nz_ = 255 * (setup_b_.noise.a(:).' + 0.4 * setup_b_.noise.b(:).');   % a, b already rescaled for the median reference
q_ = @(v) [median(abs(v)), prctile(abs(v), 90), median(v)];
fprintf('camera noise of C_obs - C_ref at C = 0.4 (counts): R %.2f  G %.2f  B %.2f\n', nz_);
fprintf('residual (counts) at converged px:   |median|  90%%   signed median\n');
eg_ = {Dg_.e(:, :, 1), Dg_.e(:, :, 2), Dg_.eB};  eb_ = {Db_.e(:, :, 1), Db_.e(:, :, 2), Db_.e(:, :, 3)};
ch_ = 'RGB';
for k_ = 1:3, fprintf('  RG solve,  %s%s   %6.2f  %6.2f  %+6.2f\n', ch_(k_), repmat(' (not in fit)', 1, k_ == 3), 255 * q_(eg_{k_}(cvg_))); end
for k_ = 1:3, fprintf('  RGB solve, %s               %6.2f  %6.2f  %+6.2f\n', ch_(k_), 255 * q_(eb_{k_}(cvb_))); end

% (b) dispersion / chromatic aberration: blue residual at the RG solution vs grad(B) . delta
CB_ = imgaussfilt(single(ref_b_(:, :, 3)) / 255, setup_b_.colour_sigma);
[gBx_, gBy_] = gradient(CB_);
[Xg_, Yg_] = meshgrid(single(cc_), single(rr_));
xs_ = Xg_ + Dg_.dx;  ys_ = Yg_ + Dg_.dy;
pr_ = interp2(gBx_, xs_, ys_, 'linear') .* Dg_.dx + interp2(gBy_, xs_, ys_, 'linear') .* Dg_.dy;   % grad(B) . delta
m_ = cvg_ & isfinite(pr_) & isfinite(Dg_.eB);
kap_ = sum(double(Dg_.eB(m_)) .* double(pr_(m_))) / sum(double(pr_(m_)) .^ 2);
expl_ = 1 - var(double(Dg_.eB(m_)) - kap_ * double(pr_(m_))) / var(double(Dg_.eB(m_)));
fprintf('blue shift = (1 + kappa) x RG shift: kappa = %+.1f%%, explains %.0f%% of the blue residual variance (water dispersion ~ +1%%)\n', ...
        100 * kap_, 100 * expl_);
fprintf('median |shift| at converged px: %.1f px (binned)\n', median(hypot(Dg_.dx(cvg_), Dg_.dy(cvg_))));

% (c) extra not-converged pixels: converged with RG, not with RGB
xtra_ = cvg_ & bitand(Db_.bits, 8) > 0;
Sg_ = hypot(sxg_, syg_);
fprintf('extra not-converged with RGB: %.1f%% of the ROI; median |S| there %.3f vs %.3f over all RG-kept px\n', ...
        100 * mean(xtra_(:)), median(Sg_(xtra_), 'omitnan'), median(Sg_(isfinite(Sg_)), 'omitnan'));
fprintf('  blue residual at the RG solution there: |median| %.2f counts vs %.2f elsewhere\n', ...
        255 * median(abs(Dg_.eB(xtra_))), 255 * median(abs(Dg_.eB(cvg_ & ~xtra_))));

figure('Name', sprintf('RGB vs RG under waves, CoreView %d frame %d', P.run_num, P.res_frame), 'Position', [30 30 1750 900]);
tl_ = tiledlayout(2, 3, 'TileSpacing', 'compact');
cm_ = interp1([0 0.5 1], [0.2 0.3 0.8; 1 1 1; 0.8 0.2 0.2], linspace(0, 1, 256));
[yy_, xx_] = find(xtra_);
nexttile(tl_);  image(obs_(rr_, cc_, :));  axis image;  hold on;  plot(xx_, yy_, 'r.', 'MarkerSize', 1);
title(sprintf('frame %d; red = not converged only with RGB (%.1f%%)', P.res_frame, 100 * mean(xtra_(:))));
lim_ = 3 * rms(sxg_, 'all', 'omitnan');
nexttile(tl_);  imagesc(sxg_, 'AlphaData', isfinite(sxg_), [-lim_ lim_]);  axis image;  colormap(gca, cm_);  colorbar;
set(gca, 'Color', [0.8 0.8 0.8]);  title('S_x, R,G solve');
d_ = sxb_ - sxg_;  ld_ = 3 * rms(d_, 'all', 'omitnan');
nexttile(tl_);  imagesc(d_, 'AlphaData', isfinite(d_), [-ld_ ld_]);  axis image;  colormap(gca, cm_);  colorbar;
set(gca, 'Color', [0.8 0.8 0.8]);  title(sprintf('S_x RGB - RG (rms %.4f)', rms(d_, 'all', 'omitnan')));
eBc_ = 255 * Dg_.eB;  eBc_(~cvg_) = NaN;
nexttile(tl_);  imagesc(eBc_, 'AlphaData', isfinite(eBc_), [-5 5]);  axis image;  colormap(gca, cm_);  colorbar;
set(gca, 'Color', [0.8 0.8 0.8]);  title('blue residual at the R,G solution (counts)');
nexttile(tl_);  ed_ = linspace(prctile(pr_(m_), 1), prctile(pr_(m_), 99), 30);
[~, ~, bi_] = histcounts(pr_(m_), ed_);  ok_ = bi_ > 0;
ev_ = double(Dg_.eB(m_));  mb_ = accumarray(bi_(ok_), ev_(ok_), [numel(ed_) - 1 1], @mean, NaN);
ce_ = (ed_(1:end-1) + ed_(2:end)) / 2;
plot(255 * ce_, 255 * mb_, 'b.-', 255 * ce_, 255 * kap_ * ce_, 'k--');  grid on;
xlabel('grad(B) . \delta (counts)');  ylabel('mean blue residual (counts)');
title(sprintf('blue vs R,G shift: kappa = %+.1f%%', 100 * kap_));  legend({'binned', 'fit'}, 'Location', 'best');
nexttile(tl_);  se_ = 0:0.05:0.6;
[~, ~, bs_] = histcounts(Sg_(cvg_), se_);  xv_ = xtra_(cvg_);  okb_ = bs_ > 0;
fx_ = accumarray(bs_(okb_), double(xv_(okb_)), [numel(se_) - 1 1], @mean, NaN);
bar((se_(1:end-1) + se_(2:end)) / 2, 100 * fx_);  grid on;
xlabel('|S| from the R,G solve');  ylabel('% not converged with RGB');  title('extra RGB failures vs slope magnitude');
clear R_ ref_b_ setup_b_ obs_ stg_ stb_ Dg_ Db_ sxg_ syg_ sxb_ r_ rr_ cc_ cvg_ cvb_ nz_ q_ eg_ eb_ ch_ k_ CB_ gBx_ gBy_ ...
      Xg_ Yg_ xs_ ys_ pr_ m_ kap_ expl_ xtra_ Sg_ tl_ cm_ yy_ xx_ lim_ d_ ld_ eBc_ ed_ bi_ ok_ mb_ ev_ ce_ se_ bs_ xv_ okb_ fx_

%% 5j. Check: are Sx and Sy consistent (sign, scale)?  curl / div vs block size (frame 2000, ~1 min, nothing saved)
% 5h gave curl / div ~0.75 for every solve.  Either fine-scale noise (differentiating at 4 px amplifies it), which
% falls with block size, or an Sx / Sy convention error: with Sy of the wrong sign, or scaled by b ~= 1, even a
% perfect surface has curl ((1 - b) d2eta/dxdy).  Per block size: score as processed, with Sy flipped, and with the
% Sy factor b that minimises the curl (b ~ +1: consistent; b ~ -1: Sy sign flipped; else a relative scale error).
% Blocks need >= 75% valid pixels.  Sy is d eta / dy with y = image up, so d/dy = -d/d(row).
P.conv_frame = 2000;
P.conv_q     = [4 8 16 32];
fprintf('\n--- 5j. Sx / Sy convention check, frame %d (%s) ---\n', P.conv_frame, P.video.colour_channels);
R_ = cisg_setup_run(struct('coreview_num', P.run_num, 'roi_mode', P.video.roi_mode, 'dye_model', P.video.dye_model, ...
                           'bin', P.video.bin, 'ref_frames', P.video.ref_frames, ...
                           'setup', struct('drift_corr', P.video.drift_corr, 'colour_channels', P.video.colour_channels)));
[ref_b_, setup_b_] = cisg_bin_setup(R_.ref_image, R_.setup, P.video.bin);
[sx_, sy_] = cisg_cache_frame(P.conv_frame, ref_b_, setup_b_, R_.frame_file, '', P.video.bin, R_.dt, R_.img_size);
mmx_ = 10 * setup_b_.cm_per_pixel_x;  mmy_ = 10 * setup_b_.cm_per_pixel_y;
fprintf('pixel size %.4f x %.4f mm (x, y, binned); rms Sx %.4f, Sy %.4f\n', mmx_, mmy_, rms(sx_, 'all', 'omitnan'), rms(sy_, 'all', 'omitnan'));
fprintf('block  (mm)    blocks   as processed   Sy flipped   best Sy factor b   score at b\n');
sc_ = @(cu, dv) rms(cu) / rms(dv);
res_ = nan(numel(P.conv_q), 4);
for j_ = 1:numel(P.conv_q)
    q_ = P.conv_q(j_);
    [nr_, nc_] = size(sx_);  nr_ = floor(nr_ / q_);  nc_ = floor(nc_ / q_);
    rs_ = @(A) reshape(double(A(1:nr_ * q_, 1:nc_ * q_)), q_, nr_, q_, nc_);
    v_  = squeeze(mean(rs_(isfinite(sx_) & isfinite(sy_)), [1 3]));            % valid fraction per block
    bx_ = squeeze(mean(rs_(sx_), [1 3], 'omitnan'));  by_ = squeeze(mean(rs_(sy_), [1 3], 'omitnan'));
    bx_(v_ < 0.75) = NaN;  by_(v_ < 0.75) = NaN;
    [axc_, axr_] = gradient(bx_, q_ * mmx_, q_ * mmy_);                          % d/dcol, d/drow
    [ayc_, ayr_] = gradient(by_, q_ * mmx_, q_ * mmy_);
    A_ = -axr_;  B_ = -ayc_;                                                     % curl(b) = dSx/dy - b dSy/dx = A + b B
    m_ = isfinite(A_) & isfinite(B_) & isfinite(axc_) & isfinite(ayr_);
    s0_ = sc_(A_(m_) + B_(m_), axc_(m_) - ayr_(m_));                             % as processed (b = 1)
    s1_ = sc_(A_(m_) - B_(m_), axc_(m_) + ayr_(m_));                             % Sy flipped (b = -1)
    b_  = -sum(A_(m_) .* B_(m_)) / sum(B_(m_) .^ 2);                             % least-squares b
    sb_ = sc_(A_(m_) + b_ * B_(m_), axc_(m_) - b_ * ayr_(m_));
    res_(j_, :) = [s0_ s1_ b_ sb_];
    fprintf('%4d  (%5.2f)  %7d      %6.3f        %6.3f          %+6.3f          %6.3f\n', q_, q_ * mmx_, nnz(m_), s0_, s1_, b_, sb_);
end
fprintf(['reading: score falling with block size and b ~ +1 at the larger blocks -> fine-scale noise, conventions OK;\n' ...
         '         flipped much lower or b ~ -1 -> Sy sign error; b far from 1 at every size -> Sx / Sy scale mismatch\n']);
figure('Name', sprintf('Sx / Sy convention check, CoreView %d frame %d', P.run_num, P.conv_frame), 'Position', [100 100 1100 420]);
subplot(1, 2, 1);  semilogx(P.conv_q * mmx_, res_(:, 1), 'o-', P.conv_q * mmx_, res_(:, 2), 's-', P.conv_q * mmx_, res_(:, 4), 'd--');
grid on;  xlabel('block size (mm)');  ylabel('rms(curl) / rms(div)');  legend({'as processed', 'Sy flipped', 'best b'}, 'Location', 'best');
title('lower = more integrable');
subplot(1, 2, 2);  semilogx(P.conv_q * mmx_, res_(:, 3), 'o-');  yline(1, 'k--');  yline(-1, 'k:');
grid on;  xlabel('block size (mm)');  ylabel('best Sy factor b');  title('+1 = consistent, -1 = Sy sign flipped');
clear R_ ref_b_ setup_b_ sx_ sy_ mmx_ mmy_ sc_ res_ j_ q_ nr_ nc_ rs_ v_ bx_ by_ axc_ axr_ ayc_ ayr_ A_ B_ m_ s0_ s1_ b_ sb_

%% 6. Slope cache + video (long: run it when nothing else uses MATLAB)
P.video.frames  = 1000:3246;
P.video.out_fps = 50;
if isempty(P.video.frames)                         % every frame on disk (bad ones become blank video frames)
    S_ = cisg_scan_run(P.run_num);  P.video.frames = S_.frames(:).';  clear S_
end
fprintf('\n--- 6. slope cache + video, %d frames (%d to %d) ---\n', numel(P.video.frames), P.video.frames(1), P.video.frames(end));
cfg_override = P.video;  cfg_override.coreview_num = P.run_num;
cisg_slope_video;                                 % clears the workspace except cfg_override
diary off;

%% 6b. Re-render the video from the slope cache (no solving: run section 0 first)
% Same settings as section 6, so the same cache folder.  Colour limits: [] = 99.5th percentile of |S| and
% |eta| over 30 frames spread through the record (the old default, 3 x rms of the first frame, saturated
% once the waves grew); or set them, e.g. 0.4 and 8 (mm).  A missing cache frame is solved serially.
P.video.frames    = 1000:3246;
P.video.out_fps   = 50;
P.video.lim_s     = [];
P.video.lim_eta   = [];
P.video.lim_pct   = 99.5;
P.video.use_parfor = false;                       % no pool: everything is cached
P.video.render_only = true;                       % never solve: an uncached frame becomes a blank video frame
delete(gcp('nocreate'));                          % the section-6 workers can hold tens of GB: free the RAM for the render
P.video.out_path  = fullfile(P.video.out_root, sprintf('CISG_slopes_CoreView_%d_frames_%d-%d_lim.mp4', ...
                             P.run_num, P.video.frames(1), P.video.frames(end)));
fprintf('\n--- 6b. re-render from cache, %d frames ---\n', numel(P.video.frames));
cfg_override = P.video;  cfg_override.coreview_num = P.run_num;
cisg_slope_video;
diary off;

%% 7. OPTIONAL, skipped by default: outlier pass over the cache (masks stored in the cache files, slopes untouched)
% Flagged only 0.3% of run 101 frame 2000 (section 5c), so it is not part of the normal run.
% Section 6 cleared the workspace: run section 0 first.  This finds the newest cache folder of the run.
L_ = dir(fullfile(P.video.out_root, sprintf('CISG_slopes_CoreView_%d_bin%d*_yup*', P.run_num, P.video.bin)));   % Sy y-up caches only
L_ = L_([L_.isdir]);
if isempty(L_), error('cisg_pipeline:cache', 'No slope cache for CoreView %d in %s', P.run_num, P.video.out_root); end
[~, k_] = max([L_.datenum]);  P.cache_dir = fullfile(L_(k_).folder, L_(k_).name);
fprintf('\n--- 7. outlier pass on %s ---\n', P.cache_dir);
diary(P.log);  diary on;
P.outlier_result = cisg_outlier_pass(P.cache_dir, [], P.outlier);
diary off;  clear L_ k_
