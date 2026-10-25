function R = cisg_setup_run(cfg)
%CISG_SETUP_RUN Calibration, still-water reference, apparent depth and ROI for one CISG run.
%
%   R = cisg_setup_run(cfg)
%
% The same setup as cisg_main.m (kept in step with it by hand), without the
% display and without writing anything.  Used by cisg_slope_video.m.
%
% cfg fields (defaults):
%   coreview_num  95
%   ref_frames    2:4:58     still-water reference = per-pixel MEDIAN of these frames of the run (flat water at the
%                            start).  Not frame 1: on run 101 it is ~0.8 counts brighter in R, G, B than frames 2-200
%                            (first-frame exposure), which reads as a uniform slope of ~0.01 (cisg_pipeline 5d); a
%                            median also drops anything in one frame only (a speck, a drip).  A scalar = that frame
%                            alone (ref_frame, the old field, is still accepted).
%   roi_mode      'default'  'default' (central 60%), 'full' (whole image, dark / off-target pixels of the reference masked), 'auto' (from the reference, see cisg_auto_roi; auto_margin px inset),
%                            'pick' (interactive) or a roi struct
%   roi_file      ''         'pick': where the picked ROI is saved / reused
%                            ('' = calibration/cisg_roi_CoreView_<N>.mat)
%   roi_repick    false      'pick': true = pick again even if roi_file exists
%   cal_path      ''         '' = by run (cisg_run_defaults: 95, 96 -> cisg_cal_CoreView_97.mat, 101 -> ..._106.mat)
%   dye_model     ''         '' = by run (cisg_run_defaults: 'joint' for the R6G runs 95, 96, else 'none'); or 'none' | 'joint'
%   raw_root      'D:\HLAB_2026\Core'
%   fs            50         frame rate (Hz)
%   bin           2          bin factor the noise estimate (precision mask) is made for (cisg_slope_video passes cfg.bin)
%   setup         struct     fields overriding the defaults below (colour_map, colour_sigma, ...)
%
% R fields: setup, cal, depth, ref_image, obs_folder, frame_file (@(n) path), dt, img_size,
%   cam_delay_s (camera start after wind start, cisg_wind_delay; NaN = unknown), t_wind (@(n) s after wind start)

    here = fileparts(mfilename('fullpath'));
    sg = fileparts(here);
    addpath(fullfile(sg, 'calibration'), fullfile(fileparts(sg), 'IR Camera', 'calibration'));

    def = struct('coreview_num', 95, 'ref_frames', 2:4:58, 'roi_mode', 'default', ...
                 'roi_file', '', 'roi_repick', false, 'auto_margin', 0, 'bin', 2, 'cal_path', '', 'dye_model', '', ...
                 'raw_root', 'D:\HLAB_2026\Core', 'fs', 50, 'setup', struct());
    if nargin < 1, cfg = struct(); end
    if isfield(cfg, 'ref_frame') && ~isfield(cfg, 'ref_frames'), cfg.ref_frames = cfg.ref_frame; end
    for f = fieldnames(def).'
        if ~isfield(cfg, f{1}), cfg.(f{1}) = def.(f{1}); end
    end

    cam = 'Flare 12M125 CCL C1112A00';
    R.obs_folder = fullfile(cfg.raw_root, sprintf('CoreView_%d', cfg.coreview_num), cam);
    R.frame_file = @(n) fullfile(R.obs_folder, sprintf('CoreView_%d_%s_%04d.raw', cfg.coreview_num, cam, n));
    R.dt = 1 / cfg.fs;

    setup = struct('water_depth', 48, 'use_apparent_depth', true, 'colour_map', 'local', ...
                   'colour_norm', 'none', 'colour_sigma', 12, 'max_disp_px', 600, 'drift_corr', 'none', 'colour_channels', 'RG', ...
                   'n_water', 1.333, 'n_air', 1.000, 'frame_rate', cfg.fs, 'dt', R.dt);
    for f = fieldnames(cfg.setup).'
        setup.(f{1}) = cfg.setup.(f{1});
    end

    cal_path = cfg.cal_path;
    rd = cisg_run_defaults(cfg.coreview_num);
    if isempty(cal_path)                                  % plate frame used for each run (same camera setup)
        cal_path = fullfile(sg, 'calibration', sprintf('cisg_cal_CoreView_%d.mat', rd.plate));
    end
    setup.dye_model = cfg.dye_model;  if isempty(setup.dye_model), setup.dye_model = rd.dye_model; end
    fprintf('CoreView %d: calibration %s, dye_model = %s\n', cfg.coreview_num, cal_path(find(cal_path == filesep, 1, 'last') + 1:end), setup.dye_model);
    S = load(cal_path);
    cal = S.cal;
    [setup.cm_per_pixel_x, setup.cm_per_pixel_y] = cisg_local_scale(cal);

    % still-water reference: per-pixel median of cfg.ref_frames (unreadable frames skipped)
    sz0 = [3072 4100];                                                      % some runs (101) are all headerless
    stack = zeros([sz0 3 numel(cfg.ref_frames)], 'uint8');  nref = 0;
    for n = cfg.ref_frames
        try
            nref = nref + 1;  stack(:, :, :, nref) = cisg_load_coreview(R.frame_file(n), sz0);
        catch ME
            if ~strcmp(ME.identifier, 'cisg:badFrame'), rethrow(ME); end
            nref = nref - 1;  fprintf('reference: frame %d unreadable, skipped\n', n);
        end
    end
    if nref == 0, error('cisg_setup_run:ref', 'No readable reference frame in %s', mat2str(cfg.ref_frames)); end
    if nref == 1
        ref_image = stack(:, :, :, 1);
    else
        ref_image = zeros([sz0 3], 'uint8');
        for c = 1:3, ref_image(:, :, c) = median(stack(:, :, c, 1:nref), 4); end
    end
    clear stack
    fprintf('reference: median of %d frame(s) %s\n', nref, mat2str(cfg.ref_frames));
    R.ref_frames = cfg.ref_frames;  R.n_ref = nref;
    [h, w, ~] = size(ref_image);

    depth = cisg_apparent_depth(ref_image, cal, struct('n_w', setup.n_water));
    setup.L_eff_cm      = depth.L_eff_mm / 10;
    setup.pattern_scale = depth.r;
    setup.target_cols   = depth.edge_cols;
    % channel centre line (y = 0 in the outputs) from the wall edges of this reference; plate centre if not detectable
    cc = cisg_channel_center(ref_image, cal, setup);
    setup.y0_cm = cc.y0_cm;
    if ~isfinite(setup.y0_cm)
        warning('cisg_setup_run:y0', 'Wall edges not both detected in the reference (top %d, bottom %d): y = 0 stays at the plate-lattice centre.', cc.found_top, cc.found_bottom);
        setup.y0_cm = 0;
    else
        fprintf('channel centre: wall edges rows %.0f / %.0f, apparent midpoint %.2f cm, y0 = %.2f cm (plate-lattice y = 0)\n', cc.top_row, cc.bottom_row, cc.mid_apparent_cm, setup.y0_cm);
    end

    roi_default = struct('row_min', round(0.2*h), 'row_max', round(0.8*h), ...
                         'col_min', round(0.2*w), 'col_max', round(0.8*w));
    if isstruct(cfg.roi_mode)
        setup.roi = cfg.roi_mode;
    elseif strcmp(cfg.roi_mode, 'auto')
        % ROI that avoids what is already unusable in the reference (wall shadow, flat colour map)
        setup.roi = cisg_auto_roi(ref_image, setup, struct('margin_px', cfg.auto_margin));
        fprintf('auto ROI: rows %d-%d, cols %d-%d\n', setup.roi.row_min, setup.roi.row_max, setup.roi.col_min, setup.roi.col_max);
    elseif strcmp(cfg.roi_mode, 'full')
        % the whole image; only data-limited areas stay empty (dark / outside the target in the reference + solver mask)
        setup.roi = struct('row_min', 1, 'row_max', h, 'col_min', 1, 'col_max', w);
        setup.dark_mask = cisg_ref_dark(ref_image, setup);
    elseif strcmp(cfg.roi_mode, 'default')
        setup.roi = roi_default;
    else
        % 'pick': pick once, save, reuse on later runs (roi_repick = true to
        % pick again; the saved / default box is drawn as the starting guess)
        roi_file = cfg.roi_file;
        if isempty(roi_file)
            roi_file = fullfile(sg, 'calibration', sprintf('cisg_roi_CoreView_%d.mat', cfg.coreview_num));
        end
        if isfile(roi_file) && ~cfg.roi_repick
            S = load(roi_file, 'roi');  setup.roi = S.roi;
            fprintf('ROI from %s (roi_repick = true to pick again)\n', roi_file);
        else
            prior = roi_default;
            if isfile(roi_file), S = load(roi_file, 'roi'); prior = S.roi; end
            setup.roi = pick_roi_interactive(ref_image, prior, setup.cm_per_pixel_x, setup.cm_per_pixel_y, ...
                sprintf('CoreView %d: pick the slope ROI (two clicks)', cfg.coreview_num));
            roi = setup.roi;  save(roi_file, 'roi');
            fprintf('ROI saved to %s\n', roi_file);
        end
    end
    setup.roi_is_default = isequal(setup.roi, roi_default);
    roi_centre = [mean([setup.roi.col_min setup.roi.col_max]) mean([setup.roi.row_min setup.roi.row_max])];
    [setup.cm_per_pixel_x, setup.cm_per_pixel_y] = cisg_local_scale(cal, roi_centre);

    % camera noise of the smoothed colours (two early frames) -> physical precision mask (setup.max_slope_sigma, slope units)
    % frames 2, 4 (frame 1 is off in exposure).  The estimate is of a difference of two single frames, sqrt(2) sigma;
    % against a median of n_ref frames the reference noise is ~sqrt(pi/2 / n_ref) sigma, so rescale to that
    setup.noise = cisg_noise_estimate(R.frame_file, cfg.bin, setup.colour_sigma / cfg.bin, [2 4], [h w]);
    fr = sqrt((1 + (pi / 2) / nref) / 2);
    setup.noise.a = setup.noise.a * fr;  setup.noise.b = setup.noise.b * fr;  setup.noise.ref_factor = fr;
    setup.max_slope_sigma = 0.02;
    R.setup = setup;  R.cal = cal;  R.depth = depth;  R.ref_image = ref_image;  R.img_size = [h w];  R.channel = cc;
    [R.cam_delay_s, src] = cisg_wind_delay(cfg.coreview_num);           % frame n at cam_delay_s + (n-1) dt after wind start
    fprintf('wind clock: camera start %.3f s after wind start (%s)\n', R.cam_delay_s, src);
    R.t_wind = @(n) R.cam_delay_s + (n - 1) * R.dt;
end
