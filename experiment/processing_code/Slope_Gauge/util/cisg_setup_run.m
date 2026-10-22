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
%   ref_frame     1          still-water reference = this frame of the run
%   roi_mode      'default'  'default' (central 60%), 'auto' (from the reference, see cisg_auto_roi; auto_margin px inset),
%                            'pick' (interactive) or a roi struct
%   roi_file      ''         'pick': where the picked ROI is saved / reused
%                            ('' = calibration/cisg_roi_CoreView_<N>.mat)
%   roi_repick    false      'pick': true = pick again even if roi_file exists
%   cal_path      ''         '' = by run: 95, 96 -> calibration/cisg_cal_CoreView_97.mat, 101 -> ..._106.mat (default 97)
%   raw_root      'D:\HLAB_2026\Core'
%   fs            50         frame rate (Hz)
%   setup         struct     fields overriding the defaults below (colour_map, colour_sigma, ...)
%
% R fields: setup, cal, depth, ref_image, obs_folder, frame_file (@(n) path), dt, img_size

    here = fileparts(mfilename('fullpath'));
    sg = fileparts(here);
    addpath(fullfile(sg, 'calibration'), fullfile(fileparts(sg), 'IR Camera', 'calibration'));

    def = struct('coreview_num', 95, 'ref_frame', 1, 'roi_mode', 'default', ...
                 'roi_file', '', 'roi_repick', false, 'auto_margin', 0, 'cal_path', '', ...
                 'raw_root', 'D:\HLAB_2026\Core', 'fs', 50, 'setup', struct());
    if nargin < 1, cfg = struct(); end
    for f = fieldnames(def).'
        if ~isfield(cfg, f{1}), cfg.(f{1}) = def.(f{1}); end
    end

    cam = 'Flare 12M125 CCL C1112A00';
    R.obs_folder = fullfile(cfg.raw_root, sprintf('CoreView_%d', cfg.coreview_num), cam);
    R.frame_file = @(n) fullfile(R.obs_folder, sprintf('CoreView_%d_%s_%04d.raw', cfg.coreview_num, cam, n));
    R.dt = 1 / cfg.fs;

    setup = struct('water_depth', 48, 'use_apparent_depth', true, 'colour_map', 'local', ...
                   'colour_norm', 'none', 'colour_sigma', 12, 'max_disp_px', 600, 'drift_corr', 'none', ...
                   'n_water', 1.333, 'n_air', 1.000, 'frame_rate', cfg.fs, 'dt', R.dt);
    for f = fieldnames(cfg.setup).'
        setup.(f{1}) = cfg.setup.(f{1});
    end

    cal_path = cfg.cal_path;
    if isempty(cal_path)                                  % plate frame used for each run (same camera setup)
        lut = containers.Map({95, 96, 101}, {97, 97, 106});
        plate = 97;  if isKey(lut, cfg.coreview_num), plate = lut(cfg.coreview_num); end
        cal_path = fullfile(sg, 'calibration', sprintf('cisg_cal_CoreView_%d.mat', plate));
    end
    S = load(cal_path);
    cal = S.cal;
    [setup.cm_per_pixel_x, setup.cm_per_pixel_y] = cisg_local_scale(cal);

    ref_image = cisg_load_coreview(R.frame_file(cfg.ref_frame));
    [h, w, ~] = size(ref_image);

    depth = cisg_apparent_depth(ref_image, cal, struct('n_w', setup.n_water));
    setup.L_eff_cm      = depth.L_eff_mm / 10;
    setup.pattern_scale = depth.r;
    setup.target_cols   = depth.edge_cols;

    roi_default = struct('row_min', round(0.2*h), 'row_max', round(0.8*h), ...
                         'col_min', round(0.2*w), 'col_max', round(0.8*w));
    if isstruct(cfg.roi_mode)
        setup.roi = cfg.roi_mode;
    elseif strcmp(cfg.roi_mode, 'auto')
        % ROI that avoids what is already unusable in the reference (wall shadow, flat colour map)
        setup.roi = cisg_auto_roi(ref_image, setup, struct('margin_px', cfg.auto_margin));
        fprintf('auto ROI: rows %d-%d, cols %d-%d\n', setup.roi.row_min, setup.roi.row_max, setup.roi.col_min, setup.roi.col_max);
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

    R.setup = setup;  R.cal = cal;  R.depth = depth;  R.ref_image = ref_image;  R.img_size = [h w];
end
