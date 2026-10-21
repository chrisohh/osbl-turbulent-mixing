%% CoreView 105 transverse LIF frame: raw | C | profile
% Transverse (cross-tank, y-z) laser sheet through rhodamine, Water_PIV
% camera.  This is TRUE LIF: the sheet excites the dye and the dye emits,
% so the image is black except where there is rhodamine and
%       C  propto  counts
% directly -- no back-light, no dye-free background frame, and none of the
% Beer-Lambert machinery the side views (surf_frame_calibrated.m) need.
% The helpers are therefore the lif_trans_* family, not lif_dye_front /
% lif_z99:
%   util/lif_trans_scale   mm/px from the yawed Water_PIV plate
%   util/lif_trans_surface surface line (brightest, y-averaged; the camera
%                          is rolled ~6 deg so it images as a slant)
%   util/lif_trans_conc    noise floor, sheet normalisation, profile, z99
%   util/lif_trans_figure + lif_trans_update   the figure, shared with
%                          trans_conc_video.m so the still and the video match
%
% This is the paper's own plane (eq. 3.14 averages over y), so the
% frame-wide z99 here is the directly comparable number -- unlike the
% longitudinal view, where z99 has to be taken locally in x.
%
% Time: frame n is trigger n of the camera counter,
%     t = (camStartElapsed - fanStartElapsed) + (n - 1) / fs
% seconds since wind start, from the run_experiment.m hot-wire file.  No
% hot-wire file has been recorded for CoreView 105 yet, so hw_file is
% empty and the title falls back to time since the first frame; set
% hw_file and show_hw when one exists.
clear; close all;
addpath(fullfile(fileparts(mfilename('fullpath')), 'util'));

%% ---------------- User settings ----------------
raw_root  = 'D:\HLAB_2026\Core';
cam       = 'Water_PIV';          % the transverse sheet was shot with the "PIV" camera
run_num   = 105;    frame     = 1843;
ref_frame = 1;                    % undisturbed dye sheet: sets the z = 0 datum

% ---- Scale ------------------------------------------------------------
% The only plate image for this camera (CoreView 91) is a vertical strip
% yawed ~59 deg to the sensor, so its squares image 173 px tall and 88 px
% wide.  lif_plate_scale's uniform fit would average the two; lif_trans_scale
% takes the VERTICAL step, which the yaw does not foreshorten -- right for
% depths, approximate for the cross-tank axis.  Set sc.mmpp to override
% with a hand value, or point plate_num at a plate imaged IN the sheet.
plate_num = 91;     plate_frm = 1;
sc.square_mm = 11;        % small printed plate (same one as Water_SURF)
sc.offset_mm = 0;         % plate-to-sheet offset toward the camera (CHECK)
sc.n_medium  = 1.33;
sc.lens_f_mm = 35;        % <-- CHECK
sc.mmpp      = [];        % [] = from the plate; a number overrides it

hw_file = '';             % '' = no hot-wire file for this run yet
fs      = 50;             % camera trigger rate (Hz) -- CHECK for run 105

show_hw  = ~isempty(hw_file);
hw_avg_s = 0.02;
hw_xlim  = [];

downsample = 2;           % every Nth pixel
% A bright reflection runs down the right edge near column 4000; drop it.
% The ROWS are deliberately NOT cropped: lif_trans_conc needs clear water
% well below the deepest dye to set its noise floor, and the dye reaches
% ~2 cm down in a frame that is ~19 cm tall.  Zoom with fo.zlim_cm instead.
col_crop = [1 3900];
row_crop = [1 3072];

% ---- Figure (see lif_trans_figure) ------------------------------------
fo.cmap        = 'hot';
fo.raw_log     = true;          % surface line ~150 counts, deep plumes ~1
fo.clims_log   = [0.5 300];
fo.clims       = [0 150];
fo.C_log       = true;          % same reason: log shows the tail that sets z99
fo.C_clims_log = [3e-3 1];
fo.C_clims     = [0 1];
fo.C_label     = '$C/C_{ref}$';   % 1 = the frame-1 surface dye layer
fo.profile_xlim = [];           % [] = autoscale per frame
fo.profile_zlim = [-4 0.5];     % cm about the local surface
fo.zlim_cm      = [-4 1.5];     % image panels: zoom on the dye (data stay whole)
fo.ylim_cm      = [];           % [] = the whole lit width

% ---- Surface (see lif_trans_surface) ----------------------------------
su.smooth_mm      = 1.0;
su.surf_smooth_cm = 3;     % along-y average: keeps the line, washes plumes out
su.search_cm      = [];    % [] = whole frame; set once the datum is known
su.min_counts     = 15;    % columns dimmer than this are not lit by the sheet
su.fit_order      = 1;     % camera roll -> a straight slant (order 1 matches
                           % frame 1843's own line to +-0.02 cm; order 3 only
                           % fits the meniscus)
su.track          = false; % false: hold the frame-1 line.  true lets the whole
                           % line ride with the surface, but once the dye
                           % detaches the picks sit on plumes and drag it down
su.offset_cm      = 0;     % shift the datum down if the emitting layer sits low

% ---- Concentration + z99 (see lif_trans_conc) -------------------------
co.frac        = 0.99;
co.smooth_mm   = su.smooth_mm;
co.surf_gap_cm = 0;        % raise to drop a saturated surface line
co.top_cm      = 0.3;
co.deep_cm     = 6;        % clear water starts this far below the surface
co.deep_min_rows = 200;    % skip columns with less clear water than that
co.min_peak    = 15;      % counts in the REFERENCE frame: below this the sheet
                           % barely lights the column and dividing by it amplifies noise
co.int_depth_cm = 5;       % the 99% integral stops here below the surface;
                           % [] = the whole frame.  CHECK the sensitivity.
co.normalize   = 'ref';    % divide each column by the sheet profile S(y), taken
                           % from the reference frame; 'column' | 'peak' | 'none'
co.bin_cm      = 1;
co.min_bin_frac = 0.8;     % a z99(y) bin needs this fraction of its columns lit
co.pen_frac    = 0.02;     % dye-penetration contour, fraction of the column peak
co.pen_min_counts = 5;     % ...and an absolute floor (~10 sigma of the read noise),
                           % so one faint speck deep down cannot set the contour
co.pen_smooth_cm  = 1.0;   % median filter along y
co.ref         = [];       % filled with the reference frame below: it defines
                           % which columns the sheet lights
co.dark        = [];       % counts (scalar) or a lens-capped frame; [] = none

%% ---------------- Scale ----------------
cal_dir = fullfile(fileparts(mfilename('fullpath')), 'calibration');
fplate  = fullfile(raw_root, sprintf('CoreView_%d', plate_num), cam, ...
                   sprintf('CoreView_%d_%s_%02d.raw', plate_num, cam, plate_frm));
[mmpp, c0] = lif_trans_scale_cached(fplate, sc, ...
                 fullfile(cal_dir, sprintf('scale_%s_Core%d_trans.mat', cam, plate_num)));

%% ---------------- Load frame + reference ----------------
fname = @(n) fullfile(raw_root, sprintf('CoreView_%d', run_num), cam, ...
                      sprintf('CoreView_%d_%s_%04d.raw', run_num, cam, n));
crop  = @(A) A(row_crop(1):downsample:row_crop(2), col_crop(1):downsample:col_crop(2));

I    = crop(lif_load_raw(fname(frame)));
Iref = crop(lif_load_raw(fname(ref_frame)));
if ~isempty(co.dark) && ischar(co.dark), co.dark = crop(lif_load_raw(co.dark)); end
[ny, nx] = size(I);
dx_mm = mmpp * downsample;
y_cm =  ((col_crop(1) + (0:nx-1)*downsample) - c0(1)) * mmpp / 10;   % full-res pixel index
z_cm = -((row_crop(1) + (0:ny-1)*downsample) - c0(2)) * mmpp / 10;   % row 1 (top) -> largest z

fprintf('Frame %d: %dx%d px after crop/downsample, %.5f mm/px (%.2f x %.2f cm)\n', ...
        frame, ny, nx, dx_mm, nx*dx_mm/10, ny*dx_mm/10);

%% ---------------- Surface datum ----------------
% The reference frame is the undisturbed dye sheet lying on the surface,
% so its line IS the water level: fit it there and carry the fit over, the
% way lif_surface_datum uses the first frame in the side views.
[z_ref, Sref] = lif_trans_surface(Iref, z_cm, dx_mm, su);
z0   = mean(z_ref, 'omitnan');
z_cm = z_cm - z0;
z_ref = z_ref - z0;
fprintf('Datum: frame %d water level at %.2f cm above the board centre -> z = 0\n', ref_frame, z0);
fprintf('  reference line: %d columns, RMS %.1f px, camera roll %.2f deg\n', ...
        Sref.n_used, Sref.rms_px, Sref.roll_deg);

su.prior = Sref;                       % keep the shape, let it ride with the surface
[z_surf, Ssurf] = lif_trans_surface(I, z_cm, dx_mm, su);
fprintf('  frame %d surface: %.2f..%.2f cm (mean %.2f cm, %.2f cm off the datum)\n', ...
        frame, min(z_surf), max(z_surf), mean(z_surf), mean(z_surf - z_ref));

%% ---------------- Time ----------------
if ~isempty(hw_file)
    [~, hw_name] = fileparts(hw_file);
    [camDelay, hw] = lif_hw_cached(hw_file, hw_avg_s, fullfile(cal_dir, ['hwsync_' hw_name '.mat']));
    t_frame = camDelay + (frame - 1) / fs;
    t_label = sprintf('$t = %.2f$ s since wind start', t_frame);
    fprintf('Frame %d: t = %.3f s since wind start (cameras started at %.3f s)\n', ...
            frame, t_frame, camDelay);
else
    hw = [];
    t_frame = (frame - 1) / fs;
    t_label = sprintf('$t = %.2f$ s after frame 1', t_frame);
    warning(['No hot-wire file for CoreView %d: the title shows time since frame 1, ' ...
             'not since wind start.  Set hw_file once run_experiment.m has one.'], run_num);
end

%% ---------------- Concentration, profile, z99 ----------------
co.ref = Iref;
R = lif_trans_conc(I, z_cm, z_surf, dx_mm, co);
fprintf('Sheet: %d of %d columns lit (peak %.0f..%.0f counts); noise floor %.2f..%.2f counts\n', ...
        nnz(~isnan(R.peak_col)), nx, min(R.peak_col), max(R.peak_col), ...
        min(R.floor_col), max(R.floor_col));
fprintf('Dye penetration (%.0f%% of the column peak): %.2f cm mean, %.2f cm deepest\n', ...
        100*co.pen_frac, mean(z_surf - R.z_pen, 'omitnan'), max(z_surf - R.z_pen));
fprintf('Sheet profile S(y) from frame %d: %s counts at the [5 25 50 75 95] percentiles\n', ...
        ref_frame, mat2str(round(prctile(R.ref_col, [5 25 50 75 95]))));
fprintf('z99 = %.2f cm below the surface (all y); z99(y) median %.2f, 5-95%% %.2f..%.2f cm (%d bins of %g cm)\n', ...
        R.z99, median(R.z99_y, 'omitnan'), prctile(R.z99_y, 5), prctile(R.z99_y, 95), ...
        nnz(~isnan(R.z99_y)), co.bin_cm);

% z99 is a 1% tail, so it is only as good as the floor subtraction: show how
% much it moves with the integration depth before quoting it.
for d = [3 4 5 6 8]
    co_d = co;  co_d.int_depth_cm = d;
    Rd = lif_trans_conc(I, z_cm, z_surf, dx_mm, co_d);
    fprintf('    int_depth %2g cm -> z99 = %6.2f cm\n', d, Rd.z99);
end

%% ---------------- Plot ----------------
if show_hw
    fo.hw = hw;
    fo.hw.xlim = hw_xlim;
end
h = lif_trans_figure(y_cm, z_cm, fo);
set(h.fig, 'Name', sprintf('CoreView_%d %s frame %d', run_num, cam, frame));
lif_trans_update(h, I, R, y_cm, z_surf, ...
    sprintf('CoreView %d, transverse LIF (%s), frame %d, %s', ...
            run_num, strrep(cam, '_', '\_'), frame, t_label), t_frame);
