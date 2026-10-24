function p = cisg_pattern_v1_params()
%CISG_PATTERN_V1_PARAMS Parameters of the v1 colour target (the printed 20" x 20" RGB gradient sheet).
%
%   p = cisg_pattern_v1_params()
%
% Record only: nothing here generates or changes anything.  v1 is the target used for runs up to at least
% CoreView 101 (the user confirms the same print was used throughout).  Source of truth: Slope_Gauge/CISG_sheet.m
% (git history: first version 2025-11-19 as an A3 sheet, 29.7 x 42.0 cm at 300 dpi; changed on 2025-12-11, commit
% 62f5d51, to 20" x 20" at 600 dpi; unchanged since).  Slope_Gauge/cisg_design_pattern.m (commit c3edd25,
% 2026-02-02) is a DESIGN helper (sizes the pattern from camera height / FOV, writes CISG_pattern.png) that makes the
% same gradient with the same channel assignment but a computed size and dpi; it is NOT what was printed here
% (no run-specific branches in either file).
%
% Pattern (per pixel of the sheet image, X = column index, Y = row index, both 1-based):
%   R = uint8(255 * (X - 1) / (width - 1))        linear in X, 0 -> 255, "increases left to right"
%   G = uint8(255 * (Y - 1) / (height - 1))       linear in Y, 0 -> 255 with increasing ROW index
%                                                  (down the image as displayed by imshow; the script comment says
%                                                  "bottom to top", its ylabel says down: the file is rows-down)
%   B = uint8(128)                                constant
% No gamma / tone curve is applied in the file (values are written as is).  What the camera sees differs
% (printing, lighting, illumination fall-off): in the run 95 reference R and G are non-monotonic near the edges
% and B is NOT constant (it varies by about +-30 counts with position, correlated with row).  Observed orientation in
% the CISG image: R decreases with image column and G decreases with image row (corr -0.6 / -0.7), i.e. the printed
% sheet is seen mirrored / rotated relative to the file; the inversion handles this through the measured gradients.
%
% Uncertain: (1) whether the printed sheet matches CISG_sheet.m exactly (printer colour management / ICC profile,
% paper, gamma) - the file itself has no gamma; (2) the physical placement of the sheet (which sheet edge is image
% left, flips) is not recorded in the code; (3) cisg_main.m comments mention camera_height 41.7 cm and a 50.8 cm
% pattern, i.e. the same 20" sheet.

    p = struct();
    p.version        = 'v1';
    p.source_script  = 'Slope_Gauge/CISG_sheet.m';
    p.git_history    = {'2025-11-19 a8fd823 first version: 29.7 x 42.0 cm (A3), 300 dpi'; ...
                        '2025-12-11 62f5d51 20 in x 20 in, 600 dpi (printed v1)'};
    p.print_size_in  = [20 20];
    p.print_size_cm  = 20 * 2.54 * [1 1];                    % 50.8 cm
    p.print_size_mm  = 508 * [1 1];
    p.dpi            = 600;
    p.size_px        = round(p.print_size_cm / 2.54 * p.dpi);   % [12000 12000] (height, width)
    p.px_per_mm      = p.dpi / 25.4;                            % 23.62 px/mm
    p.channels       = struct( ...
        'R', struct('axis', 'X (columns)', 'direction', 'left to right in the file', 'range', [0 255], 'shape', 'linear', ...
                    'formula', 'uint8(255*(X-1)/(width-1))'), ...
        'G', struct('axis', 'Y (rows)', 'direction', 'increasing row index (top to bottom as displayed)', 'range', [0 255], ...
                    'shape', 'linear', 'formula', 'uint8(255*(Y-1)/(height-1))'), ...
        'B', struct('axis', 'none', 'direction', 'constant', 'range', [128 128], 'shape', 'constant', 'formula', 'uint8(128)'));
    p.blue_level     = 128;
    p.gamma          = 1;                                       % none applied
    p.bit_depth      = 8;
    p.markers        = 'none (corner white squares are commented out in the script)';
    p.output_files   = {'CISG_gradient.png', 'CISG_gradient.tif (uncompressed)'};
    p.colour_gradient_per_mm = struct('R_counts_per_mm', 255 / 508, 'G_counts_per_mm', 255 / 508);   % ideal, at the pattern plane

    % inputs of the design helper cisg_design_pattern.m (NOT used for the print): its internal constants
    p.design_helper = struct('file', 'Slope_Gauge/cisg_design_pattern.m', 'sensor_mm', [8.8 6.6], 'max_slope', 0.3, ...
        'size_rule_cm', 'ceil((FOV_pattern + 2*max_displacement)/5)*5', 'dpi_rule', 'max(300, ceil(2.54/(spatial_res/4)))', ...
        'saves', {{'CISG_pattern.png', 'CISG_design.mat'}}, 'git', '2026-02-02 c3edd25');

    % camera / geometry of the runs that used this target (from the code comments and the cals; for reference)
    p.geometry = struct('camera', 'Flare 12M125 (CMV12000), 5.5 um pixels, Canon 35 mm lens', ...
        'camera_height_cm_comment', 41.7, 'camera_height_cm_from_cal', 48.3, 'water_depth_cm', [48 50.1], ...
        'glass_mm_assumed', 19, 'air_gap_mm_assumed', 76.2, 'target_apparent_width_mm', [257 260 244], ...
        'note', 'heights from cisg_apparent_depth: runs 95, 96, 101; see cisg_setup_run');
end
