function [ref_b, setup_b] = cisg_bin_setup(ref_image, setup, bin)
%CISG_BIN_SETUP Reference image and setup scaled for block-averaged (binned) processing.
%
%   [ref_b, setup_b] = cisg_bin_setup(ref_image, setup, bin)
%
% bin = 1 returns the inputs.  For bin > 1 the image is block-averaged, pixel
% sizes (cm/px) grow by bin, the ROI / target columns map to binned pixels and
% the colour smoothing and the displacement clamp (px) shrink by bin, so the
% physical scales stay the same.

    ref_b = ref_image;  setup_b = setup;
    if isfield(setup, 'noise') && ~isempty(setup.noise) && setup.noise.bin ~= bin
        warning('cisg_bin_setup:noise', 'Noise estimated at bin %d but used at bin %d: precision mask is approximate', setup.noise.bin, bin);
    end
    if bin == 1, return; end
    ref_b = uint8(imresize(ref_image, 1 / bin, 'box'));
    setup_b.cm_per_pixel_x = setup.cm_per_pixel_x * bin;
    setup_b.cm_per_pixel_y = setup.cm_per_pixel_y * bin;
    setup_b.colour_sigma   = setup.colour_sigma / bin;
    setup_b.max_disp_px    = setup.max_disp_px / bin;
    setup_b.target_cols    = (setup.target_cols - 0.5) / bin + 0.5;
    r = setup.roi;
    setup_b.roi = struct('row_min', ceil(r.row_min / bin), 'row_max', floor(r.row_max / bin), ...
                         'col_min', ceil(r.col_min / bin), 'col_max', floor(r.col_max / bin));
    % pattern_scale is the pattern/surface size ratio: unchanged
    if isfield(setup, 'dark_mask') && ~isempty(setup.dark_mask)
        setup_b.dark_mask = imresize(setup.dark_mask, [size(ref_b, 1) size(ref_b, 2)], 'nearest') > 0;
    end
end
