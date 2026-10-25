function dark = cisg_ref_dark(ref_image, setup)
%CISG_REF_DARK Pixels that are unusable in the reference itself: outside the target columns or in shadow / dark.
%
%   dark = cisg_ref_dark(ref_image, setup)     % logical, same size as ref_image
%
% Same criteria as cisg_auto_roi (mean RGB / 255 < setup.dark_min after smoothing, or outside setup.target_cols),
% evaluated at 1/4 resolution and expanded.  Used by roi_mode 'full' to leave only data-limited areas empty.
% dark_min defaults to 0.04 (black: walls, shadow, off-target), not the 0.15 of the auto ROI: dim but coloured
% pixels are kept and judged by the precision mask (sigma_S from camera noise, cisg_colour_displacement), since a
% brightness threshold cannot tell a dim usable gradient from an unusable one.

    bm = 0.04;  if isfield(setup, 'dark_min') && ~isempty(setup.dark_min), bm = setup.dark_min; end
    [~, ai] = cisg_auto_roi(ref_image, setup, struct('min_frac', 0, 'bright_min', bm));
    dark = imresize(ai.dark, [size(ref_image, 1) size(ref_image, 2)], 'nearest') > 0;
end
