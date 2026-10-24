function dark = cisg_ref_dark(ref_image, setup)
%CISG_REF_DARK Pixels that are unusable in the reference itself: outside the target columns or in shadow / dark.
%
%   dark = cisg_ref_dark(ref_image, setup)     % logical, same size as ref_image
%
% Same brightness criteria as cisg_auto_roi (mean RGB / 255 < 0.15 after smoothing, or outside setup.target_cols),
% evaluated at 1/4 resolution and expanded.  Used by roi_mode 'full' to leave only data-limited areas empty.

    [~, ai] = cisg_auto_roi(ref_image, setup, struct('min_frac', 0));
    dark = imresize(ai.dark, [size(ref_image, 1) size(ref_image, 2)], 'nearest') > 0;
end
