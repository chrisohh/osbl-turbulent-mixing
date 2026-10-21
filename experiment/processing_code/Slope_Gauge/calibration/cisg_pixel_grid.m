function [X_cm, Y_cm] = cisg_pixel_grid(cal, roi)
%CISG_PIXEL_GRID Tank coordinates (cm) of every CISG pixel in a ROI.
%
%   [X_cm, Y_cm] = cisg_pixel_grid(cal, roi)
%
% cal - cisg_plate_calibrate output (surf_px2world compatible)
% roi - struct with row_min row_max col_min col_max (as cisg_main /
%       pick_roi_interactive), or [row_min row_max col_min col_max]
%
% X_cm = along-wind, Y_cm = cross-wind (IR tank frame, origin at the plate
% centre).  Both are single, size [nrows ncols] of the ROI.

    if isstruct(roi)
        roi = [roi.row_min roi.row_max roi.col_min roi.col_max];
    end
    rr = roi(1):roi(2);
    cc = roi(3):roi(4);
    X_cm = zeros(numel(rr), numel(cc), 'single');
    Y_cm = zeros(numel(rr), numel(cc), 'single');

    % surf_px2world iterates a fixed-point inverse: do it a block of columns at a time
    nblk = max(1, floor(5e5 / numel(rr)));
    for j0 = 1:nblk:numel(cc)
        j = j0:min(j0 + nblk - 1, numel(cc));
        [C, R] = meshgrid(cc(j), rr);
        xy = surf_px2world(cal, [C(:) R(:)]);              % mm
        X_cm(:, j) = single(reshape(xy(:,1), numel(rr), numel(j)) / 10);
        Y_cm(:, j) = single(reshape(xy(:,2), numel(rr), numel(j)) / 10);
    end
end
