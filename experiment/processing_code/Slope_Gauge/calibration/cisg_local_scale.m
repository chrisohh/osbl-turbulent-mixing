function [cm_per_px_x, cm_per_px_y] = cisg_local_scale(cal, uv)
%CISG_LOCAL_SCALE Local pixel size (cm/px) at a pixel, along the image axes.
%
%   [cm_per_px_x, cm_per_px_y] = cisg_local_scale(cal, [u v])
%
% x = image horizontal (columns), y = image vertical (rows): the length on
% the surface of a one-pixel step in each direction.  Used for
% setup.cm_per_pixel_x/y in cisg_main.m (cisg_calculate_slopes, ROI picker).
% uv defaults to the image centre.

    if nargin < 2 || isempty(uv), uv = (fliplr(cal.image_size) + 1) / 2; end
    xy = surf_px2world(cal, [uv - [0.5 0]; uv + [0.5 0]; uv - [0 0.5]; uv + [0 0.5]]);   % mm
    cm_per_px_x = norm(xy(2,:) - xy(1,:)) / 10;
    cm_per_px_y = norm(xy(4,:) - xy(3,:)) / 10;
end
