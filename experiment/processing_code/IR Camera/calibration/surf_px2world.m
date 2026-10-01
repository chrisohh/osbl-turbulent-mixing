function [xy, pu] = surf_px2world(cal, pd)
%SURF_PX2WORLD  Raw (distorted) IR pixels -> water-surface tank coords (mm).
%
%   [xy, pu] = surf_px2world(cal, [u v])   % xy = [x_along_mm, y_cross_mm]
%
% Inverse of surf_world2px.m. The radial term has no closed-form inverse,
% so the undistorted point is found by fixed-point iteration (converges in
% a few iterations for the moderate barrel distortion of the 13 mm lens).
% pu is the undistorted pixel position, returned for diagnostics.

d  = (pd - cal.center) / cal.focal_px;
du = d;
for it = 1:50
    r2 = sum(du.^2, 2);
    g  = ones(size(r2));
    for i = 1:numel(cal.k), g = g + cal.k(i) * r2.^i; end
    du = d ./ g;
end
pu = cal.center + cal.focal_px * du;

q  = [pu, ones(size(pu,1),1)] / cal.T;           % inverse homography
xy = q(:,1:2) ./ q(:,3);
end
