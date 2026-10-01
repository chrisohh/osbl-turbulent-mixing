function pd = surf_world2px(cal, xy)
%SURF_WORLD2PX  Water-surface tank coords (mm) -> raw (distorted) IR pixels.
%
%   pd = surf_world2px(cal, [x_along_mm, y_cross_mm])   % N x 2 -> N x 2 [u v]
%
% cal comes from calibrate_surface_single.m. Model:
%   p_u = homography(cal.T, xy)                           (undistorted px)
%   p_d = c + (p_u - c) .* (1 + k1 r^2 + k2 r^4 + k3 r^6),  r = |p_u - c| / f

q  = [xy, ones(size(xy,1),1)] * cal.T;          % MATLAB row-vector convention
pu = q(:,1:2) ./ q(:,3);

d  = (pu - cal.center) / cal.focal_px;
r2 = sum(d.^2, 2);
pd = cal.center + cal.focal_px * d .* radial_factor(cal.k, r2);
end

function g = radial_factor(k, r2)
% 1 + k1 r^2 + k2 r^4 + k3 r^6 + ... for any number of coefficients
g = ones(size(r2));
for i = 1:numel(k)
    g = g + k(i) * r2.^i;
end
end
