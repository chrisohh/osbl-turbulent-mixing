function [uv_cisg, xy_mm] = cisg_ir_to_cisg(cal_ir, cal_cisg, uv_ir)
%CISG_IR_TO_CISG Map IR pixels to CISG pixels through the shared surface plane.
%
%   [uv_cisg, xy_mm] = cisg_ir_to_cisg(cal_ir, cal_cisg, uv_ir)
%
% uv_ir   - N x 2 IR pixels [u v] (e.g. an outline of an IR structure)
% xy_mm   - N x 2 CAMERA-aligned tank coordinates [x_downwind y_up] in mm
%           (origin at the plate centre): the IR plate-frame coordinates
%           (surf_px2world of the IR cal) rotated by cal_cisg.plate_to_tank
% uv_cisg - N x 2 CISG pixels [u v] (surf_world2px of the CISG cal)
%
% Both cals describe the water-surface plane (plate at the surface), so the
% map is exact for anything ON the surface; anything above or below it is
% displaced by parallax.  Points outside the trusted footprint
% (cal.bbox_mm) are extrapolated by the homography.  The inverse map is
% cisg_cisg_to_ir.m.

    xy_plate = surf_px2world(cal_ir, uv_ir);                       % IR plate frame
    xy_mm    = xy_plate * plate_to_tank(cal_cisg)';                % tank = R * plate (row vectors)
    uv_cisg  = surf_world2px(cal_cisg, xy_mm);
end

function R = plate_to_tank(cal)
    if isfield(cal, 'plate_to_tank'), R = cal.plate_to_tank; else, R = eye(2); end
end
