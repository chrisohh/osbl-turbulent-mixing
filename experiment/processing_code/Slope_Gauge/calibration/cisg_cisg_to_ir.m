function [uv_ir, xy_mm] = cisg_cisg_to_ir(cal_ir, cal_cisg, uv_cisg)
%CISG_CISG_TO_IR Map CISG pixels to IR pixels through the shared surface plane.
%
%   [uv_ir, xy_mm] = cisg_cisg_to_ir(cal_ir, cal_cisg, uv_cisg)
%
% Inverse of cisg_ir_to_cisg.m (same arguments, same plane assumption).
% xy_mm is the camera-aligned tank position; the IR cal is applied in the
% plate frame (cal_cisg.plate_to_tank undone).

    xy_mm = surf_px2world(cal_cisg, uv_cisg);                      % tank
    if isfield(cal_cisg, 'plate_to_tank'), R = cal_cisg.plate_to_tank; else, R = eye(2); end
    uv_ir = surf_world2px(cal_ir, xy_mm * R);                      % plate = R' * tank (row vectors)
end
