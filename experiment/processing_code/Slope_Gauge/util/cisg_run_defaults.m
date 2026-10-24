function rd = cisg_run_defaults(coreview_num)
%CISG_RUN_DEFAULTS Per-run defaults shared by cisg_setup_run.m and cisg_main.m.
%
%   rd = cisg_run_defaults(coreview_num)
%
% rd.plate       CoreView number of the plate frame used for the spatial calibration
%                (calibration/cisg_cal_CoreView_<plate>.mat): runs 95, 96 -> 97; 101 -> 106; others -> 97
% rd.dye_model   'joint' for the rhodamine 6G runs (95, 96: per-pixel dye unmixing, cisg_dye_unmix.m),
%                'none' otherwise (101 is dye-free)
% rd.ir_snap     IR snapshot of the same plate placement: 95, 96 -> Snap-000096; 101 -> Snap-000101
% rd.ablated_parity  parity of (m + n) of the bright (ablated) squares in that IR snapshot's tank frame (measured from the
%                IR images: odd for Snap-000096, even for Snap-000101) - input cfg.ablated_parity of cisg_plate_calibrate
% rd.ir_cal_path IR Camera/calibration/surface_cal_Snap-<ir_snap>.mat (may not exist yet: check isfile);
%                use with cisg_ir_to_cisg(load(rd.ir_cal_path).cal, cal_cisg, uv_ir)

    plate = 97;  dye = 'none';  snap = 96;  parity = 1;
    switch coreview_num
        case {95, 96},  plate = 97;  dye = 'joint';  snap = 96;
        case 101,       plate = 106;  snap = 101;  parity = 0;
    end
    pc = fileparts(fileparts(fileparts(mfilename('fullpath'))));        % processing_code
    rd = struct('plate', plate, 'dye_model', dye, 'ir_snap', snap, 'ablated_parity', parity, ...
                'ir_cal_path', fullfile(pc, 'IR Camera', 'calibration', sprintf('surface_cal_Snap-%06d.mat', snap)));
end
