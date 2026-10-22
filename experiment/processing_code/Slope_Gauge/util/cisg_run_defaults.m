function rd = cisg_run_defaults(coreview_num)
%CISG_RUN_DEFAULTS Per-run defaults shared by cisg_setup_run.m and cisg_main.m.
%
%   rd = cisg_run_defaults(coreview_num)
%
% rd.plate       CoreView number of the plate frame used for the spatial calibration
%                (calibration/cisg_cal_CoreView_<plate>.mat): runs 95, 96 -> 97; 101 -> 106; others -> 97
% rd.dye_model   'joint' for the rhodamine 6G runs (95, 96: per-pixel dye unmixing, cisg_dye_unmix.m),
%                'none' otherwise (101 is dye-free)

    plate = 97;  dye = 'none';
    switch coreview_num
        case {95, 96},  plate = 97;  dye = 'joint';
        case 101,       plate = 106;
    end
    rd = struct('plate', plate, 'dye_model', dye);
end
