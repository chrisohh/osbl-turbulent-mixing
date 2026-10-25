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
% rd.cal_path    Slope_Gauge/calibration/cisg_cal_CoreView_<plate>.mat
% rd.plate_frame CISG plate frame (frame 1 of CoreView_<plate>)
% rd.ir_snap_path IR snapshot file (D:\HLAB_2026\IR_camera\Snap-<ir_snap>.ats)
% rd.ir_flat_path board-free IR recording for the narcissus flat of calibrate_surface_single
%                (101 -> Rec-000100; '' = not recorded here, the existing IR cal is kept)

    plate = 97;  dye = 'none';  snap = 96;  parity = 1;  flat = '';
    switch coreview_num
        case {95, 96},  plate = 97;  dye = 'joint';  snap = 96;
        case 101,       plate = 106;  snap = 101;  parity = 0;  flat = 'Rec-000100.ats';
    end
    pc = fileparts(fileparts(fileparts(mfilename('fullpath'))));        % processing_code
    ir_root = 'D:\HLAB_2026\IR_camera';
    if ~isempty(flat), flat = fullfile(ir_root, flat); end
    rd = struct('plate', plate, 'dye_model', dye, 'ir_snap', snap, 'ablated_parity', parity, ...
                'ir_cal_path', fullfile(pc, 'IR Camera', 'calibration', sprintf('surface_cal_Snap-%06d.mat', snap)), ...
                'cal_path', fullfile(pc, 'Slope_Gauge', 'calibration', sprintf('cisg_cal_CoreView_%d.mat', plate)), ...
                'plate_frame', first_frame(plate), ...
                'ir_snap_path', fullfile(ir_root, sprintf('Snap-%06d.ats', snap)), 'ir_flat_path', flat);
end

function f = first_frame(n)
% frame 1 of CoreView_<n>; the zero padding follows the recording length (_01, _001, _0001)
    d = fullfile('D:\HLAB_2026\Core', sprintf('CoreView_%d', n), 'Flare 12M125 CCL C1112A00');
    stem = sprintf('CoreView_%d_Flare 12M125 CCL C1112A00_', n);
    L = dir(fullfile(d, [stem '*1.raw']));
    k = find(~cellfun(@isempty, regexp({L.name}, ['^' regexptranslate('escape', stem) '0*1\.raw$'], 'once')), 1);
    if isempty(k), f = fullfile(d, [stem '001.raw']); else, f = fullfile(d, L(k).name); end
end
