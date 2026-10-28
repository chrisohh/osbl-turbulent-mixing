function R = lif_trans_config(run_num)
%LIF_TRANS_CONFIG Transverse (laser-sheet) LIF settings of a run, from the registry.
%
%   R = LIF_TRANS_CONFIG(run_num)
%
% Reads processing_code/common/run_registry.csv (edit that, not this) for
% the transverse view: trans_cam (the Water_PIV / bottom camera of the
% run), trans_plate + trans_plate_cam (plate run and its camera folder),
% trans_ref_frame (frame with the undisturbed dye layer on the surface: the
% z = 0 datum and the sheet profile S(y), see lif_trans_conc), and the same
% wind-clock fields as the side view (lif_wind_clock).
%
% Laser-sheet runs: 104, 105 (unstrat, plate 103) and 119, 121 (strat,
% plate = the laser-sheet calibration image in run 121, camera C0170400).
% In other runs the Water_PIV camera is the LED-lit angled view.

    here = fileparts(mfilename('fullpath'));
    addpath(fullfile(here, '..', '..', 'common'));
    G = run_registry(run_num);
    if isempty(G.trans_cam)
        error('lif_trans_config:none', 'Run %d has no transverse laser-sheet view in the registry.', run_num);
    end
    core = 'D:\HLAB_2026\Core';
    pf = arrayfun(@(k) fullfile(core, sprintf('CoreView_%d', G.trans_plate), G.trans_plate_cam, ...
        sprintf('CoreView_%d_%s_%02d.raw', G.trans_plate, G.trans_plate_cam, k)), 1:5, 'UniformOutput', false);
    lbl = sprintf('%d', run_num);
    if strcmp(G.condition, 'strat'), lbl = [lbl ' strat']; end
    R = struct('label', lbl, 'run', run_num, 'condition', G.condition, ...
        'cam', G.trans_cam, 'plate_num', G.trans_plate, 'plate_cam', G.trans_plate_cam, ...
        'plate_file', {pf}, 'ref_frame', G.trans_ref_frame, ...
        'hw_file', G.wind_path, 'wind_kind', G.wind_kind, 'cam_delay_s', G.cam_delay_s, ...
        'standin_file', G.standin_path, 'density_file', G.ct_pre_path, ...
        'wave_file', G.wave_path, 'wave_standin_file', G.wave_standin_path, ...
        'stratified', strcmp(G.condition, 'strat'), 'h_cm', G.h_cm, 'dye', G.dye, ...
        'best', G.best, 'onset_s', G.onset_s, 'check', G.check);
    if isempty(R.density_file), R.density_file = []; end
end
