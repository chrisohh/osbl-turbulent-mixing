function R = lif_run_config(run_num)
%LIF_RUN_CONFIG Side-view (Water_SURF) LIF settings of a run, from the registry.
%
%   R = LIF_RUN_CONFIG(run_num)
%
% Reads processing_code/common/run_registry.csv (edit that file, not this
% one) and maps it to the fields the LIF scripts use.
%
% Reference plane: the side view integrates absorption over the whole tank
% width, so depths are referenced to the tank CENTRE plane (the width
% average).  side_offset_mm = how far the plate stood toward the camera
% from the centre.  With the 60 mm lens the camera is ~1.35 m from the
% plate, so dye at the front glass appears ~16% deeper and at the back glass
% ~12% shallower than at the centre; z99 (tail-driven) leans toward the
% front-glass value, z_rms barely.
%
% Fields: label, run, condition, cam, plate_num, plate_file (plate frames,
%   first usable is taken), square_mm, offset_mm, lens_f_mm, hw_file,
%   wind_kind, cam_delay_s, standin_file, wave_file, wave_standin_file,
%   bg_frames, density_file,
%   stratified, line_bin_cm, h_cm, dye, best, onset_s, check

    here = fileparts(mfilename('fullpath'));
    addpath(fullfile(here, '..', '..', 'common'));
    G = run_registry(run_num);
    if isempty(G.side_cam)
        error('lif_run_config:noside', 'Run %d has no side-view camera in the registry.', run_num);
    end

    lbl = sprintf('%d', run_num);
    if ~isempty(G.dye) && ~any(strcmp(G.dye, {'normal', '?'})), lbl = sprintf('%s (%s dye)', lbl, G.dye); end
    if strcmp(G.condition, 'strat'), lbl = [lbl ' strat']; end

    R = struct('label', lbl, 'run', run_num, 'condition', G.condition, ...
        'cam', G.side_cam, 'plate_num', G.side_plate, 'plate_file', {G.side_plate_files}, ...
        'square_mm', 11, 'offset_mm', G.side_offset_mm, 'lens_f_mm', 60, ...
        'hw_file', G.wind_path, 'wind_kind', G.wind_kind, 'cam_delay_s', G.cam_delay_s, ...
        'standin_file', G.standin_path, 'wave_file', G.wave_path, ...
        'wave_standin_file', G.wave_standin_path, 'bg_frames', 1:10:100, ...
        'density_file', G.ct_pre_path, 'stratified', strcmp(G.condition, 'strat'), ...
        'line_bin_cm', G.line_bin_cm, 'h_cm', G.h_cm, 'dye', G.dye, 'best', G.best, ...
        'onset_s', G.onset_s, 'check', G.check);
    if isnan(R.line_bin_cm), R.line_bin_cm = 0; end
    if isempty(R.density_file), R.density_file = []; end
end
