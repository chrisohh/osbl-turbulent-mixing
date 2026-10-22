function R = lif_run_config(run_num)
%LIF_RUN_CONFIG Per-run settings for the Water_SURF (side view) LIF runs.
%
%   R = LIF_RUN_CONFIG(run_num)
%
% One entry per CoreView run, so each run is processed with its own plate
% calibration and hot-wire file and runs can be compared (compare_depth_runs.m).
%
% Reference plane: the side view integrates absorption over the whole tank
% width, so depths are referenced to the tank CENTRE plane (the width
% average).  offset_mm = how far the plate stood toward the camera from the
% centre plane.  With the 60 mm lens the camera is ~1.35 m from the plate,
% so dye at the front glass appears ~16% deeper and at the back glass ~12%
% shallower than at the centre; z99 (tail-driven) is biased toward the
% front-glass value, z_rms barely.
%
% Fields:
%   label, cam            - run label, camera folder of the run
%   plate_num, plate_file - calibration run and full path of its plate image
%   square_mm             - plate square size (small plate: 11 mm)
%   offset_mm             - plate toward the camera from the tank centre (mm)
%   lens_f_mm             - side-view lens (Nikon 60 mm)
%   hw_file               - hot-wire file of the run (time since wind start)
%   bg_frames             - dye-free frames averaged as the background
%                           (CHECK per run: all must be before the dye arrives)
%   density_file          - CT profile (depth cm positive down from the surface,
%                           conductivity voltage), [] = unstratified
%   stratified            - logical

    root = 'D:\HLAB_2026\Core';
    hw   = 'D:\HLAB_2026\hotwire';
    % plate frames _01.._05: lif_plate_scale uses the first one that works
    plate = @(n, cam) arrayfun(@(k) fullfile(root, sprintf('CoreView_%d', n), cam, ...
                               sprintf('CoreView_%d_%s_%02d.raw', n, cam, k)), 1:5, 'UniformOutput', false);

    R = struct('label', '', 'cam', 'Water_Surf', 'plate_num', [], 'plate_file', '', ...
               'square_mm', 11, 'offset_mm', 25, 'lens_f_mm', 60, 'hw_file', '', ...
               'bg_frames', 1:10:100, 'density_file', [], 'stratified', false);

    switch run_num
        case 92
            R.label = '92';      R.plate_num = 93;
            R.plate_file = plate(93, 'Water_Surf');
            R.hw_file = fullfile(hw, 'hotwire_20260923_081422.mat');
        case 95
            R.label = '95';      R.plate_num = 94;
            R.plate_file = plate(94, 'Water_Surf');
            R.hw_file = fullfile(hw, 'hotwire_20260923_113939.mat');
        case 96
            R.label = '96 (lower dye concentration)';  R.plate_num = 94;
            R.plate_file = plate(94, 'Water_Surf');
            R.hw_file = fullfile(hw, 'hotwire_20260923_123939.mat');
        case 118
            % Plate run 122 has no Water_Surf folder and its folder names do
            % not match the cameras: the side (Water_Surf) camera is the one
            % saved as 'CCL C1112A00' (25 MB mono frames); 'MCL C0170400' is
            % the angled Water_PIV view and 'MCL C0100800' the nadir colour
            % camera (37.8 MB frames).
            R.label = '118 (stratified)';  R.plate_num = 122;
            R.plate_file = plate(122, 'Flare 12M125 CCL C1112A00');
            R.offset_mm = 0;                       % plate at the tank centre
            R.hw_file = fullfile(hw, 'hotwire_20261006_004837.mat');
            R.density_file = 'D:\HLAB_2026\CT\densityProfileData_20261008_0033.mat';
            R.stratified = true;
        otherwise
            error('lif_run_config:run', 'No settings for run %d -- add it to lif_run_config.m.', run_num);
    end
end
