function out = run_registry(arg, varargin)
%RUN_REGISTRY Per-run metadata for all instruments, from run_registry.csv.
%
%   R = RUN_REGISTRY(run)                one run as a struct (paths resolved)
%   T = RUN_REGISTRY()                   the whole table
%   runs = RUN_REGISTRY('list', 'condition', 'unstrat')   run numbers by filter
%   runs = RUN_REGISTRY('list', 'condition', 'strat', 'best', true)
%   runs = RUN_REGISTRY('list', 'has', 'side_cam')        runs with a side view
%
% run_registry.csv (next to this file) is the single source of truth: edit
% it in Excel.  One row per CoreView run; empty cell = not available.
% Columns:
%   run, condition (unstrat|strat), best (1 = representative run), h_cm,
%   dye (none|low|normal|...),
%   side_cam / side_plate / side_plate_cam / side_offset_mm
%       side-view (Water_SURF) camera folder of the run; plate run, the
%       camera folder in the plate run (folder names differ between runs),
%       plate distance toward the camera from the tank centre (mm)
%   trans_cam / trans_plate      transverse (bottom, Water_PIV) view
%   nadir_cam / nadir_plate      nadir colour camera (CISG)
%   ir_rec / ir_snap             IR recording Rec-N and its calibration Snap-N
%   wind_file, wind_kind         wind clock: hotwire | ref_probe | hotwire+wave
%                                | wave | none; cam_delay_s = camera start
%                                after wind start when there is no record
%   wind_standin                 another run's wind record shown (labelled)
%                                when this run has none
%   wave_file, wave_ref_file     wave gauge during the run / dye-free repeat
%   ct_pre                       CT profile taken BEFORE the run (strat)
%   line_bin_cm                  z99(x) line window for the side view
%   onset_s                      picked onset time (s since wind start)
%   check, notes                 open questions / remarks

    persistent T
    here = fileparts(mfilename('fullpath'));
    f = fullfile(here, 'run_registry.csv');
    d = dir(f);
    if isempty(T) || ~isequal(T.Properties.UserData, d.datenum)   % re-read when edited
        opts = detectImportOptions(f, 'TextType', 'char', 'Delimiter', ',');
        T = readtable(f, opts);
        T.Properties.UserData = d.datenum;
    end

    if nargin == 0, out = T; return; end

    if ischar(arg) && strcmp(arg, 'list')
        keep = true(height(T), 1);
        for k = 1:2:numel(varargin)
            key = varargin{k};  val = varargin{k+1};
            if strcmp(key, 'has')
                c = T.(val);
                if iscell(c), keep = keep & ~cellfun(@isempty, c); else, keep = keep & ~isnan(c); end
            elseif ischar(val)
                keep = keep & strcmp(T.(key), val);
            else
                keep = keep & (T.(key) == val);
            end
        end
        out = T.run(keep).';
        return
    end

    i = find(T.run == arg, 1);
    if isempty(i), error('run_registry:run', 'Run %d is not in %s.', arg, f); end
    R = table2struct(T(i, :));

    % resolved paths
    core = 'D:\HLAB_2026\Core';
    R.wind_path    = wind_path(R.wind_file);
    R.standin_path = wind_path(R.wind_standin);
    R.wave_path    = wind_path(R.wave_file);
    R.wave_ref_path = wind_path(R.wave_ref_file);
    R.ct_pre_path  = '';
    if ~isempty(R.ct_pre)
        R.ct_pre_path = fullfile('D:\HLAB_2026\CT', [R.ct_pre '.mat']);
    end
    R.side_plate_files = {};
    if ~isempty(R.side_plate_cam) && ~isnan(R.side_plate)
        R.side_plate_files = arrayfun(@(k) fullfile(core, sprintf('CoreView_%d', R.side_plate), ...
            R.side_plate_cam, sprintf('CoreView_%d_%s_%02d.raw', R.side_plate, R.side_plate_cam, k)), ...
            1:5, 'UniformOutput', false);
    end
    out = R;
end

function p = wind_path(name)
    p = '';
    if isempty(name), return; end
    if startsWith(name, 'wave')
        p = fullfile('D:\HLAB_2026\WaveGauge', [name '.mat']);
    else
        p = fullfile('D:\HLAB_2026\hotwire', [name '.mat']);
    end
end
