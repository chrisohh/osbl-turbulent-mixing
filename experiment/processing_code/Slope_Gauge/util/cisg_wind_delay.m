function [cam_delay_s, src] = cisg_wind_delay(coreview_num)
%CISG_WIND_DELAY Camera start after wind (fan) start for a run, from common/run_registry.csv.
%
%   [cam_delay_s, src] = cisg_wind_delay(coreview_num)
%
% Frame n is at cam_delay_s + (n - 1)/fs after wind start.  Same source as the LIF videos
% (lif_hw_cached): the run's wind record (hot-wire / reference probe / wave gauge: recorded camera
% trigger, else software start times), cached in Slope_Gauge/calibration/hwsync_<file>_avg1.mat so
% it also works where the record is not on disk; else the registry cam_delay_s; else NaN.
% src: text saying where the value came from.

    sg = fileparts(fileparts(mfilename('fullpath')));
    pc = fileparts(sg);
    addpath(fullfile(pc, 'common'), fullfile(pc, 'LIF', 'util'));      % run_registry, lif_hw_cached
    cam_delay_s = NaN;  src = 'unknown';
    try
        R = run_registry(coreview_num);
    catch ME
        warning('cisg_wind_delay:registry', 'Run %d not in run_registry (%s): no wind clock', coreview_num, ME.message);
        return;
    end
    if ~isempty(R.wind_path)
        [~, nm] = fileparts(R.wind_path);
        cache = fullfile(sg, 'calibration', sprintf('hwsync_%s_avg1.mat', nm));
        if isfile(R.wind_path) || isfile(cache)
            cam_delay_s = lif_hw_cached(R.wind_path, 1, cache);
            src = nm;
            return;
        end
    end
    if isfield(R, 'cam_delay_s') && ~isempty(R.cam_delay_s) && isfinite(R.cam_delay_s)
        cam_delay_s = R.cam_delay_s;  src = 'registry cam_delay_s (assumed)';
    end
end
