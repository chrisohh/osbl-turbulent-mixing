function cisg_fix_path(sg)
%CISG_FIX_PATH Make the Slope_Gauge folders win on the path and drop stale worktree copies.
%
%   cisg_fix_path(sg)      sg = the Slope_Gauge folder
%
% A genpath on processing_code also picks up copies under '.claude\worktrees', whose older
% cisg_*.m files can shadow the real ones (e.g. a cisg_load_coreview without expected_size).
% Removes every path entry containing '\.claude\', puts Slope_Gauge, util, calibration and
% IR Camera/calibration first, and clears loaded functions (a long-lived parallel pool
% otherwise keeps old versions).  Also run on the pool workers by cisg_slope_video.m.

    p = strsplit(path, pathsep);
    stale = p(contains(p, [filesep '.claude' filesep]) | endsWith(p, [filesep '.claude']));
    if ~isempty(stale), rmpath(stale{:}); end
    addpath(fullfile(fileparts(sg), 'IR Camera', 'calibration'), fullfile(sg, 'calibration'), fullfile(sg, 'util'), sg);
    clear functions
end
