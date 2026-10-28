%% Process several Water_SURF runs in ONE MATLAB session, one after another
% Runs surf_theta_video.m for each run number in run_queue (settings per
% run: util/lif_run_config.m), then compare_depth_runs.m.  Start it from
% your own MATLAB desktop and leave that session open until it finishes.
%
% surf_theta_video clears everything except run_num and run_queue, so the
% queue is a variable that is consumed from the front.
if ~exist('run_queue', 'var'), run_queue = [92 95 89 102 120 88 116 117]; end   % 96 and 118 are done
% (102 needs its wave-gauge file set in util/lif_run_config.m first; without
%  it, it fails with a message and the queue moves on)
diary(fullfile(fileparts(mfilename('fullpath')), sprintf('process_runs_%s.log', datestr(now, 'yyyymmdd_HHMM'))));
while ~isempty(run_queue)
    run_num = run_queue(1);  run_queue(1) = [];
    fprintf('\n===== run %d  (%s; still queued: %s) =====\n', run_num, datestr(now, 'HH:MM:SS'), mat2str(run_queue));
    try
        surf_theta_video;
    catch ME
        fprintf(2, 'run %d FAILED: %s\n', run_num, ME.message);
    end
    close all;
end
fprintf('\n===== comparison (%s) =====\n', datestr(now, 'HH:MM:SS'));
try
    compare_depth_runs;                 % all processed runs; unprocessed ones are skipped
catch ME
    fprintf(2, 'comparison FAILED: %s\n', ME.message);
end
diary off;
