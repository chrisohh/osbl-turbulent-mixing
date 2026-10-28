%% Process several transverse (laser-sheet) LIF runs in ONE MATLAB session
% Runs trans_conc_video.m for each run in run_queue (settings per run:
% processing_code/common/run_registry.csv), then fig_trans_vs_side.m.
% Laser-sheet runs: 104, 105 (unstrat), 119, 121 (strat).  Start it from
% your own MATLAB desktop and leave that session open until it finishes.
%
% trans_conc_video clears everything except run_num and run_queue, so the
% queue is consumed from the front.
if ~exist('run_queue', 'var'), run_queue = [105 104 121 119]; end
diary(fullfile(fileparts(mfilename('fullpath')), sprintf('process_trans_runs_%s.log', datestr(now, 'yyyymmdd_HHMM'))));
while ~isempty(run_queue)
    run_num = run_queue(1);  run_queue(1) = [];
    fprintf('\n===== transverse run %d  (%s; still queued: %s) =====\n', run_num, datestr(now, 'HH:MM:SS'), mat2str(run_queue));
    try
        trans_conc_video;
    catch ME
        fprintf(2, 'run %d FAILED: %s\n', run_num, ME.message);
    end
    close all;
end
fprintf('\n===== transverse vs side (%s) =====\n', datestr(now, 'HH:MM:SS'));
try
    fig_trans_vs_side;
catch ME
    fprintf(2, 'comparison FAILED: %s\n', ME.message);
end
diary off;
