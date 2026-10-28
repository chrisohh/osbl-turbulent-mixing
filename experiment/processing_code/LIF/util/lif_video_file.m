function vid_file = lif_video_file(out_file)
%LIF_VIDEO_FILE The video file to write: out_file, or a time-stamped copy if it is locked.
%
%   vid_file = LIF_VIDEO_FILE(out_file)
%
% An existing video that is open somewhere (a player, Explorer's preview
% pane) cannot be overwritten: VideoWriter fails with "Permission Denied"
% after the set-up, and process_runs loses the run.  Write next to it
% instead, as <name>_<yyyymmdd_HHMM>.mp4, and say so.  Only the video moves;
% the _z99.mat time series keeps its usual name.
    vid_file = out_file;
    if ~isfile(out_file), return; end
    fid = fopen(out_file, 'a');
    if fid >= 0, fclose(fid); return; end
    [p, b, e] = fileparts(out_file);
    vid_file = fullfile(p, sprintf('%s_%s%s', b, datestr(now, 'yyyymmdd_HHMM'), e));
    warning('lif_video_file:locked', '%s is open elsewhere (player / preview?) -- writing %s instead.', out_file, vid_file);
end
