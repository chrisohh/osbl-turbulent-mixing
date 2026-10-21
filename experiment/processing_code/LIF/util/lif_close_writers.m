function lif_close_writers(writers)
%LIF_CLOSE_WRITERS  close() every VideoWriter in a cell array, never throwing.
%
%   lif_close_writers({vid.v})
%
% For onCleanup: if the frame loop errors (or you ctrl-C it) every writer
% still gets closed, so the .mp4 files are playable instead of truncated.
% Each is wrapped on its own, so one bad writer cannot stop the rest from
% being closed.

    for k = 1:numel(writers)
        try
            close(writers{k});
        catch err
            warning('lif_close_writers:failed', 'Could not close writer %d: %s', k, err.message);
        end
    end
end
