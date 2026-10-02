function [camDelay, hw] = lif_hw_cached(hw_file, avg_s, cache_file)
%LIF_HW_CACHED Camera start time + hot-wire series, with a small cache file.
%
%   [camDelay, hw] = LIF_HW_CACHED(hw_file, avg_s, cache_file)
%
% camDelay - camera start after wind start (s): frame n is at
%            camDelay + (n-1)/fs
% hw       - lif_hw_series(H, avg_s): block-averaged U / U_ref vs time
%
% The hot-wire file lives on the lab computer, the time series on another
% one (no network between them).  Where hw_file exists both are computed
% and written to cache_file (~100 kB); elsewhere cache_file is loaded.

    if isfile(hw_file)
        H = load(hw_file);
        if isfield(H, 'camStartElapsed') && ~isnan(H.camStartElapsed)
            camDelay = H.camStartElapsed - H.fanStartElapsed;
        else
            camDelay = H.runConfig.DELAY_BEFORE_TRIG;
            warning('No camStartElapsed in %s -- using DELAY_BEFORE_TRIG = %g s.', hw_file, camDelay);
        end
        hw = lif_hw_series(H, avg_s);
        if ~isfolder(fileparts(cache_file)), mkdir(fileparts(cache_file)); end
        save(cache_file, 'camDelay', 'hw', 'avg_s', 'hw_file');
        fprintf('Hot-wire sync cached to %s\n', cache_file);
    elseif isfile(cache_file)
        S = load(cache_file);
        camDelay = S.camDelay;  hw = S.hw;
        fprintf('Hot-wire file not on this computer; sync from %s (camDelay %.3f s)\n', cache_file, camDelay);
        if S.avg_s ~= avg_s
            warning('Cached hot-wire series uses avg_s = %g s, not %g s.', S.avg_s, avg_s);
        end
    else
        error('Neither the hot-wire file %s nor the cache %s exists.', hw_file, cache_file);
    end
end
