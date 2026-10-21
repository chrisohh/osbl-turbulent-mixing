function [mmpp, c0, info] = lif_trans_scale_cached(fplate, sc, cache_file)
%LIF_TRANS_SCALE_CACHED lif_trans_scale with a small cache file.
%
%   [mmpp, c0, info] = LIF_TRANS_SCALE_CACHED(fplate, sc, cache_file)
%
% Same arrangement as lif_scale_cached: the plate images live on the lab
% computer, the time series on another one (no network between them).
% Where the plate image exists the scale is computed and written to
% cache_file; elsewhere cache_file is loaded.  Commit cache_file (a few
% hundred bytes) with the code.
%
% sc may instead carry a field mmpp (and optionally c0): the plate is then
% skipped entirely and those values are returned.  Use that while no plate
% has been imaged in the laser sheet.

    if isfield(sc, 'mmpp') && ~isempty(sc.mmpp)
        mmpp = sc.mmpp;
        if isfield(sc, 'c0') && ~isempty(sc.c0), c0 = sc.c0; else, c0 = [2048 1536]; end
        info = struct('source', 'manual');
        fprintf('Scale: MANUAL %.5f mm/px, origin [%.0f %.0f] (no plate used)\n', mmpp, c0);
        return
    end

    if isfile(fplate)
        [mmpp, c0, info] = lif_trans_scale(fplate, sc);
        if ~isfolder(fileparts(cache_file)), mkdir(fileparts(cache_file)); end
        save(cache_file, 'mmpp', 'c0', 'info', 'sc', 'fplate');
        fprintf('  scale cached to %s\n', cache_file);
    elseif isfile(cache_file)
        S = load(cache_file);
        mmpp = S.mmpp;  c0 = S.c0;  info = S.info;
        fprintf('Plate image not on this computer; scale from %s: %.5f mm/px (square %g mm)\n', ...
                cache_file, mmpp, S.sc.square_mm);
        if ~isequal(S.sc, sc)
            warning('Scale settings (sc) differ from the cached ones -- rerun where the plate image is.');
        end
    else
        error('Neither the plate image %s nor the cache %s exists.', fplate, cache_file);
    end
end
