function [mmpp, c0] = lif_scale_cached(fplate, sc, cache_file)
%LIF_SCALE_CACHED lif_plate_scale with a small cache file.
%
%   [mmpp, c0] = LIF_SCALE_CACHED(fplate, sc, cache_file)
%
% The plate images live on the lab computer, the time series on another
% one (no network between them).  Where the plate image exists, the scale
% is computed and written to cache_file; elsewhere cache_file is loaded.
% Copy / commit cache_file (a few bytes) with the code.

    if isfile(fplate)
        [mmpp, c0] = lif_plate_scale(fplate, sc);
        if ~isfolder(fileparts(cache_file)), mkdir(fileparts(cache_file)); end
        save(cache_file, 'mmpp', 'c0', 'sc', 'fplate');
        fprintf('  scale cached to %s\n', cache_file);
    elseif isfile(cache_file)
        S = load(cache_file);
        mmpp = S.mmpp;  c0 = S.c0;
        fprintf('Plate image not on this computer; scale from %s: %.4f mm/px (square %g mm)\n', ...
                cache_file, mmpp, S.sc.square_mm);
        if ~isequal(S.sc, sc)
            warning('Scale settings (sc) differ from the cached ones -- rerun where the plate image is.');
        end
    else
        error('Neither the plate image %s nor the cache %s exists.', fplate, cache_file);
    end
end
