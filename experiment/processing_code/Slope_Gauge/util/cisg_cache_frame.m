function [Sx, Sy] = cisg_cache_frame(n, ref_b, setup_b, frame_file, cache_dir, bin, dt, full_size)
%CISG_CACHE_FRAME Slopes of one frame inside the ROI, cached as frame_NNNN.mat (same format as cisg_main).
%
%   [Sx, Sy] = cisg_cache_frame(n, ref_b, setup_b, frame_file, cache_dir, bin, dt, full_size)
%
% Loads the cached Sx, Sy (single, ROI-cropped) if frame_NNNN.mat exists, else
% computes them with cisg_calculate_slopes and, when cache_dir is non-empty,
% saves them.  Safe inside parfor (each frame writes its own file).  full_size = [h w] of a
% good full-resolution frame (default: from ref_b); a bad file raises 'cisg:badFrame'.

    f = '';
    if ~isempty(cache_dir), f = fullfile(cache_dir, sprintf('frame_%04d.mat', n)); end
    if ~isempty(f) && isfile(f)
        D = load(f, 'Sx', 'Sy');  Sx = D.Sx;  Sy = D.Sy;
        return;
    end
    if nargin < 8 || isempty(full_size), full_size = bin * [size(ref_b, 1) size(ref_b, 2)]; end
    obs = cisg_load_coreview(frame_file(n), full_size);          % throws 'cisg:badFrame' on a bad file
    if bin > 1, obs = uint8(imresize(obs, 1 / bin, 'box')); end
    [sx, sy] = cisg_calculate_slopes(ref_b, obs, setup_b);
    r = setup_b.roi;
    Sx = single(sx(r.row_min:r.row_max, r.col_min:r.col_max));
    Sy = single(sy(r.row_min:r.row_max, r.col_min:r.col_max));
    if ~isempty(f)
        t = (n - 1) * dt;
        save(f, 'Sx', 'Sy', 'n', 't', '-v7');
    end
end
