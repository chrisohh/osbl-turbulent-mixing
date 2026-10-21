function [Sx, Sy] = cisg_cache_frame(n, ref_b, setup_b, frame_file, cache_dir, bin, dt)
%CISG_CACHE_FRAME Slopes of one frame inside the ROI, cached as frame_NNNN.mat (same format as cisg_main).
%
%   [Sx, Sy] = cisg_cache_frame(n, ref_b, setup_b, frame_file, cache_dir, bin, dt)
%
% Loads the cached Sx, Sy (single, ROI-cropped) if frame_NNNN.mat exists, else
% computes them with cisg_calculate_slopes and, when cache_dir is non-empty,
% saves them.  Safe inside parfor (each frame writes its own file).

    f = '';
    if ~isempty(cache_dir), f = fullfile(cache_dir, sprintf('frame_%04d.mat', n)); end
    if ~isempty(f) && isfile(f)
        D = load(f, 'Sx', 'Sy');  Sx = D.Sx;  Sy = D.Sy;
        return;
    end
    obs = cisg_load_coreview(frame_file(n));
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
