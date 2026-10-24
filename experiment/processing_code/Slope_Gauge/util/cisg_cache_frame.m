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
    [sx, sy, info] = cisg_calculate_slopes(ref_b, obs, setup_b);
    r = setup_b.roi;
    Sx = single(sx(r.row_min:r.row_max, r.col_min:r.col_max));
    Sy = single(sy(r.row_min:r.row_max, r.col_min:r.col_max));
    sgx4 = [];  sgy4 = [];
    if ~isempty(info) && isfield(info, 'sig_sx') && ~isempty(info.sig_sx)           % slope uncertainty (precision mask), every 4th px
        sgx4 = single(info.sig_sx(1:4:end, 1:4:end));  sgy4 = single(info.sig_sy(1:4:end, 1:4:end));
    end
    if ~isempty(f)
        t = (n - 1) * dt;
        vars = {'Sx', 'Sy', 'n', 't'};
        if ~isempty(info) && isfield(info, 'dye_c') && ~isempty(info.dye_c)      % dye amount (joint model), every 4th px
            c4 = single(info.dye_c(r.row_min:4:r.row_max, r.col_min:4:r.col_max));
            a4 = uint8(info.dye_applied(r.row_min:4:r.row_max, r.col_min:4:r.col_max));   % where the joint correction was applied
            vars = [vars {'c4', 'a4'}];
        end
        if ~isempty(sgx4), vars = [vars {'sgx4', 'sgy4'}]; end
        save(f, vars{:}, '-v7');
    end
end
