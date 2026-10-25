function [Sx, Sy, info] = cisg_cache_frame(n, ref_b, setup_b, frame_file, cache_dir, bin, dt, full_size, warm)
%CISG_CACHE_FRAME Slopes of one frame inside the ROI, cached as frame_NNNN.mat (same format as cisg_main).
%
%   [Sx, Sy] = cisg_cache_frame(n, ref_b, setup_b, frame_file, cache_dir, bin, dt, full_size)
%   [Sx, Sy, info] = cisg_cache_frame(..., warm)        % warm start from the previous frame
%
% Loads the cached Sx, Sy (single, ROI-cropped) if frame_NNNN.mat exists, else
% computes them with cisg_calculate_slopes and, when cache_dir is non-empty,
% saves them.  Safe inside parfor (each frame writes its own file).  full_size = [h w] of a
% good full-resolution frame (default: from ref_b); a bad file raises 'cisg:badFrame'.
% A cached outlier mask (variable 'out', from cisg_outlier_pass) is applied on load.  setup_b.bias_Sx / bias_Sy
% (cisg_flat_bias) are subtracted from computed slopes before saving.
%
% warm (optional struct; [] = cold start, the old behaviour):
%   prev   struct n, Sx, Sy of the last frame the caller solved ([] = none); used when prev.n = n - 1
%   lead   frames solved (not saved) before n when no frame n - 1 is available, each warm-started from the
%          one before: n - lead ... n - 1 -> n.  Frame n - 1 from the cache is used first if it exists.
% The colour solver is started both from the previous frame's shift and from 0, per pixel the better
% solution is kept (cisg_colour_displacement, init_dx).  info: warm_fraction, lead_used, prior_from.

    if nargin < 9, warm = []; end
    info = struct('warm_fraction', 0, 'lead_used', 0, 'prior_from', 'none');
    f = '';
    if ~isempty(cache_dir), f = fullfile(cache_dir, sprintf('frame_%04d.mat', n)); end
    if ~isempty(f) && isfile(f)
        [Sx, Sy] = load_cached(f);
        info.prior_from = 'cache';
        return;
    end
    if nargin < 8 || isempty(full_size), full_size = bin * [size(ref_b, 1) size(ref_b, 2)]; end

    % previous-frame slopes for the warm start
    prior = [];
    if ~isempty(warm)
        if isfield(warm, 'prev') && ~isempty(warm.prev) && warm.prev.n == n - 1
            prior = warm.prev;  info.prior_from = 'caller';
        else
            fp = '';  if ~isempty(cache_dir), fp = fullfile(cache_dir, sprintf('frame_%04d.mat', n - 1)); end
            if ~isempty(fp) && isfile(fp)
                [px, py] = load_cached(fp);  prior = struct('n', n - 1, 'Sx', px, 'Sy', py);  info.prior_from = 'cache';
            elseif isfield(warm, 'lead') && warm.lead > 0 && n > 1
                for m = max(1, n - warm.lead):n - 1                  % lead-in chain, not saved
                    try
                        [px, py] = solve(m, prior, ref_b, setup_b, frame_file, bin, full_size);
                        prior = struct('n', m, 'Sx', px, 'Sy', py);
                    catch ME
                        if ~strcmp(ME.identifier, 'cisg:badFrame'), rethrow(ME); end
                        prior = [];                                  % bad lead frame: restart the chain
                    end
                end
                info.lead_used = min(warm.lead, n - 1);
                if ~isempty(prior), info.prior_from = 'lead'; end
            end
        end
    end

    [Sx, Sy, ci] = solve(n, prior, ref_b, setup_b, frame_file, bin, full_size);
    if isfield(ci, 'warm_fraction'), info.warm_fraction = ci.warm_fraction; end
    r = setup_b.roi;
    sgx4 = [];  sgy4 = [];
    if ~isempty(ci) && isfield(ci, 'sig_sx') && ~isempty(ci.sig_sx)           % slope uncertainty (precision mask), every 4th px
        sgx4 = single(ci.sig_sx(1:4:end, 1:4:end));  sgy4 = single(ci.sig_sy(1:4:end, 1:4:end));
    end
    if ~isempty(f)
        t = (n - 1) * dt;
        warm_fraction = single(info.warm_fraction);
        vars = {'Sx', 'Sy', 'n', 't', 'warm_fraction'};
        if ~isempty(ci) && isfield(ci, 'dye_c') && ~isempty(ci.dye_c)      % dye amount (joint model), every 4th px
            c4 = single(ci.dye_c(r.row_min:4:r.row_max, r.col_min:4:r.col_max));
            a4 = uint8(ci.dye_applied(r.row_min:4:r.row_max, r.col_min:4:r.col_max));   % where the joint correction was applied
            vars = [vars {'c4', 'a4'}];
        end
        if ~isempty(sgx4), vars = [vars {'sgx4', 'sgy4'}]; end
        save(f, vars{:}, '-v7');
    end
end

function [Sx, Sy, info] = solve(n, prior, ref_b, setup_b, frame_file, bin, full_size)
    obs = cisg_load_coreview(frame_file(n), full_size);          % throws 'cisg:badFrame' on a bad file
    if bin > 1, obs = uint8(imresize(obs, 1 / bin, 'box')); end
    if ~isempty(prior), setup_b.init_Sx = prior.Sx;  setup_b.init_Sy = prior.Sy; end
    [sx, sy, info] = cisg_calculate_slopes(ref_b, obs, setup_b);
    r = setup_b.roi;
    Sx = single(sx(r.row_min:r.row_max, r.col_min:r.col_max));
    Sy = single(sy(r.row_min:r.row_max, r.col_min:r.col_max));
    if isfield(setup_b, 'bias_Sx') && ~isempty(setup_b.bias_Sx)          % fixed flat-water bias (cisg_flat_bias)
        Sx = Sx - setup_b.bias_Sx;  Sy = Sy - setup_b.bias_Sy;
    end
end

function [Sx, Sy] = load_cached(f)
    D = load(f);
    Sx = D.Sx;  Sy = D.Sy;
    if isfield(D, 'out') && ~isempty(D.out), Sx(D.out) = NaN;  Sy(D.out) = NaN; end
end
