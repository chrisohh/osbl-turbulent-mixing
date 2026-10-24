function nz = cisg_noise_estimate(frame_file, bin, sigma_px, frames, img_size)
%CISG_NOISE_ESTIMATE Camera noise in the smoothed colours used by the inversion, from two early frames.
%
%   nz = cisg_noise_estimate(frame_file, bin, sigma_px, frames, img_size)
%
% Two still-ish frames (default [1 3]: the start of the run) are block-averaged by bin, smoothed with the
% Gaussian of the inversion (sigma_px at that resolution) and subtracted.  The robust spread (1.4826 MAD about
% the median, so slow lighting offsets and the odd wave do not count) of the difference, per channel and in 8
% brightness bins, is exactly the noise of C_obs - C_ref in cisg_colour_displacement (two independent images).
% The bins are fitted by sigma(C) = a + b * C (C = smoothed channel, 0..1).
%
% nz fields: a, b (3 x 1: R G B, colour units 0..1), levels, table (8 x 3 bin sigmas), counts, bin, sigma_px,
%   sigma0 (3 x 1: fit at C = 0.4), frames.

    if nargin < 4 || isempty(frames), frames = [1 3]; end
    if nargin < 5, img_size = [3072 4100]; end
    I = cell(1, 2);
    for k = 1:2
        x = double(cisg_load_coreview(frame_file(frames(k)), img_size)) / 255;
        if bin > 1, x = imresize(x, 1 / bin, 'box'); end
        I{k} = imgaussfilt(x, sigma_px);
    end
    D = I{1} - I{2};  L = (I{1} + I{2}) / 2;
    edges = linspace(0, 1, 9);
    table = nan(8, 3);  counts = zeros(8, 3);
    for c = 1:3
        d = D(:,:,c);  l = L(:,:,c);
        d = d - median(d(:));
        for b = 1:8
            m = l >= edges(b) & l < edges(b + 1);
            counts(b, c) = nnz(m);
            if counts(b, c) > 2000, table(b, c) = 1.4826 * median(abs(d(m))); end
        end
    end
    levels = (edges(1:end-1) + edges(2:end)) / 2;
    a = zeros(3, 1);  bb = zeros(3, 1);
    for c = 1:3
        v = isfinite(table(:, c));
        p = polyfit(levels(v)', table(v, c), 1);
        bb(c) = p(1);  a(c) = p(2);
    end
    nz = struct('a', a, 'b', bb, 'levels', levels, 'table', table, 'counts', counts, 'bin', bin, 'sigma_px', sigma_px, ...
                'sigma0', a + 0.4 * bb, 'frames', frames);
end
