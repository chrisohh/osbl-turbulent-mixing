function A = lif_img_filter(A, kind, sigma_px)
%LIF_IMG_FILTER Noise filter for LIF frames before theta / surface detection.
%
%   A = LIF_IMG_FILTER(A, kind, sigma_px)
%
% kind:
%   'none'     - unchanged
%   'median3'  - 3x3 median: removes isolated noisy pixels, keeps edges and
%                the amplitude of anything wider than ~3 px
%   'median5'  - 5x5 median
%   'wiener5'  - 5x5 adaptive Wiener (wiener2): strong where the image is
%                flat, weak where there is structure, from the local variance
%   'gauss'    - Gaussian low-pass, sigma = sigma_px (keeps the integral,
%                softens edges); sigma_px <= 0 -> unchanged

    switch kind
        case 'none'
        case 'median3', A = medfilt2(A, [3 3], 'symmetric');
        case 'median5', A = medfilt2(A, [5 5], 'symmetric');
        case 'wiener5', A = wiener2(A, [5 5]);
        case 'gauss'
            if nargin >= 3 && sigma_px > 0, A = imgaussfilt(A, sigma_px); end
        otherwise
            error('lif_img_filter:kind', 'Unknown filter ''%s''.', kind);
    end
end
