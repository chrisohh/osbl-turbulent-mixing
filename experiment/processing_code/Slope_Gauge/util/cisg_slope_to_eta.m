function [eta_mm, info] = cisg_slope_to_eta(Sx, Sy, mm_per_px_x, mm_per_px_y, q)
%CISG_SLOPE_TO_ETA Least-squares integration of one slope frame to a relative elevation.
%
%   [eta_mm, info] = cisg_slope_to_eta(Sx, Sy, mm_per_px_x, mm_per_px_y, q)
%
% Sx, Sy        - [rows x cols] slopes, NaN allowed.  Sx = d eta/dx with x = image RIGHT
%                 (columns), Sy = d eta/dy with y = image UP (decreasing row).
% mm_per_px_*   - surface mm per pixel along columns / rows (cisg_local_scale * 10)
% q             - block size (px) for averaging before the solve (default 16)
%
% eta_mm        - elevation (mm) on the coarse grid [floor(rows/q) x floor(cols/q)], mean removed
% info          - Sx_c, Sy_c (block-averaged slopes, NaN-filled), dx_mm, dy_mm (coarse spacing),
%                 nan_fraction (of coarse cells that had no valid pixel), valid_frac (per coarse cell:
%                 fraction of pixels with both slopes valid).  eta_mm is defined everywhere (the gaps are
%                 filled for the solve); blank the filled cells with eta_mm(info.valid_frac < 0.5) = NaN.
%
% Non-periodic Poisson solve with Neumann boundaries: the slope field is
% mirror-extended (eta even, so the slope normal to a mirror is odd), then
% integrated spectrally, eta_hat = -i (kx Sx_hat + kr S_row_hat) / k^2 (Frankot-
% Chellappa), which is a DCT solve on the original quadrant.  Missing slopes are
% filled by Laplace relaxation first.
%
% eta is RELATIVE (mean removed) and scales longer than the ROI are unreliable.
% (cisg_reconstruct_eta.m integrates a periodic cube; this one is for single,
% non-periodic frames with gaps.)

    if nargin < 5, q = 16; end
    [nr0, nc0] = size(Sx);
    nr = floor(nr0 / q);  nc = floor(nc0 / q);
    blk = @(A) squeeze(mean(reshape(double(A(1:nr*q, 1:nc*q)), q, nr, q, nc), [1 3], 'omitnan'));
    Sxc = blk(Sx);  Syc = blk(Sy);

    bad = isnan(Sxc) | isnan(Syc);
    nan_frac = mean(bad(:));
    Sxc = relax_fill(Sxc, bad);  Syc = relax_fill(Syc, bad);

    dx = q * mm_per_px_x;  dr = q * mm_per_px_y;
    Srow = -Syc;                                         % d eta / d row (rows run down, y up)

    % mirror extension to 2nr x 2nc: eta even => S normal to the mirror odd
    ext = @(A, sr, sc) [A, sc * A(:, end:-1:1); sr * A(end:-1:1, :), sr * sc * A(end:-1:1, end:-1:1)];
    Sxe = ext(Sxc, 1, -1);
    Sre = ext(Srow, -1, 1);
    [ny, nx] = size(Sxe);
    kx = 2*pi * [0:nx/2, -nx/2+1:-1] / (nx * dx);
    kr = 2*pi * [0:ny/2, -ny/2+1:-1]' / (ny * dr);
    K2 = kx.^2 + kr.^2;  K2(1,1) = 1;
    eh = (-1i * kx .* fft2(Sxe) - 1i * kr .* fft2(Sre)) ./ K2;
    eh(1,1) = 0;
    e = real(ifft2(eh));
    eta_mm = e(1:nr, 1:nc);
    eta_mm = eta_mm - mean(eta_mm(:));

    vf = blk(isfinite(Sx) & isfinite(Sy));
    info = struct('Sx_c', Sxc, 'Sy_c', Syc, 'dx_mm', dx, 'dy_mm', dr, 'nan_fraction', nan_frac, 'valid_frac', vf);
end

function A = relax_fill(A, bad)
    % Laplace relaxation of the NaN cells (valid cells fixed)
    if ~any(bad(:)), return; end
    A(bad) = mean(A(~bad));
    for it = 1:400
        P = padarray(A, [1 1], 'replicate');
        N = (P(1:end-2, 2:end-1) + P(3:end, 2:end-1) + P(2:end-1, 1:end-2) + P(2:end-1, 3:end)) / 4;
        A(bad) = N(bad);
    end
end
