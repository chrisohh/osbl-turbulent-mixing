function eta = cisg_reconstruct_eta(Sx, Sy, dx, dy, varargin)
% CISG_RECONSTRUCT_ETA  Spectral integration of (Sx, Sy) -> eta(y, x, t).
%
%   eta = cisg_reconstruct_eta(Sx, Sy, dx, dy)
%   eta = cisg_reconstruct_eta(Sx, Sy, dx, dy, 'verify', true)
%
% Inputs:
%   Sx, Sy : [Ny x Nx x Nt] dimensionless surface slopes.
%   dx, dy : grid spacing in metres.
%
% Output:
%   eta    : [Ny x Nx x Nt] surface elevation in metres (single).
%
% Method:
%   eta_hat(kx, ky) = (-i kx * Sx_hat - i ky * Sy_hat) / (kx^2 + ky^2),
%   with the k = 0 mode set to zero. Per-frame DC of Sx, Sy is removed
%   before the FFT.

p = inputParser;
addParameter(p, 'verify', false, @islogical);
parse(p, varargin{:});
verify_flag = p.Results.verify;

[Ny, Nx, Nt] = size(Sx);

kx = 2*pi * ifftshift(((-Nx/2):(Nx/2-1)) / (Nx*dx));   % rad/m
ky = 2*pi * ifftshift(((-Ny/2):(Ny/2-1)) / (Ny*dy));
[KX, KY] = meshgrid(kx, ky);
K2 = KX.^2 + KY.^2;
denom_inv = 1 ./ K2;
denom_inv(1, 1) = 0;                                   % zero out DC

eta = zeros(Ny, Nx, Nt, 'single');
for it = 1:Nt
    sxf = double(Sx(:,:,it));  sxf = sxf - mean(sxf(:));
    syf = double(Sy(:,:,it));  syf = syf - mean(syf(:));
    Sxh = fft2(sxf);
    Syh = fft2(syf);
    eta_h = (-1i*KX .* Sxh + -1i*KY .* Syh) .* denom_inv;
    eta(:,:,it) = single(real(ifft2(eta_h)));
end

if verify_flag
    it = round(Nt/2);
    eta_v = double(eta(:,:,it));
    [eta_x, eta_y] = gradient(eta_v, dx, dy);
    figure('Name', 'cisg_reconstruct_eta verify', 'Position', [100 100 1400 400]);
    subplot(1,3,1);
    imagesc(double(Sx(:,:,it))); axis equal tight; colorbar;
    title(sprintf('Original S_x (frame %d)', it));
    subplot(1,3,2);
    imagesc(eta_x); axis equal tight; colorbar;
    title('\partial\eta/\partialx');
    subplot(1,3,3);
    diff_xx = double(Sx(:,:,it)) - eta_x;
    imagesc(diff_xx); axis equal tight; colorbar;
    title(sprintf('S_x - \\partial\\eta/\\partialx  (RMS=%.3g)', rms(diff_xx(:))));
    drawnow;
end
end
