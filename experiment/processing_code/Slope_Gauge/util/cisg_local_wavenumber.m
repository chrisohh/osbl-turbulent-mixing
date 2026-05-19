function [a, kx_loc, ky_loc, phi] = cisg_local_wavenumber(eta, dx, dy)
% CISG_LOCAL_WAVENUMBER  2D Riesz transform → local amplitude, wavenumber.
%
%   [a, kx_loc, ky_loc, phi] = cisg_local_wavenumber(eta, dx, dy)
%
% Inputs:
%   eta : [Ny x Nx x Nt] surface elevation (m).
%   dx, dy : grid spacing (m).
%
% Outputs (all [Ny x Nx x Nt]):
%   a       : envelope amplitude (m)
%   kx_loc  : local along-x wavenumber (rad/m), floored at 20
%   ky_loc  : local along-y wavenumber (rad/m)
%   phi     : monogenic phase (rad), unwrapped along x then y
%
% Notes:
%   Riesz kernels in spectral space:
%       R_x_hat = i*kx/|k| * eta_hat
%       R_y_hat = i*ky/|k| * eta_hat
%   The monogenic signal envelope is sqrt(eta^2 + R_x^2 + R_y^2);
%   orientation-aware phase = atan2(sqrt(R_x^2 + R_y^2), eta).
%   Local wavenumber = grad(phi) (unwrapped, lightly smoothed).

[Ny, Nx, Nt] = size(eta);

kx = 2*pi * ifftshift(((-Nx/2):(Nx/2-1)) / (Nx*dx));
ky = 2*pi * ifftshift(((-Ny/2):(Ny/2-1)) / (Ny*dy));
[KX, KY] = meshgrid(kx, ky);
Kmag = sqrt(KX.^2 + KY.^2);
Kmag(1, 1) = 1;                              % avoid 0/0 at DC

Hx = 1i * KX ./ Kmag;  Hx(1,1) = 0;          % Riesz kernels
Hy = 1i * KY ./ Kmag;  Hy(1,1) = 0;

a      = zeros(Ny, Nx, Nt, 'single');
kx_loc = zeros(Ny, Nx, Nt, 'single');
ky_loc = zeros(Ny, Nx, Nt, 'single');
phi    = zeros(Ny, Nx, Nt, 'single');

for it = 1:Nt
    e = double(eta(:,:,it));
    e = e - mean(e(:));
    Eh = fft2(e);
    Rx = real(ifft2(Hx .* Eh));
    Ry = real(ifft2(Hy .* Eh));
    a_t   = sqrt(e.^2 + Rx.^2 + Ry.^2);
    phi_t = atan2(sqrt(Rx.^2 + Ry.^2), e);
    % unwrap along x then y
    phi_t = unwrap(phi_t, [], 2);
    phi_t = unwrap(phi_t, [], 1);
    % light smoothing of phase for derivatives
    phi_s = local_smooth2(phi_t, 5);
    [kx_t, ky_t] = gradient(phi_s, dx, dy);
    kx_t = max(kx_t, 20);                    % floor at 20 rad/m
    a(:,:,it)      = single(a_t);
    phi(:,:,it)    = single(phi_t);
    kx_loc(:,:,it) = single(kx_t);
    ky_loc(:,:,it) = single(ky_t);
end
end


function ys = local_smooth2(y, n)
% Lightweight 2D moving average smoother (n-point box) — avoids dependency
% on the smoothn helper that lives in the PIV folder.
if n < 2, ys = y; return; end
k = ones(n) / n^2;
ys = conv2(y, k, 'same');
end
