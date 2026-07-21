function [u_orb, w_orb] = compute_orbitals_linear(eta, x, zeta_vec, params)
%COMPUTE_ORBITALS_LINEAR  Method A: wave orbital velocities from linear theory.
%
%   Uses FFT of the surface elevation to reconstruct orbital velocities
%   in wave-following coordinates. Each spectral component decays as
%   exp(-k * zeta) and oscillates at omega(k) = sqrt(g*k + gamma*k^3).
%
%   Same spectral machinery as generate_transfo_water.m (one-sided FFT,
%   wavenumber vector, exp decay), but multiplied by the dispersion
%   relation factor omega(k).
%
%   Inputs:
%     eta       [1 x Nx] or [Nx x 1]  mean-zero surface elevation (m)
%     x         [Nx x 1]  along-wind coordinate (m)
%     zeta_vec  [Nz x 1]  depth below surface (m, positive downward)
%     params    struct with:
%       .g       gravitational acceleration (m/s^2)
%       .gamma   surface tension / density (m^3/s^2), e.g. 0.073/998
%
%   Outputs:
%     u_orb  [Nz x Nx]  horizontal orbital velocity (m/s, positive in +x)
%     w_orb  [Nz x Nx]  vertical orbital velocity (m/s, positive upward)
%
%   See also: generate_transfo_water

Nx = length(x);
Nz = length(zeta_vec);
dx = mean(diff(x));

eta = eta(:)';  % row [1 x Nx]

% --- One-sided FFT (identical to generate_transfo_water) ---
f = fft(eta);
Sym_Vec = zeros(1, Nx);
Sym_Vec(1) = 1;
Sym_Vec(fix(Nx/2)+1) = 1;
Sym_Vec(2:fix(Nx/2)) = 2;
F = f .* Sym_Vec;

% --- Wavenumber and dispersion ---
k_vec = (0:Nx-1) / Nx * (2*pi/dx);
omega_vec = sqrt(k_vec .* (params.g + params.gamma .* k_vec.^2));

% --- Orbital velocities at each depth ---
u_orb = zeros(Nz, Nx);
w_orb = zeros(Nz, Nx);

for iz = 1:Nz
    zeta = zeta_vec(iz);
    orb = ifft(F .* omega_vec .* exp(-k_vec * zeta), 'nonsymmetric');
    u_orb(iz,:) = real(orb);
    w_orb(iz,:) = imag(orb);
end

end
