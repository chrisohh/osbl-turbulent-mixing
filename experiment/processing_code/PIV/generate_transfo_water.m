function transfo = generate_transfo_water(eta, x, zeta_vec)
%GENERATE_TRANSFO_WATER  Wave-following coordinate transform for the WATER side.
%
%   Adapted from GenerateTransfo_Fabio_v2.m (air side).
%   Air side:  zeta > 0 upward,   z = +zeta + sum(a*exp(ikx)*exp(-k*zeta))
%   Water side: zeta > 0 downward, z = -zeta + sum(a*exp(ikx)*exp(-k*zeta))
%
%   The exponential decay exp(-k*zeta) is the same: coordinate lines follow
%   the wave at zeta=0 (surface) and flatten as zeta increases (deeper).
%
%   Inputs:
%     eta       [Nx x 1]  surface elevation (m), positive up
%     x         [Nx x 1]  along-wind coordinate (m)
%     zeta_vec  [Nz x 1]  depths below surface (m, positive downward)
%                          e.g. [0; 0.001; 0.002; ... ; 0.14]
%
%   Outputs:
%     transfo.Z_grid   [Nz x Nx]  physical z at each (zeta, x), z=0 at mean
%                                  surface, negative downward
%     transfo.zeta     [Nz x 1]   zeta vector (depth below surface, positive)
%     transfo.x        [Nx x 1]   x vector
%     transfo.eta      [1 x Nx]   surface elevation used
%     transfo.phase    [1 x Nx]   wave phase at each x (from Hilbert)

Nx = length(x);
dx = mean(diff(x));

% --- FFT of surface elevation ---
eta = eta(:)';  % row vector [1 x Nx]
f = fft(eta);

% Symmetric vector: keep only positive frequencies (one-sided spectrum)
Sym_Vec = zeros(1, Nx);
Sym_Vec(1) = 1;
Sym_Vec(Nx/2 + 1) = 1;
Sym_Vec(2:Nx/2) = 2;
% Sym_Vec(Nx/2+2:end) = 0  (already zero)

F = f .* Sym_Vec;

% Wavenumber vector (rad/m, NOT rad/pixel — we work in physical coords)
k_vec = (0:Nx-1) / Nx * (2*pi/dx);

% --- Wave phase from Hilbert transform ---
eta_analytic = hilbert(eta);
phase = angle(eta_analytic);

% --- Build transformation grid ---
Nz = length(zeta_vec);
Z_grid = zeros(Nz, Nx);

for iz = 1:Nz
    zeta = zeta_vec(iz);
    % Water side: z = -zeta + wave undulation * exp(-k*zeta)
    %   At zeta=0: z = eta(x)         → follows the wave
    %   At zeta→∞: z = -zeta → flat   → flattens at depth
    Z_grid(iz,:) = -zeta + real(ifft(F .* exp(-k_vec * zeta), 'nonsymmetric'));
end

% --- Jacobian: J = 1 + d(eta_mapped)/d(zeta) ---
% For water side: J = 1 - sum(a*k*exp(-k*zeta)*cos(kx))
% (positive definite for small steepness)
J = zeros(Nz, Nx);
for iz = 1:Nz
    zeta = zeta_vec(iz);
    J(iz,:) = 1 + real(ifft(F .* k_vec .* exp(-k_vec * zeta), 'nonsymmetric'));
end

% --- Pack output ---
transfo.Z_grid = Z_grid;
transfo.J      = J;
transfo.zeta   = zeta_vec(:);
transfo.x      = x(:);
transfo.eta    = eta;
transfo.phase  = phase;

end
