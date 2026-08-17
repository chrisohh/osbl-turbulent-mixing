function [w_wave, info] = wave_w_from_eta(eta_x, dx, z_depth, opt)
%WAVE_W_FROM_ETA  Vertical wave orbital velocity profile from a surface line.
%
% Implements the eta-based (Hilbert) wave-removal of OSBL-notes/wave_removal.tex
% for TRANSVERSE PIV: the wave field is characterised from the LONGITUDINAL
% surface eta(x) at the matched time, then evaluated as a depth profile w_wave(z)
% that is uniform across the cross-wind (y) direction.
%
% Steps (per wave_removal.tex §"Compute w_wave"):
%   1. Hilbert transform of eta(x) -> analytic signal  eta_hat = a e^{i phi}.
%   2. Smooth the unwrapped phase, k(x) = d phi / dx  (local wavenumber).
%   3. Gravity-capillary dispersion (deep water):
%         omega = sqrt( k (g + gamma k^2) ),  gamma = 7.2e-5 m^3 s^-2.
%   4. w_wave(z) = a * omega * exp(k z) * sin(phi),  evaluated at x = opt.x_eval.
%
% Inputs
%   eta_x    surface elevation eta(x) at one time, metres, sampled along x
%   dx       x sample spacing, metres
%   z_depth  depths to evaluate at (<= 0, relative to mean surface), metres [Nz x 1]
%   opt      optional struct:
%              .x_eval  x location to evaluate, metres (default: 12 m, clamped)
%              .gamma   surface-tension param, m^3 s^-2 (default 7.2e-5)
%              .g       gravity, m/s^2 (default 9.81)
%              .phase_smooth  movmean window (samples) for phase (default 11)
%
% Outputs
%   w_wave   vertical orbital velocity at each z_depth, m/s [Nz x 1]
%   info     struct with a, k, omega, phi, ix actually used (diagnostics)

if nargin < 4, opt = struct(); end
if ~isfield(opt, 'x_eval'),       opt.x_eval = 12;       end   % m (wave_removal.tex)
if ~isfield(opt, 'gamma'),        opt.gamma  = 7.2e-5;   end   % m^3 s^-2
if ~isfield(opt, 'g'),            opt.g      = 9.81;     end   % m/s^2
if ~isfield(opt, 'phase_smooth'), opt.phase_smooth = 11; end   % samples

eta_x   = eta_x(:);
z_depth = z_depth(:);

% 1. analytic signal
eta_d    = eta_x - mean(eta_x, 'omitnan');
analytic = hilbert(eta_d);
a        = abs(analytic);
phi      = unwrap(angle(analytic));

% 2. smooth phase, local wavenumber
phi_s = movmean(phi, opt.phase_smooth);
k     = gradient(phi_s, dx);          % rad/m

% evaluation point (x = opt.x_eval), clamped to the available range
N  = numel(eta_x);
ix = round(opt.x_eval / dx) + 1;
ix = min(max(ix, 1), N);

a0   = a(ix);
k0   = k(ix);
phi0 = phi_s(ix);

if ~(k0 > 0)
    warning('wave_w_from_eta:nonPositiveK', ...
            'Local wavenumber at x_eval is %.3g (<=0); w_wave set to 0.', k0);
    w_wave = zeros(size(z_depth));
else
    % 3. gravity-capillary dispersion
    omega0 = sqrt( k0 * (opt.g + opt.gamma * k0^2) );
    % 4. depth profile (z_depth <= 0 -> exp(k z) decays with depth)
    w_wave = a0 * omega0 * exp(k0 * z_depth) * sin(phi0);
end

info = struct('a', a0, 'k', k0, 'omega', ...
              (k0>0) * sqrt(max(k0,0) * (opt.g + opt.gamma * max(k0,0)^2)), ...
              'phi', phi0, 'ix', ix);
end
