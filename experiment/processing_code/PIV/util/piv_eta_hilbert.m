function info = piv_eta_hilbert(eta, x, t)
% PIV_ETA_HILBERT  Wrapper around compute_wave_orbital_velocities that
% returns only the Hilbert-extracted (a, k, omega, phi) without invoking
% the orbital-velocity reconstruction.
%
% Inputs:
%   eta : [Nx x Nt]  surface elevation (m)
%   x   : [Nx x 1]   along-wind coordinate (m)
%   t   : [Nt x 1]   time (s)
%
% Output:
%   info.a, info.k, info.omega, info.phi  each [Nx x Nt]
%
% Calls compute_wave_orbital_velocities with z = 0 (single depth) and
% discards the velocity outputs.

params.g     = 9.81;
params.gamma = 7.4e-5;

[~, ~, info] = compute_wave_orbital_velocities(eta, x, t, 0, params);
end
