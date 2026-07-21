function omega = dispersion_curves(k, mode, U_E, U_S)
% DISPERSION_CURVES  Linear gravity-capillary dispersion ± Doppler shift.
%
%   omega = dispersion_curves(k, 'still')
%   omega = dispersion_curves(k, 'eulerian',   U_E)
%   omega = dispersion_curves(k, 'lagrangian', U_E, U_S)
%
% Inputs:
%   k    : wavenumber array (rad/m), positive
%   mode : 'still' | 'eulerian' | 'lagrangian'
%   U_E  : Eulerian surface drift (m/s); pass NaN to skip Doppler
%   U_S  : Stokes drift at surface (m/s)
%
% Output:
%   omega : angular frequency (rad/s) on the same grid as k.
%
% gamma = sigma/rho = 7.4e-5 m^3/s^2 (clean water at ~20 C).

g     = 9.81;
gamma = 7.4e-5;

if nargin < 3, U_E = NaN; end
if nargin < 4, U_S = NaN; end

omega_still = sqrt(g .* k + gamma .* k.^3);

switch lower(mode)
    case 'still'
        omega = omega_still;
    case 'eulerian'
        if isnan(U_E), omega = nan(size(k)); return; end
        omega = omega_still + k .* U_E;
    case 'lagrangian'
        if isnan(U_E) || isnan(U_S), omega = nan(size(k)); return; end
        omega = omega_still + k .* (U_E + U_S);
    otherwise
        error('dispersion_curves: unknown mode "%s"', mode);
end
end
