function [U, V, W, Ucal1, Ucal2, Ucal3] = convert_E2U_fn(E1, E2, E3, cal_file, Ta)
% CONVERT_E2U_FN  Convert tri-axial hot-wire voltages to velocities.
%
%   [U, V, W] = convert_E2U_fn(E1, E2, E3, cal_file)
%   [U, V, W, Ucal1, Ucal2, Ucal3] = convert_E2U_fn(E1, E2, E3, cal_file, Ta)
%
%   Function form of the pipeline in convertE2U.m, for use from acquisition
%   scripts (run_experiment / a hot-wire test) so live-acquired voltages are
%   reduced to velocities exactly like the offline processing.
%
%   E1,E2,E3 - column vectors of raw sensor voltages (V), one per sensor.
%              From setup_hotwire_daq these are data(:,1), data(:,2), data(:,3).
%   cal_file - path to the StreamWare calibration/header .txt for THIS probe.
%   Ta       - (optional) ambient/fluid temperature during acquisition, degC.
%              If given, applies the Dantec temperature correction; if omitted
%              or empty, no correction is applied (raw voltages used directly).
%
%   U,V,W    - velocity components in probe coordinates (m/s), same length as E1.
%   Ucal1..3 - per-sensor calibration velocities (m/s), for debugging/validation.
%
% NOTE: depends on parse_calibration.m (same CTA folder). When calling from
% another folder, add this folder to the path first, e.g.:
%   addpath('D:\Chris\osbl-turbulent-mixing\experiment\processing_code\CTA');
%
% NOT VERIFIED: OVERHEAT_RATIO / TCR_ALPHA below are hardcoded for this probe
% (a = 0.8, alpha = 0.46 %/degC). Update if you use a different probe/overheat.

    % Column vectors
    E1 = E1(:); E2 = E2(:); E3 = E3(:);

    cal = parse_calibration(cal_file);

    %% Optional temperature correction (Dantec: Ecorr = ((Tw-T0)/(Tw-Ta))^0.5 * Ea)
    if nargin >= 5 && ~isempty(Ta)
        OVERHEAT_RATIO = 0.8;      % a
        TCR_ALPHA      = 0.46/100; % alpha, 1/degC
        T_w = 20 + OVERHEAT_RATIO / TCR_ALPHA;      % Tw = 20degC + a/alpha
        f   = sqrt((T_w - cal.T_ref) ./ (T_w - Ta)); % same factor for all sensors
        E1  = E1 .* f;  E2 = E2 .* f;  E3 = E3 .* f;
    end

    %% Linearization (polynomial in-range, origin-to-floor line below range)
    U_floor = cal.U_floor;   % per-sensor min calibration velocity (m/s)
    E_floor = cal.E_floor;   % per-sensor min calibration voltage (V)

    Ucal1 = cal.C0_1 + cal.C1_1*E1 + cal.C2_1*E1.^2 + cal.C3_1*E1.^3 + cal.C4_1*E1.^4;
    b1 = E1 < E_floor(1);  Ucal1(b1) = (U_floor(1) / E_floor(1)) * E1(b1);

    Ucal2 = cal.C0_2 + cal.C1_2*E2 + cal.C2_2*E2.^2 + cal.C3_2*E2.^3 + cal.C4_2*E2.^4;
    b2 = E2 < E_floor(2);  Ucal2(b2) = (U_floor(2) / E_floor(2)) * E2(b2);

    Ucal3 = cal.C0_3 + cal.C1_3*E3 + cal.C2_3*E3.^2 + cal.C3_3*E3.^3 + cal.C4_3*E3.^4;
    b3 = E3 < E_floor(3);  Ucal3(b3) = (U_floor(3) / E_floor(3)) * E3(b3);

    %% Decomposition into wire-coordinate velocities (U1, U2, U3)
    cos_angle = cosd(35.3);
    A = [cal.k1_sq, 1,         cal.h1_sq;
         cal.h2_sq, cal.k2_sq, 1;
         1,         cal.h3_sq, cal.k3_sq];

    % Solve A * Usq = B for all samples at once (A is constant; B is 3 x N)
    B = [ (Ucal1.^2 * (1 + cal.k1_sq + cal.h1_sq) * cos_angle^2)';
          (Ucal2.^2 * (1 + cal.k2_sq + cal.h2_sq) * cos_angle^2)';
          (Ucal3.^2 * (1 + cal.k3_sq + cal.h3_sq) * cos_angle^2)' ];
    Usq = max(A \ B, 0);   % 3 x N, clamp negatives

    U1 = sqrt(Usq(1, :));
    U2 = sqrt(Usq(2, :));
    U3 = sqrt(Usq(3, :));

    %% Transform to probe coordinates via the probe matrix Mp
    if isfield(cal, 'Mp') && isequal(size(cal.Mp), [3 3])
        Mp = cal.Mp;
    else
        % Fallback: 55R95 default (from probe library header)
        warning('convert_E2U_fn:noMp', ...
            'cal.Mp not found/parsed; using 55R95 default transform matrix.');
        Mp = [ 0.57735   0.57735   0.57735;
              -0.707107  0.707107  6.12323e-17;
               0.408248  0.408248 -0.816496];
    end

    velocity_probe = Mp * [U1; U2; U3];   % 3 x N
    U = velocity_probe(1, :)';
    V = velocity_probe(2, :)';
    W = velocity_probe(3, :)';
end
