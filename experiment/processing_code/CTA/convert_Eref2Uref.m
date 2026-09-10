function [U_ref, cal] = convert_Eref2Uref(E_ref, cert)
% CONVERT_EREF2UREF  Convert Dantec 54T29 reference-probe voltage to velocity.
%
%   U_ref = convert_Eref2Uref(E_ref)          % built-in certificate T29-202
%   U_ref = convert_Eref2Uref(E_ref, cert)    % a different certificate
%   [U_ref, cal] = convert_Eref2Uref(...)     % also return what was used
%
%   E_ref - reference-probe voltage (V), any shape; shape is preserved.
%   cert  - (optional) struct with fields E and U, the velocity/voltage table
%           from that probe's calibration certificate. Defaults to the
%           built-in T29-202 table below.
%   cal   - struct describing the conversion actually applied, including
%           counts of samples clamped at either end.
%
% METHOD
%   The 54T29 is a LINEARIZED velocity reference transducer -- its output is
%   internally conditioned, so it does not follow King's law and there is no
%   compact analytic form that reproduces the certificate. Fitting was tested
%   against the T29-202 table (18 points, 0.096-30.6 m/s):
%
%     monotone (pchip) interpolation   1.26 % RMS   (leave-one-out)
%     4th-order poly in log(U)         3.41 % RMS   (leave-one-out)
%     King's law E^2 = A + B*U^n      13.2  % RMS
%
%   So the certificate table IS the calibration and is interpolated directly
%   with PCHIP -- shape-preserving, so it cannot introduce the oscillations a
%   high-order polynomial does between widely spaced points, and it reproduces
%   every tabulated point exactly. Interpolation error is comfortably inside
%   the certificate's own stated uncertainty (+/-2% or +/-0.02 m/s, k=2).
%
% OUT-OF-RANGE HANDLING
%   Voltages outside the tabulated span are CLAMPED, never extrapolated, and
%   counted in cal.nBelow / cal.nAbove. Extrapolating this curve is what broke
%   the previous implementation (see HISTORY). Check those counts: a nonzero
%   nAbove means the flow exceeded the probe's calibrated range and those
%   samples are ceiling values, not measurements.
%
% TEMPERATURE / PRESSURE
%   NOT APPLIED. The certificate records its own test conditions (25.1 degC,
%   101.7 kPa, output referenced to 101.325 kPa standard pressure) in
%   cal.T_cal_C / cal.P_cal_kPa, but Dantec does not publish a correction
%   formula for this transducer, so E_ref is used as measured. If the lab
%   temperature differs materially from 25.1 degC, that is an uncorrected
%   systematic error -- log ambient temperature and treat it as a known
%   uncertainty rather than assuming it away.
%
% HISTORY -- why this function was rewritten
%   The previous version modelled the probe as a piecewise-sqrt curve with
%   parameters (G1, G2, U0..U3) that did not come from this probe's
%   certificate. Checked against T29-202 it was wrong by a factor of 156 at
%   the low end and returned NEGATIVE velocities above 4.0 V (RMS error
%   37.8 m/s over a 0-30 m/s range). Its two structural errors were assuming
%   zero volts at zero velocity -- the certificate shows ~0.52 V at 0.096 m/s
%   -- and fitting a 4th-order polynomial over 0..2.64 V, then evaluating it
%   out at 4.8 V where it had already turned over. Any U_ref computed before
%   this rewrite should be regarded as meaningless. Passing the old
%   (G1/G2/U0..U3) struct is therefore rejected outright rather than silently
%   reinterpreted.

    % --- Built-in certificate -------------------------------------------
    % Dantec Dynamics certificate T29-202, probe 54T29 serial 0202.
    % Certificate date 250617 (YYMMDD) = "Calibrated 170625" (DDMMYY) =
    % 17 June 2025 -- the two fields are the same date in different formats,
    % not a year discrepancy.
    % Procedure Dantec 5900F214. Approved P. Nielsen (PNN), Skovlunde, DK.
    % Range 0-30 m/s; transfer function file T29_0202.ref.
    % Test conditions: atmospheric air, 25.1 degC, 101.7 kPa, 40% RH.
    CERT = struct( ...
        'id',    'T29-202', ...
        'probe', '54T29', ...
        'serial','0202', ...
        'date',  '2025-06-17', ...
        'T_cal_C',   25.1, ...
        'P_cal_kPa', 101.7, ...
        'P_ref_kPa', 101.325, ...
        'RH_pct',    40, ...
        'U', [0.0958 0.1396 0.2013 0.2832 0.5124 0.8113 1.6710 2.4340 3.3420 ...
              4.6090 6.1020 7.9970 10.2900 12.9300 16.4000 20.3800 25.2400 30.6000], ...
        'E', [0.5196 0.6122 0.7304 0.8645 1.1244 1.3475 1.8125 2.1113 2.3849 ...
              2.6651 2.8977 3.1526 3.4219 3.6940 3.9939 4.2833 4.5752 4.8348]);

    if nargin < 2 || isempty(cert)
        cert = CERT;
    elseif isfield(cert, 'G1') || isfield(cert, 'U3')
        error('convert_Eref2Uref:legacyCal', ...
            ['This function no longer accepts the (G1,G2,U0..U3) struct -- that ' ...
             'model was wrong for this probe by up to 156x and returned negative ' ...
             'velocities above 4 V. Call convert_Eref2Uref(E_ref) to use the ' ...
             'built-in T29-202 certificate, or pass a struct with .E and .U ' ...
             'fields from your probe''s own certificate.']);
    elseif ~isfield(cert, 'E') || ~isfield(cert, 'U')
        error('convert_Eref2Uref:badCert', ...
            'cert must be a struct with fields E and U (certificate table).');
    end

    Ec = cert.E(:);
    Uc = cert.U(:);
    if numel(Ec) ~= numel(Uc)
        error('convert_Eref2Uref:certSize', 'cert.E and cert.U must be the same length.');
    end
    [Ec, ord] = sort(Ec);   % pchip needs ascending abscissae
    Uc = Uc(ord);

    sz    = size(E_ref);
    E_ref = double(E_ref(:));

    % --- Clamp, then interpolate ----------------------------------------
    below = E_ref < Ec(1);
    above = E_ref > Ec(end);
    Eclamped = min(max(E_ref, Ec(1)), Ec(end));

    U_ref = interp1(Ec, Uc, Eclamped, 'pchip');
    U_ref = reshape(U_ref, sz);

    cal = cert;
    cal.method = 'pchip interpolation of certificate table';
    cal.E_min  = Ec(1);    cal.U_min = Uc(1);
    cal.E_max  = Ec(end);  cal.U_max = Uc(end);
    cal.nBelow = sum(below);
    cal.nAbove = sum(above);
    cal.nTotal = numel(E_ref);

    if cal.nAbove > 0
        warning('convert_Eref2Uref:aboveRange', ...
            ['%d of %d samples (%.1f%%) exceed the calibrated range (%.4f V = ' ...
             '%.2f m/s) and were clamped -- these are NOT measurements.'], ...
            cal.nAbove, cal.nTotal, 100*cal.nAbove/cal.nTotal, cal.E_max, cal.U_max);
    end
    if cal.nBelow > 0
        % Expected whenever the fan is off: the probe floor is 0.096 m/s, so
        % genuine still air sits below the table. Informational, not a warning.
        fprintf(['convert_Eref2Uref: %d of %d samples (%.1f%%) below the ' ...
                 'calibrated floor (%.4f V = %.4f m/s), clamped.\n'], ...
            cal.nBelow, cal.nTotal, 100*cal.nBelow/cal.nTotal, cal.E_min, cal.U_min);
    end
end
