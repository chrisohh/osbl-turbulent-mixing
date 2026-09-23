%% plot_cal_velocity.m
% One figure, time vs velocity, holding the three 260909 calibration runs:
%
%   hotwire_cal_20260909_211425   hot-wire, fan 2.1..9 V   -> 55P95 section
%   hotwire_cal_20260909_212437   hot-wire, fan 2.0..3 V   -> 55P95 section
%   hotwire_cal_20260909_172006   reference probe          -> T29 section
%
% Both calibrations come from the same StreamWare export, probe4.txt: the
% tri-axial hot-wire is the 55P95 block (read by parse_calibration.m, which
% takes the first probe in the file), the reference probe is the T29 block
% (read by parse_probe_section.m). The T29 polynomial in probe4.txt is used
% for the reference probe instead of the certificate table built into
% convert_Eref2Uref.m, so all three traces come from one calibration file.
%
% Time base is hwT_wind: 0 = first fan write of that run.
%
% Usage: just run this script (edit the config block first).

clear; clc;

%% ---- Config -------------------------------------------------------------
HW_DIR   = 'D:\HLAB_2026\hotwire';
CAL_FILE = 'C:\Users\airsealab\Documents\GitHub\osbl-turbulent-mixing\experiment\data\260909\probe4.txt';

HW_PROBE  = '55P95';   % probe block in CAL_FILE used for the hot-wire runs
REF_PROBE = 'T29';     % probe block in CAL_FILE used for the reference probe

HW_FILES = { 'hotwire_cal_20260909_211425.mat', ...
             'hotwire_cal_20260909_212437.mat' };
REF_FILE = 'hotwire_cal_20260909_172006_partial.mat';   % the 260909 17:20 ref-probe run

HW_COMPONENT = 'U';    % 'U' = streamwise component, 'speed' = sqrt(U^2+V^2+W^2)
SMOOTH_T     = 0;      % moving-average window in seconds (0 = raw signal)

SAVE_FIG = false;
FIG_DIR  = fullfile(HW_DIR, 'figs');
% -------------------------------------------------------------------------

addpath(fileparts(mfilename('fullpath')));
if SAVE_FIG && ~isfolder(FIG_DIR)
    mkdir(FIG_DIR);
end

% parse_calibration.m reads the FIRST probe block in the file, so the hot-wire
% reduction is only correct if that block is the tri-axial probe.
firstName = regexp(fileread(CAL_FILE), 'Probe name:\s*([^\r\n]*)', 'tokens', 'once');
if ~strcmpi(strtrim(firstName{1}), HW_PROBE)
    error('plot_cal_velocity:wrongFirstProbe', ...
        ['First probe block in %s is "%s", not %s -- convert_E2U_fn would use ' ...
         'the wrong coefficients.'], CAL_FILE, strtrim(firstName{1}), HW_PROBE);
end

figure('Name', 'Calibration runs -- velocity vs time', 'Position', [60 60 1100 600]);
hold on;
colors = [0.000 0.447 0.741;
          0.850 0.325 0.098;
          0.494 0.184 0.556];

%% ---- Hot-wire runs (55P95 block of probe4.txt) --------------------------
for iRun = 1:numel(HW_FILES)
    matFile = fullfile(HW_DIR, HW_FILES{iRun});
    if ~isfile(matFile)
        warning('plot_cal_velocity:missing', 'Not found: %s -- skipped.', matFile);
        continue;
    end
    S = load(matFile);
    [~, tag] = fileparts(HW_FILES{iRun});
    runLabel = erase(tag, 'hotwire_cal_');

    [U, V, W] = convert_E2U_fn(S.E1, S.E2, S.E3, CAL_FILE);
    if strcmpi(HW_COMPONENT, 'speed')
        vel = sqrt(U.^2 + V.^2 + W.^2);
    else
        vel = U;
    end

    t   = S.hwT_wind;
    vel = smooth_trace(vel, t, SMOOTH_T);

    plot(t, vel, 'LineWidth', 0.5, 'Color', colors(iRun,:), ...
        'DisplayName', sprintf('%s  hot-wire %s (%s)', runLabel, HW_COMPONENT, HW_PROBE));

    fprintf('%s: %d samples, t = %.1f..%.1f s, %s = %.3f..%.3f m/s (mean %.3f)\n', ...
        runLabel, numel(vel), min(t), max(t), HW_COMPONENT, ...
        min(vel), max(vel), mean(vel));
end

%% ---- Reference probe (T29 block of probe4.txt) --------------------------
refPath = fullfile(HW_DIR, REF_FILE);
if isfile(refPath)
    R = load(refPath);
    [~, refTag] = fileparts(REF_FILE);
    refLabel = erase(erase(refTag, 'hotwire_cal_'), '_partial');

    calRef = parse_probe_section(CAL_FILE, REF_PROBE);
    if calRef.nSensors ~= 1
        warning('plot_cal_velocity:refSensors', ...
            '%s block has %d sensors; using the first.', REF_PROBE, calRef.nSensors);
    end

    E = R.E_ref(:);
    U_ref = polyval(fliplr(calRef.C(1,:)), E);

    % Same out-of-range handling as convert_E2U_fn: below the min calibration
    % point the polynomial diverges, so a straight line from the origin through
    % that point is used instead. Those samples are extrapolation, not
    % calibration. Above the max point the fit is unconstrained, so it is
    % clamped and counted.
    below = E < calRef.E_floor(1);
    U_ref(below) = (calRef.U_floor(1) / calRef.E_floor(1)) * E(below);
    above = E > calRef.E_ceil(1);
    U_ref(above) = calRef.U_ceil(1);

    fprintf(['%s: %d samples, t = %.1f..%.1f s, E_ref = %.4f..%.4f V, ' ...
             'U_ref = %.3f..%.3f m/s (mean %.3f)\n'], ...
        refLabel, numel(U_ref), min(R.hwT_wind), max(R.hwT_wind), ...
        min(E), max(E), min(U_ref), max(U_ref), mean(U_ref));
    fprintf(['  %s calibrated range: %.3f..%.3f m/s (%.4f..%.4f V), ' ...
             'T_ref = %.2f degC\n'], REF_PROBE, calRef.U_floor(1), calRef.U_ceil(1), ...
        calRef.E_floor(1), calRef.E_ceil(1), calRef.T_ref);
    if any(below)
        fprintf('  %.1f%% of samples below the calibrated floor (extrapolated).\n', ...
            100*nnz(below)/numel(E));
    end
    if any(above)
        fprintf('  %.1f%% of samples above the calibrated ceiling (clamped).\n', ...
            100*nnz(above)/numel(E));
    end

    U_ref = smooth_trace(U_ref, R.hwT_wind, SMOOTH_T);
    plot(R.hwT_wind, U_ref, 'LineWidth', 0.5, 'Color', colors(3,:), ...
        'DisplayName', sprintf('%s  ref probe (%s)', refLabel, REF_PROBE));
else
    warning('plot_cal_velocity:noRef', 'Reference-probe file not found: %s', refPath);
end

hold off;
xlabel('Time since first fan write (s)');
ylabel('Velocity (m/s)');
title(sprintf('260909 calibration runs -- velocities from %s', ...
    'probe4.txt (55P95 hot-wire, T29 reference)'));
legend('Interpreter', 'none', 'Location', 'best');
grid on;

if SAVE_FIG
    print(gcf, fullfile(FIG_DIR, 'cal_velocity_260909.png'), '-dpng', '-r150');
end

%% ---- local functions -----------------------------------------------------
function y = smooth_trace(y, t, winT)
% Moving average over winT seconds; pass winT = 0 to leave the trace alone.
    if winT <= 0
        return;
    end
    fs = 1 / median(diff(t), 'omitnan');
    n  = max(1, round(winT * fs));
    y  = movmean(y, n);
end
