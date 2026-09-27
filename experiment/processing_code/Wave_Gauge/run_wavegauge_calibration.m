%% run_wavegauge_calibration.m
% Stand-alone static calibration of the wave gauge on ai7. The gauge is moved
% by hand to a series of known heights (still water, no fan); at each height
% the script records a short block on ai7, discards the first part (meniscus /
% hand disturbance settling) and averages the rest into one calibration
% point. At the end a polynomial  rel_height = f(voltage)  is fitted and saved
% in the SAME struct format WaveGauge_calibration.m produced, so
% plotWaveHeight_time.m can load the result unchanged:
%   calibration.polynomial_coefficients  -> polyval(p, voltage) gives cm
%
% This replaces typing voltages into WaveGauge_calibration.m by hand: the
% voltages are measured, and every point keeps its raw trace + std so a
% point taken while the surface was still moving is visible.
%
% Same structure as Run_exp/run_hotwire_calibration.m, but no fan is touched
% -- the tunnel must be OFF and the water still for a static calibration.
%
% Sign convention (from WaveGauge_calibration.m):
%   rel_height = -(height_reading - HEIGHT_ZERO)
% i.e. the height reading is the gauge-carriage scale, and moving the gauge
% DOWN into the water reads as a positive surface elevation. Flip SIGN_FLIP if
% your scale runs the other way.
%
% ASSUMPTIONS (check before running):
%   - Wave gauge output is on DEV_ID/ai7. The old WaveGauge_receiver.m used
%     "Dev1"/ai7 SingleEnded; on the current rig Dev4 is the USB-6451 that
%     carries the probes. Confirm which device the gauge is wired to in NI MAX.
%   - Amplifier output stays within +-10 V (the old cal spanned -10..7.3 V, so
%     the bottom end was at the rail -- points pinned at -10 V are flagged).

%connected to AI3 gain 4.50 and zero at 0.50

clear; clc;

DEV_ID      = "Dev4";          % <-- confirm in NI MAX
WG_CHAN     = "ai7";
WG_TERM     = "SingleEnded";   % as in WaveGauge_receiver.m
WG_RANGE    = [-10 10];
FS          = 1000;            % Hz

POINT_SETTLE_T = 3;    % s, recorded but not averaged (surface settling after moving the gauge)
POINT_AVG_T    = 5;   % s, averaged into the calibration point

HEIGHT_ZERO = 46.5;     % scale reading (cm) at still-water level -- prompted below read at the bottom line instead of where 0 line is
SIGN_FLIP   = false;    % rel_height = -(h - HEIGHT_ZERO), as in WaveGauge_calibration.m
POLY_ORDER  = 3;       % cubic, as in WaveGauge_calibration.m
RAIL_TOL    = 0.05;    % V, points within this of the range limits are flagged as clipped

%% --- Still-water reference --------------------------------------------------
r = strtrim(input(sprintf('Scale reading at still-water level (cm) [%g]: ', HEIGHT_ZERO), 's'));
if ~isempty(r)
    HEIGHT_ZERO = str2double(r);
    if isnan(HEIGHT_ZERO), error('Not a number: "%s"', r); end
end

%% --- Set up DAQ -------------------------------------------------------------
dq = daq("ni");
ch = addinput(dq, DEV_ID, WG_CHAN, "Voltage");
ch.TerminalConfig = WG_TERM;
ch.Range = WG_RANGE;
dq.Rate = FS;
fprintf('Wave gauge on %s/%s (%s, %g..%g V) at %d Hz\n', ...
    DEV_ID, WG_CHAN, WG_TERM, WG_RANGE(1), WG_RANGE(2), FS);

% Quick live check so a dead channel shows up before the first point.
v = read(dq, seconds(1));
fprintf('Current reading: %.4f V (std %.4f V over 1 s)\n\n', mean(v.Variables), std(v.Variables));

%% --- Point loop -------------------------------------------------------------
% Heights can be taken in any order; they are sorted before fitting.
fprintf('For each point: move the gauge, let the surface settle, then enter\n');
fprintf('the scale reading (cm).  Other inputs:  r = redo last point,  q = finish\n\n');

ptHeight = [];            % scale reading (cm)
ptMean   = [];            % mean voltage over the averaging window
ptStd    = [];
ptTime   = datetime.empty;
ptRaw    = {};            % full trace per point (settle + avg), volts
ptClip   = logical([]);

nSettle = round(POINT_SETTLE_T * FS);
recT    = POINT_SETTLE_T + POINT_AVG_T;

figure('Name','Wave gauge calibration (live)');
axLive = axes; grid(axLive, 'on');
xlabel(axLive, 'Voltage (V)'); ylabel(axLive, 'Scale reading (cm)');

% First point is the still-water point itself, so its voltage is measured
% (and plotted) rather than only drawn as a height line.
r = strtrim(input(sprintf(['Set the gauge at still-water level (%g cm), then press Enter ' ...
    'to record it (s = skip): '], HEIGHT_ZERO), 's'));
firstPt = ~strcmpi(r, 's');

while true
    if firstPt
        r = sprintf('%.17g', HEIGHT_ZERO);   % exact round-trip, so the == match below finds it
        firstPt = false;
    else
        r = strtrim(input(sprintf('[pt %d] height (cm) / r / q: ', numel(ptHeight)+1), 's'));
    end
    if strcmpi(r, 'q')
        break;
    elseif strcmpi(r, 'r')
        if isempty(ptHeight)
            disp('Nothing to redo.'); continue;
        end
        fprintf('Removed point %.2f cm (%.4f V).\n', ptHeight(end), ptMean(end));
        ptHeight(end) = []; ptMean(end) = []; ptStd(end) = [];
        ptTime(end) = []; ptRaw(end) = []; ptClip(end) = [];
    else
        h = str2double(r);
        if isnan(h)
            fprintf('Not understood: "%s"\n', r); continue;
        end
        fprintf('  recording %.0f s (%.0f s settle + %.0f s avg)...', recT, POINT_SETTLE_T, POINT_AVG_T);
        d = read(dq, seconds(recT));
        vv = d.Variables;
        va = vv(nSettle+1:end);
        clipped = any(va <= WG_RANGE(1) + RAIL_TOL | va >= WG_RANGE(2) - RAIL_TOL);

        ptHeight(end+1) = h;             %#ok<SAGROW>
        ptMean(end+1)   = mean(va);      %#ok<SAGROW>
        ptStd(end+1)    = std(va);       %#ok<SAGROW>
        ptTime(end+1)   = datetime('now'); %#ok<SAGROW>
        ptRaw{end+1}    = vv;            %#ok<SAGROW>
        ptClip(end+1)   = clipped;       %#ok<SAGROW>

        fprintf(' %.4f V  (std %.4f V)%s\n', ptMean(end), ptStd(end), ...
            ternary(clipped, '  ** AT RANGE LIMIT -- excluded from fit **', ''));
    end

    cla(axLive); hold(axLive, 'on');
    errorbar(axLive, ptMean(~ptClip), ptHeight(~ptClip), ptStd(~ptClip), 'horizontal', 'o');
    if any(ptClip)
        plot(axLive, ptMean(ptClip), ptHeight(ptClip), 'rx', 'MarkerSize', 10);
    end
    yline(axLive, HEIGHT_ZERO, 'k--', 'still water');
    iZero = find(ptHeight == HEIGHT_ZERO, 1, 'last');
    if ~isempty(iZero)
        plot(axLive, ptMean(iZero), HEIGHT_ZERO, 'kp', 'MarkerSize', 14, 'MarkerFaceColor', 'y');
        xline(axLive, ptMean(iZero), 'k:', sprintf('%.3f V', ptMean(iZero)));
    end
    hold(axLive, 'off'); drawnow;
end

if isempty(ptHeight)
    disp('No points recorded -- nothing to fit.');
    return;
end

%% --- Fit --------------------------------------------------------------------
[height_calib, idx] = sort(ptHeight);
voltage_calib = ptMean(idx);
voltage_std   = ptStd(idx);
clip_calib    = ptClip(idx);
raw_calib     = ptRaw(idx);
time_calib    = ptTime(idx);

if SIGN_FLIP
    rel_height_calib = -(height_calib - HEIGHT_ZERO);
else
    rel_height_calib = height_calib - HEIGHT_ZERO;
end

use = ~clip_calib;
p = []; residuals = nan(size(rel_height_calib));
if nnz(use) <= POLY_ORDER
    warning('Only %d usable points for an order-%d fit -- add more heights.', nnz(use), POLY_ORDER);
else
    p = polyfit(voltage_calib(use), rel_height_calib(use), POLY_ORDER);
    residuals(use) = rel_height_calib(use) - polyval(p, voltage_calib(use));
    fprintf('\nOrder-%d fit, rel_height(cm) = polyval(p, V):\n  p = [%s]\n', ...
        POLY_ORDER, strjoin(compose('%.6g', p), ' '));
    fprintf('  Max residual: %.3f cm   RMS residual: %.3f cm   (%d points)\n', ...
        max(abs(residuals(use))), rms(residuals(use)), nnz(use));
end

ptTable = table(height_calib(:), rel_height_calib(:), voltage_calib(:), voltage_std(:), ...
    clip_calib(:), residuals(:), time_calib(:), ...
    'VariableNames', {'Height_cm','RelHeight_cm','Voltage_V','VoltageStd_V', ...
                      'Clipped','Residual_cm','Time'});
disp(' '); disp(ptTable);

%% --- Plots ------------------------------------------------------------------
figure('Name','Wave gauge calibration');
subplot(2,1,1); hold on;
errorbar(voltage_calib(use), rel_height_calib(use), voltage_std(use), 'horizontal', 'o', ...
    'DisplayName', 'data');
if any(clip_calib)
    plot(voltage_calib(clip_calib), rel_height_calib(clip_calib), 'rx', ...
        'MarkerSize', 10, 'DisplayName', 'clipped (not fitted)');
end
if ~isempty(p)
    vf = linspace(min(voltage_calib(use)), max(voltage_calib(use)), 200);
    plot(vf, polyval(p, vf), 'r-', 'LineWidth', 1.5, 'DisplayName', sprintf('order-%d fit', POLY_ORDER));
end
hold off; grid on; legend('Location','best');
xlabel('Voltage (V)'); ylabel('Relative height (cm)');
title(sprintf('Wave gauge %s/%s, zero = %g cm', DEV_ID, WG_CHAN, HEIGHT_ZERO));

subplot(2,1,2);
plot(voltage_calib, residuals, 'o-'); grid on; yline(0, 'k:');
xlabel('Voltage (V)'); ylabel('Residual (cm)');

% Raw trace of every point with the averaged part shaded -- check that
% POINT_SETTLE_T was long enough (no drift inside the averaging window).
figure('Name','Wave gauge calibration -- raw points');
nP = numel(raw_calib);
nc = ceil(sqrt(nP)); nr = ceil(nP / nc);
for k = 1:nP
    subplot(nr, nc, k);
    tk = (0:numel(raw_calib{k})-1) / FS;
    plot(tk, raw_calib{k}); hold on;
    xline(POINT_SETTLE_T, 'k--'); hold off; grid on;
    title(sprintf('%.1f cm', height_calib(k)));
    if k > nP - nc, xlabel('t (s)'); end
end

%% --- Save (prompted) --------------------------------------------------------
% Same fields as WaveGauge_calibration.m, plus the measured extras.
calibration.date = datetime('now');
calibration.height_zero = HEIGHT_ZERO;
calibration.height_calib = height_calib;
calibration.voltage_calib = voltage_calib;
calibration.rel_height_calib = rel_height_calib;
calibration.polynomial_coefficients = p;
calibration.polynomial_order = POLY_ORDER;
calibration.max_residual_cm = max(abs(residuals(use)));
calibration.rms_residual_cm = rms(residuals(use));
calibration.voltage_std = voltage_std;
calibration.clipped = clip_calib;
calibration.raw_voltage = raw_calib;
calibration.point_time = time_calib;
calibration.settle_t = POINT_SETTLE_T;
calibration.avg_t = POINT_AVG_T;
calibration.fs = FS;
calibration.device = DEV_ID;
calibration.channel = WG_CHAN;
calibration.terminal = WG_TERM;
calibration.sign_flip = SIGN_FLIP;

% Saved next to this script (not MATLAB's current folder), which is where
% run_wavegauge_fan_calibration.m looks for it.
saveName = fullfile(fileparts(mfilename('fullpath')), ...
    sprintf('wave_gauge_calibration_%s.mat', datestr(calibration.date, 'yyyymmdd_HHMM')));
drawnow;
resp = input(sprintf('\nSave this calibration to %s? [Y/n]: ', saveName), 's');
if isempty(resp) || strncmpi(strtrim(resp), 'y', 1)
    save(saveName, 'calibration');
    csvName = strrep(saveName, '.mat', '_points.csv');
    writetable(ptTable, csvName);
    fprintf('Saved %s and %s\n', saveName, csvName);
else
    disp('Not saved. To save later, run:');
    fprintf('  save(saveName, ''calibration'')\n');
end


%% Local functions
function out = ternary(cond, a, b)
    if cond, out = a; else, out = b; end
end
