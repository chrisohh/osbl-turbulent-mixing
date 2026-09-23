%% analyze_hotwire_cal_steps.m
% Same style as analyze_steps_from_trace.m (settle-time sweep, all-holds-
% overlaid settling plot, velocity-vs-fan-voltage), applied to a saved
% hot-wire calibration run instead of a recovered ref-probe trace.
%
% Unlike analyze_steps_from_trace.m this does not need to rebuild the step
% schedule -- hotwire_cal_*.mat already has stepV/stepStart/stepEnd/stepDir
% from run_hotwire_calibration.m -- so this script re-averages each hold at
% a chosen SETTLE_T (independent of whatever STEP_SETTLE_T was used live)
% and plots the settling behaviour directly.
%
% Velocity is U from convert_E2U_fn using CAL_FILE below (probe4.txt, 55P95
% block), not the U saved in the .mat (which used the acquisition-time
% probe3.txt).

clear; clc;

%% ---- Config -------------------------------------------------------------
MAT_FILE = 'D:\HLAB_2026\hotwire\hotwire_cal_20260909_212437.mat';
CAL_FILE = 'C:\Users\airsealab\Documents\GitHub\osbl-turbulent-mixing\experiment\data\260909\probe4.txt';
HW_PROBE = '55P95';

SETTLE_T = 20;                       % s discarded at the start of each hold before averaging
SWEEP    = [0 5 10 15 20 30 45 60];  % settle times compared in the sweep table
MIN_DWELL = 5;                       % s of settled data required to keep a truncated hold

CTA_DIR = fileparts(mfilename('fullpath'));
addpath(CTA_DIR);
% -------------------------------------------------------------------------

firstName = regexp(fileread(CAL_FILE), 'Probe name:\s*([^\r\n]*)', 'tokens', 'once');
if ~strcmpi(strtrim(firstName{1}), HW_PROBE)
    error('analyze_hotwire_cal_steps:wrongFirstProbe', ...
        'First probe block in %s is "%s", not %s.', CAL_FILE, strtrim(firstName{1}), HW_PROBE);
end

S = load(MAT_FILE);
[~, tag] = fileparts(MAT_FILE);
fprintf('Loaded %s (%d samples)\n', tag, numel(S.hwT_wind));
fprintf('Velocity from %s (%s block)\n', CAL_FILE, HW_PROBE);

%% 1. Velocity trace (probe4/55P95, NOT the U saved with the acquisition cal)
[U, ~, ~] = convert_E2U_fn(S.E1, S.E2, S.E3, CAL_FILE);
t = S.hwT_wind(:);
U = U(:);

%% 2. Step schedule -- taken directly from the saved run ---------------------
stepV = S.stepV(:)';
sS    = S.stepStart(:)';
sE    = S.stepEnd(:)';

sE = min(sE, t(end));
usable = (sE - sS - SETTLE_T) >= MIN_DWELL;
truncated = sE < S.stepEnd(:)' & usable;

if any(truncated)
    fprintf('Keeping %d truncated hold(s) (levels %s V) -- %s s of settled data left.\n', ...
        sum(truncated), strjoin(compose('%.2f', stepV(truncated)), ' '), ...
        strjoin(compose('%.0f', sE(truncated)-sS(truncated)-SETTLE_T), ' '));
end
if any(~usable)
    fprintf('Dropping %d hold(s) with <%d s of settled data (levels %s V).\n', ...
        sum(~usable), MIN_DWELL, strjoin(compose('%.2f', stepV(~usable)), ' '));
end
stepV = stepV(usable); sS = sS(usable); sE = sE(usable);
stepDir = S.stepDir(usable);
n = numel(stepV);

%% 3. Settle sweep -----------------------------------------------------------
M = nan(n, numel(SWEEP));  Sd = nan(n, numel(SWEEP));
for j = 1:numel(SWEEP)
    for k = 1:n
        m = t >= sS(k)+SWEEP(j) & t < sE(k);
        if sum(m) > 10, M(k,j) = mean(U(m)); Sd(k,j) = std(U(m)); end
    end
end

fprintf('\nMean U (m/s) by settle time discarded:\n%6s', 'FanV');
for j = 1:numel(SWEEP)
    fprintf('%9s', sprintf('%ds', SWEEP(j)));
end
fprintf('%10s\n', 'drift%');
for k = 1:n
    drift = 100*(M(k,end)-M(k,1))/M(k,1);
    fprintf('%6.2f', stepV(k)); fprintf('%9.4f', M(k,:)); fprintf('%10.2f\n', drift);
end
fprintf('\nSuggested settle: the column where drift stops changing materially.\n');

%% 4. Averages at the chosen SETTLE_T ---------------------------------------
Um = nan(n,1); Us = nan(n,1); Un = zeros(n,1);
for k = 1:n
    m = t >= sS(k)+SETTLE_T & t < sE(k);
    Un(k) = sum(m);
    if Un(k) > 10, Um(k) = mean(U(m)); Us(k) = std(U(m)); end
end
T = table(stepV(:), string(stepDir(:)), sS(:), sE(:), Un, Um, Us, Us./Um, ...
    'VariableNames', {'FanV','Direction','tStart','tEnd','nSamples','U_mean','U_std','TurbIntensity'});
fprintf('\nSettled averages (discarding %d s per hold):\n', SETTLE_T); disp(T);

%% 5. Plots -------------------------------------------------------------
fig = figure('Position',[80 80 1150 850],'Color','w', 'Name', ...
    sprintf('%s -- settling analysis (%s)', tag, HW_PROBE));
tiledlayout(2,2,'TileSpacing','compact','Padding','compact');

ax = nexttile([1 2]);
plot(t, U, 'LineWidth', 0.5, 'Color', [0.0 0.45 0.74]); hold on;
shade_windows(ax, sS+SETTLE_T, sE);
hold off;
xlabel('Time since first fan write (s)'); ylabel('U (m/s)');
title(sprintf('%s -- shaded = averaged (settle %d s discarded), cal = %s', ...
    strrep(tag,'_','\_'), SETTLE_T, HW_PROBE));
grid on; xlim([t(1) t(end)]);

nexttile; hold on;
% Every hold on a common "time since step change" axis -- the settling time
% read directly off the data rather than guessed.
cols = parula(n);
for k = 1:n
    m = t >= sS(k) & t < sE(k);
    plot(t(m)-sS(k), U(m), 'Color', [cols(k,:) 0.85], 'LineWidth', 0.6, ...
        'DisplayName', sprintf('%.1f V', stepV(k)));
end
yl = ylim; plot([SETTLE_T SETTLE_T], yl, 'r--', 'LineWidth', 1.5, 'DisplayName','settle cut');
ylim(yl); hold off;
xlabel('Time since step change (s)'); ylabel('U (m/s)');
title('Settling, all holds overlaid'); legend('Location','eastoutside','FontSize',7); grid on;

nexttile;
errorbar(stepV, Um, Us, 'o-', 'LineWidth', 1.6, 'MarkerSize', 7, ...
    'MarkerFaceColor', 'w', 'CapSize', 5);
xlabel('Fan command (V)'); ylabel('U (m/s)');
title(sprintf('Velocity vs fan voltage (settled, %d s discarded)', SETTLE_T));
grid on;

fprintf('\nDone. Settled points are in table T.\n');

%% ---- local functions -----------------------------------------------------
function shade_windows(ax, tStart, tEnd)
% Grey band behind each averaging window. Tries XREGION (R2023b+, draws
% behind the data automatically); falls back to a plain patch on older
% MATLAB, which does not reorder z-order but is visible via transparency.
    ok = ~isnan(tStart) & ~isnan(tEnd);
    if ~any(ok), return; end
    tStart = tStart(ok); tEnd = tEnd(ok);

    wasHeld = ishold(ax); hold(ax, 'on');
    try
        h = xregion(ax, tStart(:), tEnd(:), ...
            'FaceColor', [0.1 0.3 0.7], 'FaceAlpha', 0.12);
        set(h, 'HandleVisibility', 'off');
    catch
        yl = ylim(ax);
        for k = 1:numel(tStart)
            patch(ax, [tStart(k) tEnd(k) tEnd(k) tStart(k)], [yl(1) yl(1) yl(2) yl(2)], ...
                [0.1 0.3 0.7], 'EdgeColor', 'none', 'FaceAlpha', 0.12, ...
                'HandleVisibility', 'off');
        end
        ylim(ax, yl);
    end
    if ~wasHeld, hold(ax, 'off'); end
end
