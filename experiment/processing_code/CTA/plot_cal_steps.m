%% plot_cal_steps.m
% Step-averaged calibration points from run_hotwire_calibration.m runs.
%
% Figure 1 -- hot-wire: the two 260909 step runs (211425 and 212437) drawn on
%             the same axes, one subplot per sensor, up-leg and down-leg kept
%             separate so hysteresis stays visible.
% Figure 2 -- reference probe (Dantec 54T29, the probe named in probe4.txt):
%             the 260909 17:20 step run, fan command vs E_ref and U_ref.
%
% Each point is the mean over that step's averaging window
% ([stepStart+settle, stepEnd]); error bars are the within-step std.
%
% Usage: just run this script (edit the config block first).

clear; clc;

%% ---- Config -------------------------------------------------------------
HW_DIR = 'D:\HLAB_2026\hotwire';

HW_FILES = { 'hotwire_cal_20260909_211425.mat', ...
             'hotwire_cal_20260909_212437.mat' };

% The only reference-probe step run from 260909 ~17:00 in this folder. It is
% named hotwire_cal_* (the calibration script's naming) but logged RefProbe
% only -- no hot-wire channels.
REF_FILE = 'hotwire_cal_20260909_172006_partial.mat';

SAVE_FIGS = false;
FIG_DIR   = fullfile(HW_DIR, 'figs');
% -------------------------------------------------------------------------

addpath(fileparts(mfilename('fullpath')));
if SAVE_FIGS && ~isfolder(FIG_DIR)
    mkdir(FIG_DIR);
end

runColors = [0.000 0.447 0.741;
             0.850 0.325 0.098;
             0.466 0.674 0.188];

%% ===== Figure 1: hot-wire step points, both runs overlaid =================
figure('Name', 'Hot-wire calibration steps', 'Position', [60 60 1000 800]);
sensorNames = {'Probe1','Probe2','Probe3'};
legEntries  = gobjects(0);
legLabels   = {};

for iRun = 1:numel(HW_FILES)
    matFile = fullfile(HW_DIR, HW_FILES{iRun});
    if ~isfile(matFile)
        warning('plot_cal_steps:missing', 'Not found: %s -- skipped.', matFile);
        continue;
    end
    S = load(matFile);
    [~, tag] = fileparts(HW_FILES{iRun});
    runLabel = erase(tag, 'hotwire_cal_');
    col = runColors(mod(iRun-1, size(runColors,1)) + 1, :);

    isUp   = strcmp(S.stepDir(:), 'up')   & S.stepN > 0;
    isDown = strcmp(S.stepDir(:), 'down') & S.stepN > 0;
    isZero = strcmp(S.stepDir(:), 'zero') & S.stepN > 0;

    fprintf('=== %s ===\n', HW_FILES{iRun});
    fprintf('Fan levels: %s V | %d steps (%d up, %d down, %d zero) | %d Hz\n', ...
        mat2str(S.FAN_V_LEVELS), numel(S.stepV), nnz(isUp), nnz(isDown), ...
        nnz(isZero), S.HOTWIRE_FS);

    for i = 1:3
        subplot(3,1,i); hold on;
        c = strcmp(S.aiNames, sensorNames{i});

        h = errorbar(S.stepV(isUp), S.stepMean(isUp,c), S.stepStd(isUp,c), ...
            'o-', 'Color', col, 'MarkerFaceColor', col, 'LineWidth', 1.2, ...
            'MarkerSize', 5, 'CapSize', 4);
        if any(isDown)
            errorbar(S.stepV(isDown), S.stepMean(isDown,c), S.stepStd(isDown,c), ...
                's--', 'Color', col, 'LineWidth', 1.2, 'MarkerSize', 6, 'CapSize', 4);
        end
        if any(isZero)
            plot(S.stepV(isZero), S.stepMean(isZero,c), 'x', ...
                'Color', col, 'MarkerSize', 9, 'LineWidth', 1.2);
        end

        if i == 1
            legEntries(end+1) = h; %#ok<SAGROW>
            legLabels{end+1}  = runLabel; %#ok<SAGROW>
        end
        hold off;
    end

    % Step table to the console so the numbers behind the points are visible
    fprintf('Fan V   :  %s\n', num2str(S.stepV, '%8.2f'));
    for i = 1:3
        c = strcmp(S.aiNames, sensorNames{i});
        fprintf('mean E%d :  %s\n', i, num2str(S.stepMean(:,c)', '%8.4f'));
    end
    fprintf('\n');
end

for i = 1:3
    subplot(3,1,i);
    ylabel(sprintf('E_%d (V)', i)); grid on;
    if i == 1
        title('Hot-wire step averages -- solid/circle = up-leg, dashed/square = down-leg');
        if ~isempty(legEntries)
            legend(legEntries, legLabels, 'Interpreter', 'none', 'Location', 'best');
        end
    end
    if i == 3
        xlabel('Fan command (V)');
    end
end

if SAVE_FIGS
    print(gcf, fullfile(FIG_DIR, 'cal_steps_hotwire_211425_vs_212437.png'), '-dpng', '-r150');
end

%% ===== Figure 2: reference-probe step points ==============================
refPath = fullfile(HW_DIR, REF_FILE);
if ~isfile(refPath)
    warning('plot_cal_steps:noRef', 'Reference-probe file not found: %s', refPath);
    return;
end

R = load(refPath);
[~, refTag] = fileparts(REF_FILE);
fprintf('=== %s (reference probe) ===\n', REF_FILE);
if isfield(R, 'cal_ref') && isstruct(R.cal_ref)
    fprintf('Probe %s s/n %s, certificate %s (%s)\n', ...
        R.cal_ref.probe, R.cal_ref.serial, R.cal_ref.id, R.cal_ref.date);
end
fprintf('Fan levels: %s V | %d steps\n', mat2str(R.FAN_V_LEVELS), numel(R.stepV));

isUp   = strcmp(R.stepDir(:), 'up')   & R.stepN > 0;
isDown = strcmp(R.stepDir(:), 'down') & R.stepN > 0;
isZero = strcmp(R.stepDir(:), 'zero') & R.stepN > 0;
cRef   = strcmp(R.aiNames, 'RefProbe');

nValid = nnz(R.stepN > 0);
if nValid < numel(R.stepV)
    warning('plot_cal_steps:emptySteps', ...
        ['%d of %d step windows contain no samples in %s -- the logged record ' ...
         'covers t = %.0f..%.0f s while the steps ran t = %.0f..%.0f s.'], ...
        numel(R.stepV) - nValid, numel(R.stepV), REF_FILE, ...
        min(R.hwT_wind), max(R.hwT_wind), min(R.stepStart), max(R.stepEnd));
end

figure('Name', 'Reference probe calibration steps', 'Position', [1080 60 900 800]);

% Raw trace with the step-averaging windows marked, so a mismatch between the
% planned schedule and the record that actually survived is visible.
subplot(3,1,1);
yyaxis left;  plot(R.hwT_wind, R.E_ref, 'LineWidth', 0.6); ylabel('E_{ref} (V)');
yyaxis right; plot(R.hwT_wind, R.U_ref, 'LineWidth', 0.6); ylabel('U_{ref} (m/s)');
shade_windows(R.avgStart, R.avgEnd);
xlim([min([R.hwT_wind(:); R.stepStart(:)]), max([R.hwT_wind(:); R.stepEnd(:)])]);
xlabel('Time since first fan write (s)'); grid on;
title(sprintf('%s  --  reference probe (54T29), shaded = step windows', ...
    strrep(refTag, '_', '\_')));

subplot(3,1,2); hold on;
errorbar(R.stepV(isUp), R.stepMean(isUp,cRef), R.stepStd(isUp,cRef), ...
    'o-', 'Color', runColors(1,:), 'MarkerFaceColor', runColors(1,:), ...
    'LineWidth', 1.2, 'MarkerSize', 5, 'CapSize', 4, 'DisplayName', 'up');
if any(isDown)
    errorbar(R.stepV(isDown), R.stepMean(isDown,cRef), R.stepStd(isDown,cRef), ...
        's--', 'Color', runColors(1,:), 'LineWidth', 1.2, 'MarkerSize', 6, ...
        'CapSize', 4, 'DisplayName', 'down');
end
if any(isZero)
    plot(R.stepV(isZero), R.stepMean(isZero,cRef), 'kx', 'MarkerSize', 9, ...
        'LineWidth', 1.2, 'DisplayName', 'zero');
end
hold off;
ylabel('mean E_{ref} (V)'); grid on; legend('Location', 'best');
title('Step averages -- fan command vs reference voltage');

subplot(3,1,3); hold on;
errorbar(R.stepV(isUp), R.stepUrefMean(isUp), R.stepUrefStd(isUp), ...
    'o-', 'Color', runColors(2,:), 'MarkerFaceColor', runColors(2,:), ...
    'LineWidth', 1.2, 'MarkerSize', 5, 'CapSize', 4, 'DisplayName', 'up');
if any(isDown)
    errorbar(R.stepV(isDown), R.stepUrefMean(isDown), R.stepUrefStd(isDown), ...
        's--', 'Color', runColors(2,:), 'LineWidth', 1.2, 'MarkerSize', 6, ...
        'CapSize', 4, 'DisplayName', 'down');
end
if any(isZero)
    plot(R.stepV(isZero), R.stepUrefMean(isZero), 'kx', 'MarkerSize', 9, ...
        'LineWidth', 1.2, 'DisplayName', 'zero');
end
hold off;
xlabel('Fan command (V)'); ylabel('U_{ref} (m/s)'); grid on; legend('Location', 'best');
title('Fan command -> reference velocity');

fprintf('Fan V    :  %s\n', num2str(R.stepV, '%8.2f'));
fprintf('mean Eref:  %s\n', num2str(R.stepMean(:,cRef)', '%8.4f'));
fprintf('mean Uref:  %s\n', num2str(R.stepUrefMean', '%8.3f'));
fprintf('\n');

if SAVE_FIGS
    print(gcf, fullfile(FIG_DIR, 'cal_steps_refprobe_172006.png'), '-dpng', '-r150');
end

%% ---- local functions -----------------------------------------------------
function shade_windows(tStart, tEnd)
% Shade the per-step averaging windows behind whatever is already plotted.
    ax = gca;
    yl = ylim(ax);
    for k = 1:numel(tStart)
        p = patch(ax, [tStart(k) tEnd(k) tEnd(k) tStart(k)], ...
                      [yl(1) yl(1) yl(2) yl(2)], [0.9 0.9 0.9], ...
                      'EdgeColor', 'none', 'FaceAlpha', 0.25, ...
                      'HandleVisibility', 'off');
        % uistack refuses to reorder children on a yyaxis (two-ruler) axes, so
        % the patch simply stays on top there -- hence the low FaceAlpha.
        try %#ok<TRYNC>
            uistack(p, 'bottom');
        end
    end
    ylim(ax, yl);
end
