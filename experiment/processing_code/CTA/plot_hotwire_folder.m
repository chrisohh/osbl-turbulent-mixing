%% plot_hotwire_folder.m
% Plot every hotwire_*.mat run found in a data folder, re-deriving the
% velocities U,V,W from the raw sensor voltages E1,E2,E3 with the probe
% calibration set below (probe4) instead of the CAL_FILE that was baked in
% at acquisition time (probe3 for the 260909 runs).
%
% The saved U,V,W in the .mat files are NOT used for the main plots -- they
% are only kept for the probe3-vs-probe4 comparison figure.
%
% Files that logged only the reference probe (E_ref, no hot-wire channels)
% are plotted as reference velocity via convert_Eref2Uref.m.
%
% Usage: just run this script (edit the config block first).

clear; clc;

%% ---- Config -------------------------------------------------------------
HW_DIR   = 'D:\HLAB_2026\hotwire';                    % folder holding hotwire_*.mat
CAL_FILE = 'C:\Users\airsealab\Documents\GitHub\osbl-turbulent-mixing\experiment\data\260909\probe4.txt';

COMPARE_ORIGINAL = true;   % overlay the U,V,W saved with the acquisition-time cal
SAVE_FIGS        = false;  % write PNGs into FIG_DIR
FIG_DIR          = fullfile(HW_DIR, 'figs');
% -------------------------------------------------------------------------

CTA_DIR = fileparts(mfilename('fullpath'));   % parse_calibration + convert_E2U_fn live here
addpath(CTA_DIR);

if ~isfile(CAL_FILE)
    error('plot_hotwire_folder:noCal', 'Calibration file not found: %s', CAL_FILE);
end
if SAVE_FIGS && ~isfolder(FIG_DIR)
    mkdir(FIG_DIR);
end

files = dir(fullfile(HW_DIR, 'hotwire_*.mat'));
files = files(~contains({files.name}, 'hotwire_cal_'));   % skip calibration sweeps
if isempty(files)
    error('plot_hotwire_folder:noFiles', 'No hotwire_*.mat files in %s', HW_DIR);
end

[~, calName] = fileparts(CAL_FILE);
fprintf('Calibration in use: %s (%s)\n', calName, CAL_FILE);
fprintf('Found %d run file(s) in %s\n\n', numel(files), HW_DIR);

for iFile = 1:numel(files)
    matFile = fullfile(files(iFile).folder, files(iFile).name);
    [~, tag] = fileparts(files(iFile).name);
    S = load(matFile);

    fprintf('=== %s ===\n', files(iFile).name);
    if isfield(S, 'CAL_FILE')
        fprintf('Acquisition-time calibration: %s\n', S.CAL_FILE);
    end

    hasHotwire = isfield(S, 'E1') && ~isempty(S.E1);
    hasRef     = isfield(S, 'E_ref') && ~isempty(S.E_ref);
    hasFan     = isfield(S, 'sentT') && ~isempty(S.sentT);

    % Time base: 0 = wind start
    if isfield(S, 'hwT_wind') && ~isempty(S.hwT_wind)
        t = S.hwT_wind;
    elseif isfield(S, 'hwT')
        t = S.hwT;
    else
        warning('No time vector in %s -- skipped.', files(iFile).name);
        continue;
    end
    onIdx = t >= 0;

    fanT = [];
    if hasFan
        fanT = S.sentT;
        if isfield(S, 'fanStartElapsed')
            fanT = fanT - S.fanStartElapsed;
        end
    end

    %% ---- Hot-wire run: re-derive velocities with the chosen calibration --
    if hasHotwire
        [U, V, W, Uc1, Uc2, Uc3] = convert_E2U_fn(S.E1, S.E2, S.E3, CAL_FILE);
        speed = sqrt(U.^2 + V.^2 + W.^2);

        figure('Name', sprintf('%s -- %s', tag, calName), ...
               'Position', [60 60 950 800]);

        subplot(3,1,1);
        if hasFan
            plot(fanT, S.sentV, 'LineWidth', 1.5, 'Color', [0.85 0.33 0.10]);
        else
            text(0.5, 0.5, 'no fan record', 'Units', 'normalized', ...
                 'HorizontalAlignment', 'center');
        end
        ylabel('Fan (V)');
        title(sprintf('%s  --  velocities from %s', strrep(tag, '_', '\_'), calName));
        grid on;

        subplot(3,1,2);
        plot(t, S.E1, t, S.E2, t, S.E3, 'LineWidth', 0.5);
        ylabel('Sensor (V)'); legend('E1','E2','E3', 'Location', 'best'); grid on;

        subplot(3,1,3);
        plot(t, U, t, V, t, W, 'LineWidth', 0.5);
        xlabel('Time since wind start (s)'); ylabel('Velocity (m/s)');
        legend('U','V','W', 'Location', 'best'); grid on;

        linkaxes(findobj(gcf, 'Type', 'axes'), 'x');
        if SAVE_FIGS
            saveas(gcf, fullfile(FIG_DIR, sprintf('%s_%s.png', tag, calName)));
        end

        % Comparison with the velocities saved at acquisition time
        if COMPARE_ORIGINAL && isfield(S, 'U') && ~isempty(S.U)
            [~, oldName] = fileparts(S.CAL_FILE);
            figure('Name', sprintf('%s -- %s vs %s', tag, oldName, calName), ...
                   'Position', [1020 60 900 700]);
            comps   = {'U','V','W'};
            newVals = {U, V, W};
            oldVals = {S.U, S.V, S.W};
            for k = 1:3
                subplot(3,1,k);
                plot(t, oldVals{k}, 'LineWidth', 0.5); hold on;
                plot(t, newVals{k}, 'LineWidth', 0.5); hold off;
                ylabel([comps{k} ' (m/s)']); grid on;
                if k == 1
                    legend(oldName, calName, 'Interpreter', 'none', 'Location', 'best');
                    title(sprintf('%s  --  acquisition cal vs %s', ...
                        strrep(tag, '_', '\_'), calName));
                end
            end
            xlabel('Time since wind start (s)');
            linkaxes(findobj(gcf, 'Type', 'axes'), 'x');
            if SAVE_FIGS
                saveas(gcf, fullfile(FIG_DIR, ...
                    sprintf('%s_%s_vs_%s.png', tag, oldName, calName)));
            end
        end

        %% Statistics over the wind-on portion
        fprintf('\n-- Wind-on statistics (t >= 0, %d samples), cal = %s --\n', ...
            nnz(onIdx), calName);
        fprintf('Mean:  U=%7.3f  V=%7.3f  W=%7.3f  |speed|=%7.3f m/s\n', ...
            mean(U(onIdx)), mean(V(onIdx)), mean(W(onIdx)), mean(speed(onIdx)));
        fprintf('RMS:   U=%7.3f  V=%7.3f  W=%7.3f m/s\n', ...
            std(U(onIdx)), std(V(onIdx)), std(W(onIdx)));
        if mean(U(onIdx)) ~= 0
            fprintf('Tu:    u=%5.1f%%  v=%5.1f%%  w=%5.1f%%  (normalised by mean U)\n', ...
                100*std(U(onIdx))/mean(U(onIdx)), ...
                100*std(V(onIdx))/mean(U(onIdx)), ...
                100*std(W(onIdx))/mean(U(onIdx)));
        end
        if COMPARE_ORIGINAL && isfield(S, 'U') && ~isempty(S.U)
            fprintf('Acquisition cal mean: U=%7.3f  V=%7.3f  W=%7.3f m/s\n', ...
                mean(S.U(onIdx)), mean(S.V(onIdx)), mean(S.W(onIdx)));
        end

        % How much of the record sits below the calibrated range -- there the
        % polynomial is replaced by the origin-to-floor line, so those samples
        % are an extrapolation, not a calibration.
        cal   = parse_calibration(CAL_FILE);
        below = [nnz(S.E1 < cal.E_floor(1)), nnz(S.E2 < cal.E_floor(2)), ...
                 nnz(S.E3 < cal.E_floor(3))] / numel(S.E1) * 100;
        fprintf('Samples below calibrated range: E1 %.1f%%, E2 %.1f%%, E3 %.1f%%\n', below);
        fprintf('Per-sensor mean calibration velocity: %.3f / %.3f / %.3f m/s\n', ...
            mean(Uc1(onIdx)), mean(Uc2(onIdx)), mean(Uc3(onIdx)));
    end

    %% ---- Reference-probe-only run ---------------------------------------
    if hasRef
        [U_ref, cal_ref] = convert_Eref2Uref(S.E_ref);
        fprintf('\nReference probe: %s s/n %s, certificate %s (%s)\n', ...
            cal_ref.probe, cal_ref.serial, cal_ref.id, cal_ref.date);

        figure('Name', sprintf('%s -- reference probe', tag), ...
               'Position', [60 60 950 700]);
        subplot(3,1,1);
        if hasFan
            plot(fanT, S.sentV, 'LineWidth', 1.5, 'Color', [0.85 0.33 0.10]);
        end
        ylabel('Fan (V)');
        title(sprintf('%s  --  reference probe (54T29)', strrep(tag, '_', '\_')));
        grid on;

        subplot(3,1,2);
        plot(t, S.E_ref, 'LineWidth', 0.5); ylabel('E_{ref} (V)'); grid on;

        subplot(3,1,3);
        plot(t, U_ref, 'LineWidth', 0.5, 'Color', [0.49 0.18 0.56]);
        xlabel('Time since wind start (s)'); ylabel('U_{ref} (m/s)'); grid on;

        linkaxes(findobj(gcf, 'Type', 'axes'), 'x');
        if SAVE_FIGS
            saveas(gcf, fullfile(FIG_DIR, sprintf('%s_refprobe.png', tag)));
        end

        fprintf('Mean U_ref = %.3f m/s, RMS = %.3f m/s (t >= 0)\n', ...
            mean(U_ref(onIdx)), std(U_ref(onIdx)));
    end

    if ~hasHotwire && ~hasRef
        fprintf('No hot-wire or reference channels in this file -- nothing to plot.\n');
    end
    fprintf('\n');
end
