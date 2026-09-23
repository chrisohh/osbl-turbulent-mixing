%% run_fan_kick_calibration.m
% Copy of run_hotwire_calibration.m with a KICK-START: after the zero hold the
% fan gets KICK_V for KICK_T seconds (enough to break it free), then drops to
% the staircase levels -- which may sit BELOW the 2 V start-up threshold,
% since a spinning fan keeps turning at a lower voltage than it needs to
% start. The aim is a smaller 0 -> first-velocity jump at wind start.
%
% MODE picks the level list (edit the presets below):
%   'breakaway' -- no kick, slow climb from below 2 V: where does it START?
%   'stall'     -- kick, then step DOWN: where does it STOP once spinning?
%   'cal'       -- kick, drop to the lowest keep-spinning level, step UP:
%                  the V -> U transfer function under the same start-up as
%                  the experiment (feed stepTable to fit_fan_transfer.m).
% In 'stall' the levels are walked in the order listed, but still labelled
% 'up' in stepTable -- the label means "part of the staircase", not direction.
%
% Kick hold is never averaged (its settle is its whole length).
%
% ---- original header (run_hotwire_calibration.m) ----
% Stand-alone in-situ calibration run: drives the fan through a STAIRCASE of
% voltage levels (up then back down) while logging the hot-wire and/or the
% reference probe continuously. Each level is held long enough that the flow
% settles, and the settled part of each hold is averaged into one calibration
% point -- giving an E vs U_ref curve for each hot-wire sensor.
%
% This is deliberately SEPARATE from run_experiment.m:
%   - no cameras, no camera-trigger DAQ session, no DELAY_BEFORE_TRIG
%   - fan is a staircase with holds, not a continuous linear ramp
%   - the output of interest is a per-step averaged table (+ polynomial fit),
%     not a time series of the mixing event
% Mixing the two into one script would mean a mode flag threaded through
% every timing branch, so they stay apart.
%
% ASSUMPTIONS (same hardware as run_experiment.m -- check before running):
%   - USB-6451 (fan + probes) is "Dev4" in NI MAX
%   - Fan on Dev4/ao0, hot-wire on Dev4/ai0:ai2, reference probe on Dev4/ai4
%
% NOT verified: how long the flow actually takes to settle after a fan step.
% STEP_SETTLE_T below is a guess -- inspect the per-step time series the
% script plots (settle window is shaded) and increase it if the traces are
% still drifting when averaging begins.

clear; clc;

DEV_ID = "Dev4";   % <-- confirm this matches NI MAX for the USB-6451

%% --- Staircase definition -------------------------------------------------
% Levels the fan is stepped through on the way UP. The down-leg walks the
% same levels in reverse, so hysteresis (does 4V on the way down give the
% same flow as 4V on the way up?) shows up as a split in the E-U curve.
MODE = 'stall';   % 'breakaway' | 'stall' | 'cal'

KICK_V = 2.3;     % V, start-up pulse (a little above the ~2 V threshold)
KICK_T = 2;       % s, long enough to see the fan turn; 0 = no kick

switch MODE
    case 'breakaway'
        KICK_T = 0;                              % must start from rest
        FAN_V_LEVELS = 1.5:0.05:2.2;
    case 'stall'
        FAN_V_LEVELS = 2.2:-0.05:1.2;            % walked in this (descending) order
    case 'cal'
        % Start at the lowest level that kept the fan spinning in 'stall'.
        FAN_V_LEVELS = [1.6 1.8 2 2.2 2.5 3 3.5 4 5 6 7 8];
    otherwise
        error('MODE must be ''breakaway'', ''stall'' or ''cal''.');
end
INCLUDE_DOWN   = false;                  % also step back down through the same levels
SKIP_REPEAT_TOP = true;                 % don't re-hold the top level at the start of the down-leg

STEP_SETTLE_T  = 20;   % s, discarded after each step change (flow settles in ~18 s)
STEP_DWELL_T   = 30;   % s, averaged into the calibration point
STEP_TRANSITION_T = 0; % s, ramp time BETWEEN levels (0 = instant step)
STEP_TRANSITION_DT = 0.25; % s, sub-step interval used during a transition ramp

PRE_RUN_ZERO_T  = 10;  % s, fan at 0 V before the staircase -- gives the zero-flow point
POST_RUN_ZERO_T = 5;   % s, fan back at 0 V at the end (recorded, not used in the fit)

HOTWIRE_FS = 4000;   % Hz

%% --- Channel selection ----------------------------------------------------
% Same group idea as run_experiment.m: Probe1-3 are the three sensors of one
% tri-axial probe and are always logged together.
allAiConfig = struct( ...
    'group', {'Hotwire',      'Hotwire',      'Hotwire',      'RefProbe'}, ...
    'name',  {'Probe1',       'Probe2',       'Probe3',       'RefProbe'}, ...
    'chan',  {"ai0",          "ai1",          "ai2",          "ai4"}, ...
    'term',  {"Differential", "Differential", "Differential", "Differential"}, ...
    'range', {[-10 10],       [-10 10],       [-10 10],       [-10 10]});

ENABLED_AI = ask_channels();   % interactive prompt (local function, bottom of file)

aiConfig   = allAiConfig(ismember({allAiConfig.group}, ENABLED_AI));
aiNames    = {aiConfig.name};
aiGroups   = {aiConfig.group};
hasHotwire = any(strcmp(aiGroups, 'Hotwire'));
hasRef     = any(strcmp(aiGroups, 'RefProbe'));

if isempty(aiConfig)
    error('No channels selected.');
end
if ~hasRef
    warning(['No reference probe selected -- there is no velocity reference, ' ...
             'so this run can only show hot-wire step response vs FAN VOLTAGE, ' ...
             'not a real E->U calibration.']);
end

%% --- Calibration constants ------------------------------------------------
CAL_FILE = 'D:\Chris\osbl-turbulent-mixing\experiment\data\260909\probe4.txt'; % this probe's cal/header
CTA_DIR  = 'D:\Chris\osbl-turbulent-mixing\experiment\processing_code\CTA';
addpath(CTA_DIR);

% The reference-probe calibration is the certificate table built into
% convert_Eref2Uref.m (54T29 s/n 0202, cert T29-202) -- nothing to set here.
% cal_ref is returned BY that call below and saved with the run, so every data
% file records which certificate produced its velocities.
cal_ref = [];

POLY_ORDER = 4;   % matches the 4th-order form used by convert_E2U_fn / convert_Eref2Uref

%% --- Build the step schedule ---------------------------------------------
% One row per hold. Built up front so the whole run is known (and printable)
% before any hardware is touched.
levelsUp = FAN_V_LEVELS(:)';
if INCLUDE_DOWN
    levelsDown = fliplr(levelsUp);
    if SKIP_REPEAT_TOP
        levelsDown(1) = [];
    end
else
    levelsDown = [];
end

stepV   = [0, levelsUp, levelsDown, 0];
stepDir = [{'zero'}, repmat({'up'}, 1, numel(levelsUp)), ...
           repmat({'down'}, 1, numel(levelsDown)), {'zero'}];
% The two zero holds get their own durations; every staircase level gets the
% standard settle+dwell.
stepHoldT = [PRE_RUN_ZERO_T, ...
             repmat(STEP_SETTLE_T + STEP_DWELL_T, 1, numel(levelsUp) + numel(levelsDown)), ...
             POST_RUN_ZERO_T];
% Averaging window inside each hold: skip the settle, take the rest. For the
% zero holds the settle is scaled down so a short hold still yields a point.
stepSettleT = min(STEP_SETTLE_T, 0.5 * stepHoldT);
stepSettleT(2:end-1) = STEP_SETTLE_T;

% Kick: inserted right after the leading zero hold. Settle = whole hold, so
% its averaging window is empty and it never becomes a calibration point.
if KICK_T > 0
    stepV       = [stepV(1),       KICK_V, stepV(2:end)];
    stepDir     = [stepDir(1),     {'kick'}, stepDir(2:end)];
    stepHoldT   = [stepHoldT(1),   KICK_T, stepHoldT(2:end)];
    stepSettleT = [stepSettleT(1), KICK_T, stepSettleT(2:end)];
end

nSteps = numel(stepV);
totalT = sum(stepHoldT) + max(nSteps-1, 0) * STEP_TRANSITION_T;
fprintf('Calibration staircase: %d holds, ~%.1f s total (%.1f min)\n', ...
    nSteps, totalT, totalT/60);
fprintf('  levels (V): %s\n', strjoin(compose('%.2f', stepV), ' '));

%% --- Set up DAQ -----------------------------------------------------------
fprintf('Configuring analog inputs (%s):\n', strjoin(aiNames, ', '));
hw = setup_hotwire_daq(DEV_ID, HOTWIRE_FS, aiConfig);

fanD = daq("ni");
addoutput(fanD, DEV_ID, "ao0", "Voltage");

%% --- Run the staircase ----------------------------------------------------
% Same timing approach as run_experiment.m: every boundary is scheduled
% against an ABSOLUTE elapsed time from the first fan write, so a slow
% iteration is absorbed by the next wait instead of accumulating drift.
t0 = tic;
sentT = [];   % every voltage actually written, on the t0 clock
sentV = [];

start(hw, "continuous");
hwStartElapsed = toc(t0);
disp('Hot-wire acquisition started.');

% stepStart/stepEnd are on the "wind clock" (0 = first fan write), filled in
% as the run progresses.
stepStart = nan(1, nSteps);
stepEnd   = nan(1, nSteps);

try
    % First write defines the wind clock, exactly as in run_experiment.m.
    write(fanD, stepV(1));
    fanStartElapsed = toc(t0);
    sentT(end+1) = fanStartElapsed; sentV(end+1) = stepV(1); %#ok<SAGROW>
    fprintf('[%2d/%2d] %-4s %.2f V, hold %.1f s\n', 1, nSteps, stepDir{1}, stepV(1), stepHoldT(1));

    cursor = 0;                      % elapsed on the wind clock
    stepStart(1) = cursor;
    cursor = cursor + stepHoldT(1);
    stepEnd(1) = cursor;
    wait_until(t0, fanStartElapsed + cursor);

    for k = 2:nSteps
        % Transition to the next level: instant by default, or a short ramp
        % if STEP_TRANSITION_T > 0 (gentler on the fan for big jumps).
        if STEP_TRANSITION_T > 0
            nSub = max(round(STEP_TRANSITION_T / STEP_TRANSITION_DT), 1);
            subV = linspace(stepV(k-1), stepV(k), nSub + 1);
            for s = 2:numel(subV)
                wait_until(t0, fanStartElapsed + cursor + (s-1)*STEP_TRANSITION_T/nSub);
                write(fanD, subV(s));
                sentT(end+1) = toc(t0); sentV(end+1) = subV(s); %#ok<SAGROW>
            end
            cursor = cursor + STEP_TRANSITION_T;
        else
            write(fanD, stepV(k));
            sentT(end+1) = toc(t0); sentV(end+1) = stepV(k); %#ok<SAGROW>
        end

        fprintf('[%2d/%2d] %-4s %.2f V, hold %.1f s (avg last %.1f s)\n', ...
            k, nSteps, stepDir{k}, stepV(k), stepHoldT(k), stepHoldT(k) - stepSettleT(k));

        stepStart(k) = cursor;
        cursor = cursor + stepHoldT(k);
        stepEnd(k) = cursor;
        wait_until(t0, fanStartElapsed + cursor);
    end

    write(fanD, 0);
    sentT(end+1) = toc(t0); sentV(end+1) = 0; %#ok<SAGROW>
    disp('Staircase done -- fan at 0 V.');
catch ME
    disp('Interrupted -- setting fan to 0V.');
    write(fanD, 0);
    stop(hw);
    disp('Hot-wire acquisition stopped (interrupted).');
    rethrow(ME);
end

stop(hw);
hwStopElapsed = toc(t0);
disp('Hot-wire acquisition stopped.');

%% --- Read back ------------------------------------------------------------
data     = read(hw, "all");
hwT      = seconds(data.Time);        % 0 = first acquired sample
hwT_t0   = hwT + hwStartElapsed;      % master t0 clock (matches sentT)
hwT_wind = hwT_t0 - fanStartElapsed;  % 0 = first fan write

E1 = []; E2 = []; E3 = []; U = []; V = []; W = [];
E_ref = []; U_ref = [];

if hasHotwire
    E1 = data{:, strcmp(aiNames, 'Probe1')};
    E2 = data{:, strcmp(aiNames, 'Probe2')};
    E3 = data{:, strcmp(aiNames, 'Probe3')};
    % Velocities from the EXISTING calibration file -- useful as a sanity
    % check (does the old cal still reproduce U_ref?), not as the new cal.
    [U, V, W] = convert_E2U_fn(E1, E2, E3, CAL_FILE);
end
if hasRef
    E_ref = data{:, strcmp(aiNames, 'RefProbe')};
    % Same source as run_experiment.m: CAL_FILE's T29 block (parse_probe_section),
    % NOT convert_Eref2Uref's built-in factory certificate -- see
    % run_hotwire_calibration.m for why the two must not be mixed.
    calRef = parse_probe_section(CAL_FILE, 'T29');
    U_ref  = polyval(fliplr(calRef.C(1,:)), E_ref);
    belowRef = E_ref < calRef.E_floor(1);
    U_ref(belowRef) = (calRef.U_floor(1) / calRef.E_floor(1)) * E_ref(belowRef);
    aboveRef = E_ref > calRef.E_ceil(1);
    U_ref(aboveRef) = calRef.U_ceil(1);
    cal_ref = calRef;
    fprintf('Reference probe (T29, %s): %.4f..%.4f V  ->  %.3f..%.3f m/s (%.1f%% below floor, %.1f%% above ceiling)\n', ...
        CAL_FILE, min(E_ref), max(E_ref), min(U_ref), max(U_ref), ...
        100*nnz(belowRef)/numel(E_ref), 100*nnz(aboveRef)/numel(E_ref));
end

%% --- Per-step averaging ---------------------------------------------------
% One calibration point per hold: mean over [stepStart+settle, stepEnd].
avgStart = stepStart + stepSettleT;
avgEnd   = stepEnd;

nCh = numel(aiNames);
stepMean = nan(nSteps, nCh);
stepStd  = nan(nSteps, nCh);
stepN    = zeros(nSteps, 1);
stepUrefMean = nan(nSteps, 1);
stepUrefStd  = nan(nSteps, 1);
stepUmean = nan(nSteps, 1);   % hot-wire U from the OLD cal file, for comparison

allData = data{:,:};
for k = 1:nSteps
    m = hwT_wind >= avgStart(k) & hwT_wind < avgEnd(k);
    stepN(k) = sum(m);
    if stepN(k) == 0
        continue;
    end
    stepMean(k,:) = mean(allData(m,:), 1);
    stepStd(k,:)  = std(allData(m,:), 0, 1);
    if hasRef
        stepUrefMean(k) = mean(U_ref(m));
        stepUrefStd(k)  = std(U_ref(m));
    end
    if hasHotwire
        stepUmean(k) = mean(U(m));
    end
end

stepTable = table(stepV(:), string(stepDir(:)), stepStart(:), stepEnd(:), ...
    avgStart(:), avgEnd(:), stepN, ...
    'VariableNames', {'FanV','Direction','tStart','tEnd','tAvgStart','tAvgEnd','nSamples'});
for i = 1:nCh
    stepTable.(['mean_' aiNames{i}]) = stepMean(:,i);
    stepTable.(['std_'  aiNames{i}]) = stepStd(:,i);
end
if hasRef
    stepTable.U_ref     = stepUrefMean;
    stepTable.U_ref_std = stepUrefStd;
end
if hasHotwire
    stepTable.U_oldcal = stepUmean;
end

disp(' ');
disp('Per-step averages:');
disp(stepTable);

%% --- Fit E -> U_ref per sensor -------------------------------------------
% The actual deliverable when both groups are logged: a 4th-order polynomial
% per hot-wire sensor, same functional form as the C0..C4 coefficients in the
% StreamWare calibration file that convert_E2U_fn.m reads. These are NOT
% written back into the cal file automatically -- copy them in yourself after
% checking the fit residuals below.
calFit = struct([]);
if hasHotwire && hasRef
    valid = stepN > 0 & ~isnan(stepUrefMean);
    Uc = stepUrefMean(valid);
    sensorNames = {'Probe1','Probe2','Probe3'};
    fprintf('\n%d-th order fits, U_ref = C0 + C1*E + ... (valid points: %d)\n', ...
        POLY_ORDER, sum(valid));
    if sum(valid) <= POLY_ORDER
        warning(['Only %d valid steps for an order-%d fit -- add more levels ' ...
                 'or lower POLY_ORDER.'], sum(valid), POLY_ORDER);
    else
        for i = 1:numel(sensorNames)
            col = strcmp(aiNames, sensorNames{i});
            Ec  = stepMean(valid, col);
            p   = polyfit(Ec, Uc, POLY_ORDER);   % highest power first
            res = Uc - polyval(p, Ec);
            c   = fliplr(p);                     % C0..C4, lowest power first
            calFit(i).sensor = sensorNames{i}; %#ok<SAGROW>
            calFit(i).E = Ec;
            calFit(i).U = Uc;
            calFit(i).p = p;
            calFit(i).C = c;
            calFit(i).rms_residual = sqrt(mean(res.^2));
            calFit(i).max_residual = max(abs(res));
            fprintf('  %s: %s   RMS %.4f m/s, max %.4f m/s\n', sensorNames{i}, ...
                strjoin(compose('C%d=%.6g', (0:POLY_ORDER)', c(:)), ' '), ...
                calFit(i).rms_residual, calFit(i).max_residual);
        end
    end
elseif hasHotwire
    disp('No reference probe -- skipping E->U fit (fan voltage is not a velocity).');
end

%% --- Plots ----------------------------------------------------------------
% 1) Full time series with each averaging window shaded, so it's obvious
%    whether STEP_SETTLE_T was long enough.
figure('Name','Calibration run -- time series');
nSub = double(hasHotwire) + double(hasRef) + 1;
iSub = 0;

iSub = iSub + 1;
subplot(nSub,1,iSub);
stairs(sentT - fanStartElapsed, sentV, 'LineWidth', 1.2);
shade_windows(avgStart, avgEnd);
xlabel('Time since first fan write (s)'); ylabel('Fan (V)');
title('Fan command staircase (shaded = averaged)'); grid on;

if hasHotwire
    iSub = iSub + 1;
    subplot(nSub,1,iSub);
    plot(hwT_wind, E1, hwT_wind, E2, hwT_wind, E3, 'LineWidth', 0.8);
    shade_windows(avgStart, avgEnd);
    xlabel('Time since first fan write (s)'); ylabel('Voltage (V)');
    legend('E1','E2','E3'); title('Hot-wire raw voltages'); grid on;
end

if hasRef
    iSub = iSub + 1;
    subplot(nSub,1,iSub);
    yyaxis left;  plot(hwT_wind, E_ref, 'LineWidth', 0.8); ylabel('E_{ref} (V)');
    yyaxis right; plot(hwT_wind, U_ref, 'LineWidth', 0.8); ylabel('U_{ref} (m/s)');
    shade_windows(avgStart, avgEnd);
    xlabel('Time since first fan write (s)');
    title('Reference probe'); grid on;
end

% Same time axis on every subplot: the fan stairs stop at the last write while
% the probe trace runs a little longer, so the auto limits would differ.
xl = [min([sentT(:) - fanStartElapsed; hwT_wind(:)]), ...
      max([sentT(:) - fanStartElapsed; hwT_wind(:)])];
set(findobj(gcf, 'Type', 'axes'), 'XLim', xl);

% 2) The calibration curve itself, up-leg and down-leg drawn separately so
%    hysteresis is visible rather than averaged away.
if hasRef
    figure('Name','Calibration curve');
    isUp   = strcmp(stepDir(:), 'up')   & stepN > 0;
    isDown = strcmp(stepDir(:), 'down') & stepN > 0;
    isZero = strcmp(stepDir(:), 'zero') & stepN > 0;

    if hasHotwire
        sensorNames = {'Probe1','Probe2','Probe3'};
        for i = 1:3
            subplot(2,2,i);
            col = strcmp(aiNames, sensorNames{i});
            hold on;
            plot(stepMean(isUp,col),   stepUrefMean(isUp),   'o-', 'DisplayName','up');
            plot(stepMean(isDown,col), stepUrefMean(isDown), 's--','DisplayName','down');
            plot(stepMean(isZero,col), stepUrefMean(isZero), 'kx', 'DisplayName','zero');
            if numel(calFit) >= i && ~isempty(calFit(i).p)
                Efit = linspace(min(calFit(i).E), max(calFit(i).E), 200);
                plot(Efit, polyval(calFit(i).p, Efit), 'k-', 'DisplayName','fit');
            end
            hold off;
            xlabel(sprintf('E_%d (V)', i)); ylabel('U_{ref} (m/s)');
            title(sprintf('Sensor %d', i)); legend('Location','best'); grid on;
        end
        subplot(2,2,4);
    end

    hold on;
    plot(stepV(isUp),   stepUrefMean(isUp),   'o-', 'DisplayName','up');
    plot(stepV(isDown), stepUrefMean(isDown), 's--','DisplayName','down');
    hold off;
    xlabel('Fan command (V)'); ylabel('U_{ref} (m/s)');
    title('Fan voltage -> reference velocity'); legend('Location','best'); grid on;
end

%% --- Save (prompted) ------------------------------------------------------
% Named hotwire_cal_* so it still matches plot_hotwire_data.m's hotwire_*
% file filter, and the variable names it expects are all present.
saveName = sprintf('hotwire_cal_kick_%s_%s.mat', MODE, datestr(now,'yyyymmdd_HHMMSS'));
saveVars = {'data','hwT','hwT_t0','hwT_wind','CAL_FILE','cal_ref', ...
            'E1','E2','E3','U','V','W','E_ref','U_ref', ...
            'aiConfig','aiNames','aiGroups','sentT','sentV', ...
            'stepV','stepDir','stepStart','stepEnd','avgStart','avgEnd', ...
            'stepMean','stepStd','stepN','stepTable','calFit', ...
            'FAN_V_LEVELS','MODE','KICK_V','KICK_T','STEP_SETTLE_T','STEP_DWELL_T','HOTWIRE_FS', ...
            'fanStartElapsed','hwStartElapsed','hwStopElapsed'};

drawnow;
resp = input(sprintf('\nSave this calibration run to %s? [Y/n]: ', saveName), 's');
if isempty(resp) || strncmpi(strtrim(resp), 'y', 1)
    save(saveName, saveVars{:});
    csvName = strrep(saveName, '.mat', '_steps.csv');
    writetable(stepTable, csvName);
    fprintf('Saved %s and %s\n', saveName, csvName);
else
    disp('Not saved. All variables are still in the workspace -- to save later, run:');
    fprintf('  save(saveName, saveVars{:})\n');
end


%% Local functions
function groups = ask_channels()
% Prompt for which probe group(s) to log. Loops until a valid answer is
% given -- a typo here would otherwise silently log the wrong channels for
% the whole run.
    while true
        fprintf('Which channels?\n');
        fprintf('  1) Hot-wire only     (ai0:ai2)\n');
        fprintf('  2) Reference probe only (ai4)\n');
        fprintf('  3) Both  [default]\n');
        r = strtrim(input('Choice [1/2/3]: ', 's'));
        switch lower(r)
            case {'1','h','hotwire','hw'}
                groups = {'Hotwire'};      return;
            case {'2','r','ref','refprobe'}
                groups = {'RefProbe'};     return;
            case {'3','b','both',''}
                groups = {'Hotwire','RefProbe'}; return;
            otherwise
                fprintf('Not understood: "%s"\n\n', r);
        end
    end
end

function wait_until(t0, targetElapsed)
% Block until toc(t0) reaches targetElapsed. Absolute targets (not relative
% pauses) keep the staircase on schedule across the whole run.
    CHECK_INTERVAL = 0.05;
    while toc(t0) < targetElapsed
        pause(max(0, min(CHECK_INTERVAL, targetElapsed - toc(t0))));
    end
end

function shade_windows(tStart, tEnd)
% Grey band behind each averaging window on the current axes.
%
% Uses XREGION (R2023b+), which draws behind the data automatically and
% rescales with the axes. The earlier patch + uistack version errored on the
% reference-probe subplot: that axes uses YYAXIS, and uistack reorders by
% assigning to the axes' Children, which a yyaxis axes rejects ("Children may
% only be set to a permutation of itself") because each side exposes only its
% own children. The patch fallback below therefore does NOT reorder -- it
% relies on transparency instead, so it is safe on any axes type.
    ok = ~isnan(tStart) & ~isnan(tEnd);
    if ~any(ok), return; end
    tStart = tStart(ok); tEnd = tEnd(ok);

    wasHeld = ishold; hold on;
    try
        h = xregion(tStart(:), tEnd(:), ...
            'FaceColor', [0.45 0.45 0.45], 'FaceAlpha', 0.13);
        set(h, 'HandleVisibility', 'off');   % keep the bands out of legends
    catch
        yl = ylim;
        for k = 1:numel(tStart)
            patch([tStart(k) tEnd(k) tEnd(k) tStart(k)], [yl(1) yl(1) yl(2) yl(2)], ...
                [0.45 0.45 0.45], 'EdgeColor','none', 'FaceAlpha',0.13, ...
                'HandleVisibility','off');
        end
        ylim(yl);
    end
    if ~wasHeld, hold off; end
end
