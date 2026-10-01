%% run_wavegauge_fan_calibration.m
% Wave gauge + nadir Color camera run for comparing gauge elevation with the
% wave slope retrieved from the Color images. The fan follows the SAME ramp as
% Run_exp/run_experiment.m (linear-in-U through fan_transfer.mat, ease-in,
% ramp down), the wave gauge on ai7 logs the whole time, and the Color camera
% triggers from DELAY_BEFORE_TRIG after wind start until the fan ramp ends --
% one continuous clip, as in run_experiment.m.
%
% No averaging: the output is the gauge time series eta(t) plus the time of
% every camera frame on the same clock, and eta at each frame.
%
% FRAME TIMING -- two options:
%   - CAM_TRIG_CHAN = "" (default): frame times are reconstructed from the
%     software-timed counter start. run_experiment.m measured ~50-70 ms of
%     jitter on that start -- 10-20% of a 2-4 Hz wave period, so a
%     frame-by-frame slope vs eta comparison will carry that phase error.
%   - CAM_TRIG_CHAN = "ai3" (or any spare Dev4 input): tee the Dev3/ctr1
%     trigger into that input. Each frame's rising edge is then found in the
%     same acquisition as the gauge, to within one sample. FS is raised to
%     10 kHz so the 700 us pulse spans ~7 samples.
%
% This is NOT the gauge's static height calibration: that needs still water
% and is done by run_wavegauge_calibration.m (fan off). This script LOADS that
% calibration (WG_CAL_FILE) to turn volts into cm; without one, eta is in volts.
%
% ASSUMPTIONS (check before running):
%   - USB-6451 (fan + gauge) is "Dev4" in NI MAX; fan on Dev4/ao0
%   - Wave gauge on Dev4/ai7, SingleEnded (as in WaveGauge_receiver.m)
%   - PCI-6621 (camera counters) is "Dev3"; Color camera on ctr1, 50 Hz,
%     700 us exposure -- copied from run_experiment.m
%   - Variables fanD / hw / camD use the same names as run_experiment.m, so
%     Run_exp/cleanup_daq.m zeroes the fan and stops the counters after a
%     Ctrl-C.
%
% NOT verified:
%   - that the USB-6451 sustains 10 kHz on 2 channels (only needed when
%     CAM_TRIG_CHAN is set). If it errors, lower FS_TRIG to 5000 (still ~3.5
%     samples per 700 us pulse).
%   - start(camD, "continuous") call signature -- same caveat as run_experiment.m.

clear; clc;

DEV_ID      = "Dev4";   % <-- confirm this matches NI MAX for the USB-6451
DEV_CAMERAS = "Dev3";   % <-- confirm this matches NI MAX for the PCI-6621
RUN_EXP_DIR = 'D:\Chris\osbl-turbulent-mixing\experiment\processing_code\Run_exp';
addpath(RUN_EXP_DIR);   % setup_camera_triggers
RUN_TAG = '';           % optional short label for the saved filename

%% --- Wave gauge ------------------------------------------------------------
WG_CHAN  = "ai7";
WG_TERM  = "SingleEnded";   % as in WaveGauge_receiver.m
WG_RANGE = [-10 10];
FS       = 1000;            % Hz, gauge only
FS_TRIG  = 10000;           % Hz, used instead when CAM_TRIG_CHAN is set

%% --- Fan ramp (same settings as Run_exp/run_experiment.m) -----------------
FAN_V_START = 2.1;    % fan cycles on/off at 2.0-2.01 V -- start clear of that
FAN_V_END   = 9.1;    % U(9.1V) ~= 10.0 m/s. KNOWN RISK: water came out of the
                      % fan at 9V (2026-09-23) -- watch for it near the top.
FAN_RAMP_LINEAR_U = true;
FAN_TF_FILE = fullfile(RUN_EXP_DIR, 'fan_transfer.mat');
FAN_V_START_HOLD_T = 2;    % s, hold at FAN_V_START before climbing
FAN_RAMP_EASE_T    = 5;    % s, ease-in from zero slope at the start of the climb
FAN_RAMPUP_T       = 70;   % s, TOTAL wind start -> FAN_V_END, hold included
FAN_RAMP_STOP_T    = Inf;  % s since wind start: cut the same ramp early (Inf = full)
FAN_RAMPDOWN_T     = 5;    % s
FAN_HOLD_T         = 0;    % s, hold at the top
FAN_DT             = 0.5;  % s, step interval

PRE_WIND_BASELINE_T = 10;  % s of still-water gauge data before the fan starts (0 to skip)
POST_WIND_LOG_T     = 0;   % s of gauge data after the fan reaches 0 V (waves decaying)

%% --- Camera (nadir Color) -------------------------------------------------
% Same values as run_experiment.m: pulse width IS the exposure (camera
% triggers on HIGH level), 0.035 at 50 Hz = 700 us.
allCamConfig = struct( ...
    'name',      {'IR',   'Color',  'Mono'}, ...
    'ctr',       {'ctr2', 'ctr1',   'ctr3'}, ...
    'freq',      {50,     50,       50}, ...
    'dutyCycle', {0.5,    0.035,    0.035}, ...
    'delay',     {0,      0,        0});
ENABLED_CAMERAS   = {'Color'};
DELAY_BEFORE_TRIG = 10;    % s after wind start before the camera starts (as run_experiment.m)

CAM_TRIG_CHAN   = "";      % e.g. "ai3" if the ctr1 trigger is teed into Dev4 -- see header
CAM_TRIG_TERM   = "SingleEnded";
CAM_TRIG_THRESH = 2.5;     % V, TTL rising-edge threshold

%% --- Wave gauge calibration (volts -> cm) ---------------------------------
% Empty = newest wave_gauge_calibration_*.mat found next to this script, in
% Run_exp, or in MATLAB's current folder (older runs saved to the current folder).
WG_CAL_FILE = '';
if isempty(WG_CAL_FILE)
    searchDirs = unique({fileparts(mfilename('fullpath')), RUN_EXP_DIR, pwd});
    f = [];
    for i = 1:numel(searchDirs)
        f = [f; dir(fullfile(searchDirs{i}, 'wave_gauge_calibration_*.mat'))]; %#ok<AGROW>
    end
    if ~isempty(f)
        [~, iNew] = max([f.datenum]);
        WG_CAL_FILE = fullfile(f(iNew).folder, f(iNew).name);
    end
end
wgCal = [];
if ~isempty(WG_CAL_FILE)
    S = load(WG_CAL_FILE, 'calibration');
    wgCal = S.calibration;
    fprintf('Wave gauge calibration: %s\n', WG_CAL_FILE);
else
    warning('No wave gauge calibration found -- eta will be in VOLTS.');
end

%% --- Camera on/off ----------------------------------------------------------
r = strtrim(input(sprintf('Trigger camera (%s)? [Y/n]: ', strjoin(ENABLED_CAMERAS, ', ')), 's'));
USE_CAMERA = isempty(r) || strncmpi(r, 'y', 1);
camConfig = allCamConfig(ismember({allCamConfig.name}, ENABLED_CAMERAS));
logTrig = USE_CAMERA && strlength(CAM_TRIG_CHAN) > 0;
if logTrig, FS = FS_TRIG; end

%% --- Build the fan schedule (wind clock, 0 = first fan write) -------------
% Same construction as run_experiment.m, built up front so the camera window
% is known before anything starts.
T_climb = FAN_RAMPUP_T - FAN_V_START_HOLD_T;
if T_climb <= FAN_RAMP_EASE_T
    error(['FAN_RAMPUP_T (%.1fs) - FAN_V_START_HOLD_T (%.1fs) = %.1fs of climb, ' ...
           'not longer than FAN_RAMP_EASE_T (%.1fs).'], ...
           FAN_RAMPUP_T, FAN_V_START_HOLD_T, T_climb, FAN_RAMP_EASE_T);
end
nStepsUp   = round(T_climb / FAN_DT) + 1;
nStepsDown = round(FAN_RAMPDOWN_T / FAN_DT) + 1;

fanTF = [];
if FAN_RAMP_LINEAR_U
    % Equal velocity increments, mapped back to voltage through the measured
    % V->U curve, with the same quadratic ease-in as run_experiment.m.
    S = load(FAN_TF_FILE, 'fanTF'); fanTF = S.fanTF;
    if FAN_V_START < fanTF.V(1) || FAN_V_END > fanTF.V(end)
        error('Ramp %.2f..%.2f V is outside the fan transfer range %.2f..%.2f V (%s).', ...
            FAN_V_START, FAN_V_END, fanTF.V(1), fanTF.V(end), FAN_TF_FILE);
    end
    uEnds = interp1(fanTF.V, fanTF.U, [FAN_V_START FAN_V_END], 'pchip');
    U0 = uEnds(1); U1 = uEnds(2);
    Te = FAN_RAMP_EASE_T;
    rU = (U1 - U0) / (T_climb - Te/2);   % m/s per s, plateau rate
    tSteps = linspace(0, T_climb, nStepsUp);
    rampUpU = nan(size(tSteps));
    easeMask = tSteps <= Te;
    rampUpU(easeMask)  = U0 + rU * tSteps(easeMask).^2 / (2*Te);
    rampUpU(~easeMask) = U0 + rU*Te/2 + rU*(tSteps(~easeMask) - Te);
    rampUpV = interp1(fanTF.U, fanTF.V, rampUpU, 'pchip');
    fprintf('Ramp linear in U: %.2f -> %.2f m/s over %.1f s climb (%.1f s ease-in), %.4f m/s/s\n', ...
        U0, U1, T_climb, Te, rU);
else
    rampUpV = linspace(FAN_V_START, FAN_V_END, nStepsUp);
end

tUp = FAN_V_START_HOLD_T + (0:nStepsUp-1) * FAN_DT;
tUp(1) = 0;                                   % first write IS wind start; the climb resumes after the hold
keepUp = [true, tUp(2:end) <= FAN_RAMP_STOP_T + 1e-9];
rampUpV = rampUpV(keepUp); tUp = tUp(keepUp);
rampEndT = FAN_V_START_HOLD_T + (numel(rampUpV)-1) * FAN_DT;
V_top    = rampUpV(end);
holdEndT = rampEndT + FAN_HOLD_T;

rampDownV = linspace(V_top, 0, nStepsDown);
tDown     = holdEndT + (0:nStepsDown-1) * FAN_DT;

schedT = [tUp, tDown];
schedV = [rampUpV, rampDownV];
fanEndT = schedT(end);

camOnSched  = DELAY_BEFORE_TRIG;
camOffSched = fanEndT;          % camera stops when the fan ramp ends, as in run_experiment.m
if USE_CAMERA
    fprintf('Camera (%s): %.1f -> %.1f s after wind start, %.1f s continuous, ~%d frames at %g Hz\n', ...
        strjoin({camConfig.name}, ', '), camOnSched, camOffSched, camOffSched - camOnSched, ...
        round((camOffSched - camOnSched) * camConfig(1).freq), camConfig(1).freq);
    if logTrig
        fprintf('Camera trigger logged on %s/%s -> frame times from the recorded pulses.\n', DEV_ID, CAM_TRIG_CHAN);
    else
        disp('Camera trigger NOT logged -- frame times from the software start (~50-70 ms uncertainty).');
    end
end
fprintf('Fan: %.2f -> %.2f V, ramp ends %.1f s, fan at 0 V at %.1f s after wind start.\n', ...
    FAN_V_START, V_top, rampEndT, fanEndT);

%% --- Set up DAQ -----------------------------------------------------------
hw = daq("ni");
ch = addinput(hw, DEV_ID, WG_CHAN, "Voltage");
ch.TerminalConfig = WG_TERM;
ch.Range = WG_RANGE;
if logTrig
    ch = addinput(hw, DEV_ID, CAM_TRIG_CHAN, "Voltage");
    ch.TerminalConfig = CAM_TRIG_TERM;
    ch.Range = [-10 10];
end
hw.Rate = FS;
fprintf('Wave gauge on %s/%s (%s) at %d Hz\n', DEV_ID, WG_CHAN, WG_TERM, FS);

fanD = daq("ni");
addoutput(fanD, DEV_ID, "ao0", "Voltage");
write(fanD, 0);   % known starting state (ao0 holds its last value between sessions)

camD = [];
if USE_CAMERA
    camD = setup_camera_triggers(DEV_CAMERAS, camConfig);
    input('Arm the camera software (waiting for trigger), then press Enter to start the run...', 's');
end

% Camera state carried through the waits; events fire against the wind clock.
cam = struct('d', {camD}, 'evT', {[camOnSched; camOffSched]}, 'evOn', {[true; false]}, ...
             'next', 1, 'origin', NaN, 'running', false, 'onT', [], 'offT', []);

%% --- Run --------------------------------------------------------------------
t0 = tic;
sentT = []; sentV = [];

start(hw, "continuous");
hwStartElapsed = toc(t0);
disp('Acquisition started.');

if PRE_WIND_BASELINE_T > 0
    fprintf('Collecting %.1f s of still-water baseline...\n', PRE_WIND_BASELINE_T);
    pause(PRE_WIND_BASELINE_T);
end

try
    write(fanD, schedV(1));
    fanStartElapsed = toc(t0);
    cam.origin = fanStartElapsed;
    sentT(end+1) = fanStartElapsed; sentV(end+1) = schedV(1); %#ok<SAGROW>
    fprintf('Wind start: %.2f V\n', schedV(1));

    for k = 2:numel(schedT)
        cam = wait_until(t0, fanStartElapsed + schedT(k), cam);
        write(fanD, schedV(k));
        sentT(end+1) = toc(t0); sentV(end+1) = schedV(k); %#ok<SAGROW>
        if mod(k, 20) == 0 || k == numel(schedT)
            fprintf('  t = %5.1f s   %.2f V\n', schedT(k), schedV(k));   % throttled, as in run_experiment.m
        end
    end
    disp('Fan ramp done -- fan at 0 V.');

    % Camera off event is at fanEndT, i.e. now.
    cam = service_cam(cam, t0);
    if cam.running   % e.g. DELAY_BEFORE_TRIG longer than the whole ramp -- close it
        stop(cam.d); cam.offT(end+1) = toc(t0) - cam.origin; cam.running = false;
    end

    if POST_WIND_LOG_T > 0
        fprintf('Logging %.1f s after the fan stops...\n', POST_WIND_LOG_T);
        wait_until(t0, fanStartElapsed + fanEndT + POST_WIND_LOG_T, cam);
    end
catch ME
    disp('Interrupted -- setting fan to 0V.');
    write(fanD, 0);
    if ~isempty(camD), stop(camD); disp('Camera trigger counters stopped (interrupted).'); end
    stop(hw);
    disp('Acquisition stopped (interrupted).');
    rethrow(ME);
end

stop(hw);
hwStopElapsed = toc(t0);
disp('Acquisition stopped.');

camOnT  = cam.onT;    % wind clock, actual (software) start/stop of the trigger
camOffT = cam.offT;
if ~isempty(camOnT)
    fprintf('Camera start: %.4f s after wind start (target %.3f s, jitter %+.4f s)\n', ...
        camOnT(1), camOnSched, camOnT(1) - camOnSched);
end

%% --- Read back ------------------------------------------------------------
data     = read(hw, "all");
hwT      = seconds(data.Time);
hwT_t0   = hwT + hwStartElapsed;
hwT_wind = hwT_t0 - fanStartElapsed;   % 0 = wind start; baseline samples negative

E_wg = data{:, 1};
if ~isempty(wgCal)
    eta = polyval(wgCal.polynomial_coefficients, E_wg);   % cm
    etaUnit = 'cm';
    vc = wgCal.voltage_calib;
    if isfield(wgCal, 'clipped'), vc = vc(~wgCal.clipped); end   % older cal files lack it
    outOfCal = E_wg < min(vc) | E_wg > max(vc);
    fprintf('Wave gauge: %.1f%% of samples outside the calibrated voltage range (extrapolated).\n', ...
        100*nnz(outOfCal)/numel(E_wg));
else
    eta = E_wg;
    etaUnit = 'V';
end
clipped = E_wg <= WG_RANGE(1) + 0.05 | E_wg >= WG_RANGE(2) - 0.05;
if any(clipped)
    warning('%.2f%% of wave gauge samples are at the +-10 V limit.', 100*nnz(clipped)/numel(E_wg));
end

%% --- Frame times --------------------------------------------------------------
% frameT: wind-clock time of each camera trigger (frame n <-> image n in the
% camera's saved sequence, assuming the camera saved every trigger).
frameT = []; frameSource = 'none'; E_trig = [];
if USE_CAMERA && ~isempty(camOnT)
    camFreq = camConfig(1).freq;
    nExpected = floor((camOffT(1) - camOnT(1)) * camFreq);
    if logTrig
        E_trig = data{:, 2};
        hi = E_trig > CAM_TRIG_THRESH;
        iRise = find(diff(hi) == 1) + 1;
        frameT = hwT_wind(iRise);
        frameSource = 'recorded trigger';
        dT = diff(frameT);
        fprintf('Frames: %d rising edges recorded (expected ~%d); interval %.4f +- %.4f s (max %.4f)\n', ...
            numel(frameT), nExpected, mean(dT), std(dT), max(dT));
        if any(abs(dT - 1/camFreq) > 0.5/camFreq)
            warning('Irregular trigger intervals -- check the tee / threshold (CAM_TRIG_THRESH).');
        end
    else
        frameT = camOnT(1) + (0:nExpected-1)' / camFreq;
        frameSource = 'software start (+-~50-70 ms)';
        fprintf('Frames: ~%d, times reconstructed from the software start.\n', nExpected);
    end
    etaAtFrame = interp1(hwT_wind, eta, frameT);
else
    etaAtFrame = [];
end
frameTable = table((1:numel(frameT))', frameT(:), etaAtFrame(:), ...
    'VariableNames', {'Frame', 't_wind_s', ['eta_' etaUnit]});

%% --- Plots ----------------------------------------------------------------
figure('Name','Wave ramp run');
nSub = 3 + double(logTrig);
ax = gobjects(nSub, 1);

ax(1) = subplot(nSub,1,1);
% plot (not stairs), as plot_experiment_signals.m draws run_experiment.m's
% fan -- the AO actually holds each write for FAN_DT, but the fan's inertia
% smooths the 0.9 V ramp-down steps anyway.
plot(sentT - fanStartElapsed, sentV, 'LineWidth', 1.2);
shade_camera(camOnT, camOffT);
ylabel('Fan (V)'); grid on;
title('Fan command (magenta = camera triggering)');

ax(2) = subplot(nSub,1,2);
plot(hwT_wind, eta, 'LineWidth', 0.5);
shade_camera(camOnT, camOffT);
ylabel(sprintf('\\eta (%s)', etaUnit)); title('Wave gauge'); grid on;

ax(3) = subplot(nSub,1,3);
% Frequency content through the ramp -- shows the peak moving as the waves grow.
winS = 4;   % s
% 4th output is the one-sided PSD (etaUnit^2/Hz); the first output |s|^2
% would be raw FFT power with no physical units.
[~, fS, tS, pS] = spectrogram(detrend(eta), hann(round(winS*FS)), round(winS*FS/2), [], FS);
imagesc(tS + hwT_wind(1), fS, 10*log10(pS + eps)); axis xy;
cb = colorbar;
cb.Label.String = sprintf('10 log_{10} S_\\eta(f)  (dB re 1 %s^2/Hz)', etaUnit);
ylim([0 10]); ylabel('f (Hz)'); title(sprintf('\\eta spectrogram (%g s windows)', winS));
shade_camera(camOnT, camOffT);

if logTrig
    ax(4) = subplot(nSub,1,4);
    plot(hwT_wind, E_trig); hold on;
    plot(frameT, CAM_TRIG_THRESH * ones(size(frameT)), 'r.', 'MarkerSize', 4); hold off;
    ylabel('Trigger (V)'); title('Recorded camera trigger (red = detected frames)'); grid on;
end
xlabel('Time since wind start (s)');
linkaxes(ax, 'x');
xlim(ax(1), [hwT_wind(1), hwT_wind(end)]);

%% --- Save (prompted) ------------------------------------------------------
if isempty(RUN_TAG)
    saveName = sprintf('wave_ramp_%s.mat', datestr(now,'yyyymmdd_HHMMSS'));
else
    saveName = sprintf('wave_ramp_%s_%s.mat', RUN_TAG, datestr(now,'yyyymmdd_HHMMSS'));
end
runConfig = struct( ...
    'RUN_TAG', RUN_TAG, 'DEV_ID', DEV_ID, 'WG_CHAN', WG_CHAN, 'WG_TERM', WG_TERM, 'FS', FS, ...
    'FAN_V_START', FAN_V_START, 'FAN_V_END', FAN_V_END, ...
    'FAN_RAMP_LINEAR_U', FAN_RAMP_LINEAR_U, 'FAN_TF_FILE', FAN_TF_FILE, ...
    'FAN_V_START_HOLD_T', FAN_V_START_HOLD_T, 'FAN_RAMP_EASE_T', FAN_RAMP_EASE_T, ...
    'FAN_RAMPUP_T', FAN_RAMPUP_T, 'FAN_RAMP_STOP_T', FAN_RAMP_STOP_T, ...
    'FAN_RAMPDOWN_T', FAN_RAMPDOWN_T, 'FAN_HOLD_T', FAN_HOLD_T, 'FAN_DT', FAN_DT, ...
    'rampEndT', rampEndT, 'V_top', V_top, 'fanEndT', fanEndT, ...
    'PRE_WIND_BASELINE_T', PRE_WIND_BASELINE_T, 'POST_WIND_LOG_T', POST_WIND_LOG_T, ...
    'USE_CAMERA', USE_CAMERA, 'ENABLED_CAMERAS', {ENABLED_CAMERAS}, ...
    'DELAY_BEFORE_TRIG', DELAY_BEFORE_TRIG, 'CAM_TRIG_CHAN', CAM_TRIG_CHAN, ...
    'CAM_TRIG_THRESH', CAM_TRIG_THRESH, 'fanTF', fanTF);
saveVars = {'data','hwT','hwT_t0','hwT_wind','E_wg','eta','etaUnit','WG_CAL_FILE','wgCal', ...
            'E_trig','frameT','frameSource','etaAtFrame','frameTable', ...
            'sentT','sentV','schedT','schedV','camConfig','camOnT','camOffT', ...
            'fanStartElapsed','hwStartElapsed','hwStopElapsed','runConfig'};

drawnow;
resp = input(sprintf('\nSave this run to %s? [Y/n]: ', saveName), 's');
if isempty(resp) || strncmpi(strtrim(resp), 'y', 1)
    save(saveName, saveVars{:});
    if ~isempty(frameT)
        csvName = strrep(saveName, '.mat', '_frames.csv');
        writetable(frameTable, csvName);
        fprintf('Saved %s and %s\n', saveName, csvName);
    else
        fprintf('Saved %s\n', saveName);
    end
else
    disp('Not saved. To save later, run:');
    fprintf('  save(saveName, saveVars{:})\n');
end


%% Local functions
function cam = wait_until(t0, targetElapsed, cam)
% Block until toc(t0) reaches targetElapsed, firing any camera start/stop
% that falls due in the meantime (same idea as run_experiment.m's
% waitUntilWithCameraCheck). Absolute targets keep the ramp on schedule.
    CHECK_INTERVAL = 0.02;
    while toc(t0) < targetElapsed
        pause(max(0, min(CHECK_INTERVAL, targetElapsed - toc(t0))));
        cam = service_cam(cam, t0);
    end
    cam = service_cam(cam, t0);
end

function cam = service_cam(cam, t0)
% Fire every camera event whose wind-clock time has passed.
    if isempty(cam.d) || isnan(cam.origin), return; end
    while cam.next <= numel(cam.evT) && toc(t0) - cam.origin >= cam.evT(cam.next)
        if cam.evOn(cam.next) && ~cam.running
            start(cam.d, "continuous");   % <-- NOT VERIFIED call signature (as in run_experiment.m)
            cam.onT(end+1) = toc(t0) - cam.origin;
            cam.running = true;
            fprintf('  camera ON  at %.2f s\n', cam.onT(end));
        elseif ~cam.evOn(cam.next) && cam.running
            stop(cam.d);
            cam.offT(end+1) = toc(t0) - cam.origin;
            cam.running = false;
            fprintf('  camera OFF at %.2f s\n', cam.offT(end));
        end
        cam.next = cam.next + 1;
    end
end

function shade_camera(tOn, tOff)
% Magenta band behind each camera window on the current axes.
    if isempty(tOn), return; end
    wasHeld = ishold; hold on;
    try
        h = xregion(tOn(:), tOff(:), 'FaceColor', [0.85 0.2 0.85], 'FaceAlpha', 0.12);
        set(h, 'HandleVisibility', 'off');
    catch
        yl = ylim;
        for k = 1:numel(tOn)
            patch([tOn(k) tOff(k) tOff(k) tOn(k)], [yl(1) yl(1) yl(2) yl(2)], ...
                [0.85 0.2 0.85], 'EdgeColor','none', 'FaceAlpha',0.12, 'HandleVisibility','off');
        end
        ylim(yl);
    end
    if ~wasHeld, hold off; end
end
