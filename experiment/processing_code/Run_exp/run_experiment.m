%% run_experiment.m
% Top-level script: starts fan ramp + hot-wire at t=0 on the USB-6451 (Dev4).
% Camera counters (PCI-6621, Dev1) start DELAY_BEFORE_TRIG seconds after the
% fan starts, and stop as soon as the fan ramp itself ends -- no external
% trigger wire, no LabVIEW hand-off. Both DAQ sessions are owned by MATLAB;
% timing between them is just software delay (pause), same mechanism already
% used for the fan ramp's internal timing.
%
% ASSUMPTIONS (check these before running):
%   - USB-6451 (fan/hot-wire) is "Dev4" in NI MAX -- CHANGE if different
%   - PCI-6621 (camera counters) is "Dev1" in NI MAX -- CHANGE if different
%   - Fan on Dev4/ao0, hot-wire on Dev4/ai0:ai2
%   - Fan ramp: 1.5->8V over 60s, hold 30s, ramp down 5s -- adjust as needed
%   - Hot-wire: 4 kHz sample rate
%   - Camera counters: IR/ctr2 (50 Hz, duty arbitrary) + Color/ctr1 (50 Hz,
%     0.035 duty = 700us pulse width, confirmed requirement) + Mono/ctr3
%     (50 Hz, 0.25 duty -- confirm this value, currently a placeholder)
%
% NOT verified -- confirm before relying on this:
%   - Whether start(camD, "continuous") or plain start(camD) is correct for
%     this DAQ Toolbox version (see test_camera_triggers.m notes)

clear; clc;

DEV_ID = "Dev4";   % <-- confirm this matches NI MAX for the USB-6451

FAN_V_START = 1.5;
FAN_V_END   = 8;
FAN_RAMPUP_T  = 60;   % seconds, ramp up duration
FAN_RAMPDOWN_T  = 5;   % seconds, ramp up duration
FAN_HOLD_T  = 10;    % seconds, hold at V_END
FAN_DT      = 0.5;  % seconds, step interval

HOTWIRE_FS       = 1000;   % Hz
DELAY_BEFORE_TRIG = 10;    % seconds after wind start before starting camera counters
PRE_WIND_BASELINE_T = 0;   % seconds of zero-flow hot-wire data before the fan starts (0 to skip)

% Analog input channels on the USB-6451. Hot-wire probes are differential
% (ai0/ai1/ai2, hardware-paired with ai8/ai9/ai10). RefProbe is the velocity
% reference transducer.
% CONFIRM: RefProbe is SingleEnded here because only one pin (ai4) was given.
% If it is actually wired differentially it uses ai4/ai12 and 'term' must be
% "Differential" instead.
% Channels are selected by GROUP, not individually: Probe1-3 are the three
% sensors of one tri-axial probe and are always logged together (the U,V,W
% decomposition needs all three), so they cannot be split apart here.
allAiConfig = struct( ...
    'group', {'Hotwire',      'Hotwire',      'Hotwire',      'RefProbe'}, ...
    'name',  {'Probe1',       'Probe2',       'Probe3',       'RefProbe'}, ...
    'chan',  {"ai0",          "ai1",          "ai2",          "ai4"}, ...
    'term',  {"Differential", "Differential", "Differential", "Differential"}, ...
    'range', {[-10 10],       [-10 10],       [-10 10],       [-10 10]});

% Valid groups: 'Hotwire', 'RefProbe'. Use {'RefProbe'} for reference only,
% {'Hotwire'} to skip the reference, or both.
ENABLED_AI = {'Hotwire', 'RefProbe'};%'Hotwire',
aiConfig   = allAiConfig(ismember({allAiConfig.group}, ENABLED_AI));
aiNames    = {aiConfig.name};
aiGroups   = {aiConfig.group};
hasHotwire = any(strcmp(aiGroups, 'Hotwire'));
hasRef     = any(strcmp(aiGroups, 'RefProbe'));

if isempty(aiConfig)
    error('ENABLED_AI selected no channels -- valid groups are ''Hotwire'' and ''RefProbe''.');
end

DEV_CAMERAS = "Dev3";   % <-- confirm this matches NI MAX for the PCI-6621
% Color needs a specific 700us trigger pulse width (confirmed) -> at 50 Hz
% (20000us period), dutyCycle = 700/20000 = 0.035. IR's duty cycle doesn't
% matter (pulse width not critical for that camera), so 0.5 is an arbitrary
% default.
allCamConfig = struct( ...
    'name',      {'IR',   'Color',  'Mono'}, ...
    'ctr',       {'ctr2', 'ctr1', 'ctr3'}, ...
    'freq',      {50,     50,   50}, ...
    'dutyCycle', {0.5,    0.035, 0.035}, ...
    'delay',     {0,      0,    0});

ENABLED_CAMERAS = {'IR','Color','Mono'};   % <-- edit to trigger only certain cameras, e.g. {'IR'} or {} for none
camConfig = allCamConfig(ismember({allCamConfig.name}, ENABLED_CAMERAS));

CAL_FILE = 'D:\Chris\osbl-turbulent-mixing\experiment\data\260721\probe3.txt'; % <-- this probe's cal/header
CTA_DIR  = 'D:\Chris\osbl-turbulent-mixing\experiment\processing_code\CTA';     % parse_calibration + convert_E2U_fn
addpath(CTA_DIR);

%% Set up hot-wire AI task (background acquisition, starts immediately)
% No total-time is fed in: totalTime was never actually wired into a scan
% count in setup_hotwire_daq.m, so a guessed duration wasn't really
% controlling anything. Instead, hw runs "continuous" and is stopped
% explicitly at the same real-world moment as the cameras -- when the fan
% ramp ends -- same mechanism, no guessed duration needed.
fprintf('Configuring analog inputs (%s):\n', strjoin(aiNames, ', '));
hw = setup_hotwire_daq(DEV_ID, HOTWIRE_FS, aiConfig);

%% Set up camera trigger counters (configured now, started later at the DELAY_BEFORE_TRIG mark)
if ~isempty(camConfig)
    camD = setup_camera_triggers(DEV_CAMERAS, camConfig);
    fprintf('Camera trigger counters configured (%s).\n', strjoin({camConfig.name}, ', '));
else
    camD = [];
    disp('No cameras enabled (ENABLED_CAMERAS empty) -- skipping camera trigger setup.');
end

%% Set up fan AO channel
fanD = daq("ni");
addoutput(fanD, DEV_ID, "ao0", "Voltage");

%% Run fan ramp (blocking, foreground) while hot-wire logs in background.
% Cameras must start DELAY_BEFORE_TRIG seconds after t=0 -- since that falls
% DURING the ramp-up (e.g. 10s into a 60s ramp-up), run_fan_ramp.m is NOT
% used here (it blocks until the whole ramp finishes). The ramp logic is
% inlined instead so elapsed time can be checked between steps.
disp('Starting fan ramp...');
t0 = tic;
camerasStarted = false;
camStartElapsed = NaN;
camStopElapsed = NaN;
sentT = [];
sentV = [];

%start hot-wire
start(hw, "continuous");   % <-- NOT VERIFIED call signature, see header
hwStartElapsed = toc(t0);
disp('Hot-wire acquisition started.');

% Pre-wind baseline: the hot-wire is now genuinely acquiring, so this pause
% records PRE_WIND_BASELINE_T seconds of zero-flow reference before the fan
% starts. These samples appear as negative hwT_wind values.
if PRE_WIND_BASELINE_T > 0
    fprintf('Collecting %.1f s of pre-wind baseline...\n', PRE_WIND_BASELINE_T);
    pause(PRE_WIND_BASELINE_T);
end

try
    nStepsUp = round(FAN_RAMPUP_T / FAN_DT) + 1;
    nStepsDown = round(FAN_RAMPDOWN_T / FAN_DT) + 1;
    nHoldSteps = max(round(FAN_HOLD_T / FAN_DT), 1);

    % Every step is scheduled against an ABSOLUTE time from t0 (rather than
    % pausing a relative duration each iteration), so a slow step is absorbed
    % by the next one instead of pushing the whole ramp later. Without this
    % the ~10ms per-step overshoot accumulated to ~2s over the full ramp.
    holdDt    = FAN_HOLD_T / nHoldSteps;
    holdEndT  = FAN_RAMPUP_T + FAN_HOLD_T;

    disp('Ramping up...');
    rampUpV = linspace(FAN_V_START, FAN_V_END, nStepsUp);

    % The FIRST fan write is the true "wind start" -- it lands ~0.19 s after
    % t0 because start(hw,...) and setup run first. Anchoring both the ramp
    % schedule and the camera delay here (instead of t0) makes
    % DELAY_BEFORE_TRIG a true delay from wind start, removing a ~143 ms
    % systematic offset.
    write(fanD, rampUpV(1));
    fanStartElapsed = toc(t0);
    sentT(end+1) = fanStartElapsed; sentV(end+1) = rampUpV(1); %#ok<SAGROW>
    camTargetElapsed = fanStartElapsed + DELAY_BEFORE_TRIG;
    fprintf('%.2f V\n', rampUpV(1));

    for k = 2:nStepsUp
        [camerasStarted, camStartElapsed] = waitUntilWithCameraCheck(fanStartElapsed + (k-1)*FAN_DT, camerasStarted, camStartElapsed, t0, camTargetElapsed, camD);
        v = rampUpV(k);
        write(fanD, v);
        sentT(end+1) = toc(t0); sentV(end+1) = v; %#ok<SAGROW> logged at the write
        if mod(k, 10) == 0 || k == nStepsUp
            fprintf('%.2f V\n', v);   % throttled: printing every step blocked the camera-start check
        end
    end

    fprintf('Holding at %.1f V...\n', FAN_V_END);
    for i = 1:nHoldSteps
        [camerasStarted, camStartElapsed] = waitUntilWithCameraCheck(fanStartElapsed + FAN_RAMPUP_T + i*holdDt, camerasStarted, camStartElapsed, t0, camTargetElapsed, camD);
        sentT(end+1) = toc(t0); sentV(end+1) = FAN_V_END; %#ok<SAGROW>
    end

    disp('Ramping down...');
    rampDownV = linspace(FAN_V_END, FAN_V_START, nStepsDown);
    for j = 1:nStepsDown
        [camerasStarted, camStartElapsed] = waitUntilWithCameraCheck(fanStartElapsed + holdEndT + (j-1)*FAN_DT, camerasStarted, camStartElapsed, t0, camTargetElapsed, camD);
        v = rampDownV(j);
        write(fanD, v);
        sentT(end+1) = toc(t0); sentV(end+1) = v; %#ok<SAGROW> logged at the write
        fprintf('%.2f V\n', v);
    end

    disp('Fan ramp done.');
catch ME
    disp('Interrupted -- setting fan to 0V.');
    write(fanD, 0);
    if ~isempty(camD) && camerasStarted
        stop(camD);
        disp('Camera trigger counters stopped (interrupted).');
    end
    stop(hw);
    disp('Hot-wire acquisition stopped (interrupted).');
    rethrow(ME);
end

% Safety net: if DELAY_BEFORE_TRIG is longer than the whole ramp, start now.
[camerasStarted, camStartElapsed] = maybeStartCameras(camerasStarted, camStartElapsed, t0, camTargetElapsed, camD);

% Report camera-start accuracy. Times are quoted as delay from the FIRST FAN
% WRITE (true wind start), which is what DELAY_BEFORE_TRIG is meant to mean.
if ~isnan(camStartElapsed)
    camDelayActual = camStartElapsed - fanStartElapsed;
    camStartJitter = camStartElapsed - camTargetElapsed;
    fprintf('Camera start: %.4f s after wind start (target %.3f s, jitter %+.4f s)\n', ...
        camDelayActual, DELAY_BEFORE_TRIG, camStartJitter);
else
    camDelayActual = NaN;
    camStartJitter = NaN;
end

% Note: residual camera-start jitter (~50-70 ms) traces to occasional
% write(fanD,...) calls blocking that long -- measured at mean 2.4 ms /
% max 73 ms. That is the floor for software-timed starts here.

% Wind (fan ramp) has ended -- stop the cameras AND the hot-wire now (both
% end at the same real-world moment: when the fan ramp finishes). No
% acquisition of meaningful data happens past this point, so precise
% trigger cutoff timing isn't critical here -- clicking Stop Record in the
% camera software and getting it one more trigger is handled manually,
% outside this script.
camStopElapsed = toc(t0);
if ~isempty(camD)
    stop(camD);
    disp('Camera trigger counters stopped (fan ramp ended).');
end

stop(hw);
hwStopElapsed = toc(t0);
disp('Hot-wire acquisition stopped (fan ramp ended).');

%% Read back the buffered hot-wire data
data = read(hw, "all");
% data.Time starts at 0 at the FIRST ACQUIRED SAMPLE, not at t0 and not at
% wind start -- so it needs shifting onto the shared timeline. Sample spacing
% itself is hardware-clocked (verified: 380000 samples / 4000 Hz = 95.0000 s).
hwT      = seconds(data.Time);        % raw: 0 = first sample
hwT_t0   = hwT + hwStartElapsed;      % on the master t0 clock (matches sentT)
hwT_wind = hwT_t0 - fanStartElapsed;  % 0 = wind start; pre-wind samples negative
disp('Hot-wire acquisition complete.');

%% Pull out the enabled channels by name (column order follows aiConfig)
E_ref = [];
if hasRef
    E_ref = data{:, strcmp(aiNames, 'RefProbe')};
    fprintf('Reference probe: mean %.4f V, range %.4f..%.4f V\n', ...
        mean(E_ref), min(E_ref), max(E_ref));
end

%% Convert hot-wire voltages to velocity (same pipeline as convertE2U.m)
% Left empty when the Hotwire group is not enabled, rather than saving filler
% that could later be mistaken for real data.
U = []; V = []; W = [];
E1 = []; E2 = []; E3 = [];
if hasHotwire
    E1 = data{:, strcmp(aiNames, 'Probe1')};
    E2 = data{:, strcmp(aiNames, 'Probe2')};
    E3 = data{:, strcmp(aiNames, 'Probe3')};
    [U, V, W] = convert_E2U_fn(E1, E2, E3, CAL_FILE);
    disp('Hot-wire voltages converted to velocity (U,V,W).');
else
    disp('Hotwire group not enabled -- skipping velocity conversion.');
end

% Save raw voltages, time bases, converted velocity, and event timings.
% hwT_wind is the one to use for analysis: 0 = wind start. aiConfig/aiNames
% record which channels were actually logged in this run.
save(sprintf('hotwire_%s.mat', datestr(now,'yyyymmdd_HHMMSS')), ...
    'data', 'hwT', 'hwT_t0', 'hwT_wind', 'CAL_FILE','E1', 'E2', 'E3','U', 'V', 'W', 'E_ref', ...
    'aiConfig', 'aiNames', 'aiGroups', 'sentT', 'sentV', ...
    'fanStartElapsed', 'hwStartElapsed', 'hwStopElapsed', ...
    'camStartElapsed', 'camStopElapsed', 'camDelayActual', 'camStartJitter');

disp('Experiment complete.');

%% Plot actual experiment signals
% Trigger window reconstructed from the real recorded camStartElapsed /
% camStopElapsed (not a guess) -- drawn as a representative 50 Hz square
% wave using the first enabled camera's frequency. Per plot_experiment_signals.m's
% own caveat, this does NOT capture per-camera pulse-width differences
% (e.g. Color's exposure vs IR) -- it's a visual "cameras were active here"
% window, not a literal per-channel waveform.
% All series below are built on the master t0 clock, then shifted by
% fanStartElapsed at plot time so t=0 means wind start for every subplot.
tEnd = max(sentT(end), hwT_t0(end));

if ~isempty(camConfig) && ~isnan(camStartElapsed)
    camFreq = camConfig(1).freq;
    dtTrig = 1 / (camFreq * 200);
    trigT = 0 : dtTrig : tEnd;
    inWindow = trigT >= camStartElapsed & trigT <= camStopElapsed;
    trigV = double(mod(trigT, 1/camFreq) < (1/camFreq/2)) .* inWindow;
else
    trigT = [0, tEnd];
    trigV = [0, 0];
end

% Hot-wire acquisition indicator: 0 = not taking data, 1 = acquiring.
% Built as an exact square pulse from the recorded start/stop times (corner
% points only, so it renders as a clean step) rather than the raw voltage.
hwLogicT = [0, hwStartElapsed, hwStartElapsed, hwStopElapsed, hwStopElapsed, tEnd];
hwLogicV = [0, 0,              1,              1,             0,             0];

plot_experiment_signals('actual', ...
    'FanT', sentT - fanStartElapsed, 'FanV', sentV, ...
    'HotwireT', hwLogicT - fanStartElapsed, 'HotwireV', hwLogicV, ...
    'TrigT', trigT - fanStartElapsed, 'TrigV', trigV);
%% Per-channel traces (only for the channels actually logged this run)
hwT_wind=hwT;
nSub = 2*double(hasHotwire) + double(hasRef);
if nSub > 0
    figure;
    iSub = 0;

    if hasHotwire
        iSub = iSub + 1;
        subplot(nSub,1,iSub)
        plot(hwT_wind, E1, hwT_wind, E2, hwT_wind, E3, 'LineWidth', 1);
        xlabel('Time since wind start (s)');
        ylabel('Voltage (V)');
        legend('E1', 'E2', 'E3');
        title('Hot-wire raw voltages');
        grid on;

        iSub = iSub + 1;
        subplot(nSub,1,iSub)
        plot(hwT_wind, U, hwT_wind, V, hwT_wind, W, 'LineWidth', 1);
        xlabel('Time since wind start (s)');
        ylabel('Velocity (m/s)');
        legend('U', 'V', 'W');
        title('Converted velocity (probe coordinates)');
        grid on;
    end

    if hasRef
        iSub = iSub + 1;
        subplot(nSub,1,iSub)
        plot(hwT_wind, E_ref, 'LineWidth', 1);
        xlabel('Time since wind start (s)');
        ylabel('Voltage (V)');
        legend('E_{ref}');
        title('Reference probe (ai4)');
        grid on;
    end
end


%% Local functions
function [started, startElapsed] = waitUntilWithCameraCheck(targetElapsed, started, startElapsed, t0, delay, camD)
% Wait until toc(t0) reaches targetElapsed, checking the camera-start
% condition every ~20ms. Waiting to an absolute target (rather than pausing
% a relative duration) keeps the ramp on schedule -- a step that runs long
% is absorbed by the next wait instead of accumulating drift.
    CHECK_INTERVAL = 0.02; % seconds
    while toc(t0) < targetElapsed
        pause(max(0, min(CHECK_INTERVAL, targetElapsed - toc(t0))));
        [started, startElapsed] = maybeStartCameras(started, startElapsed, t0, delay, camD);
    end
    % Also check when the target had already passed (no loop iterations ran).
    [started, startElapsed] = maybeStartCameras(started, startElapsed, t0, delay, camD);
end

function [started, startElapsed] = maybeStartCameras(started, startElapsed, t0, delay, camD)
    if ~started && toc(t0) >= delay
        if ~isempty(camD)
            start(camD, "continuous");   % <-- NOT VERIFIED call signature, see header
            disp('Camera trigger counters started.');
        end
        started = true;
        startElapsed = toc(t0);
    end
end
