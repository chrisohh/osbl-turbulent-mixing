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
%     0.035 duty = 700us pulse width, confirmed requirement)
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
FAN_HOLD_T  = 30;    % seconds, hold at V_END
FAN_DT      = 0.5;  % seconds, step interval

HOTWIRE_FS       = 4000;   % Hz
DELAY_BEFORE_TRIG = 10;    % seconds after t=0 before starting camera counters

DEV_CAMERAS = "Dev3";   % <-- confirm this matches NI MAX for the PCI-6621
% Color needs a specific 700us trigger pulse width (confirmed) -> at 50 Hz
% (20000us period), dutyCycle = 700/20000 = 0.035. IR's duty cycle doesn't
% matter (pulse width not critical for that camera), so 0.5 is an arbitrary
% default.
allCamConfig = struct( ...
    'name',      {'IR',   'Color'}, ...
    'ctr',       {'ctr2', 'ctr1'}, ...
    'freq',      {50,     50}, ...
    'dutyCycle', {0.5,    0.035}, ...
    'delay',     {0,      0});

ENABLED_CAMERAS = {'IR','Color'};   % <-- edit to trigger only certain cameras, e.g. {'IR'} or {} for none
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
hw = setup_hotwire_daq(DEV_ID, HOTWIRE_FS, []);

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
disp('Hot-wire acquisition started.');

try
    nStepsUp = round(FAN_RAMPUP_T / FAN_DT) + 1;
    nStepsDown = round(FAN_RAMPDOWN_T / FAN_DT) + 1;
    nHoldSteps = max(round(FAN_HOLD_T / FAN_DT), 1);

    disp('Ramping up...');
    for v = linspace(FAN_V_START, FAN_V_END, nStepsUp)
        write(fanD, v);
        fprintf('%.2f V\n', v);
        pause(FAN_DT);
        sentT(end+1) = toc(t0); sentV(end+1) = v; %#ok<SAGROW>
        [camerasStarted, camStartElapsed] = maybeStartCameras(camerasStarted, camStartElapsed, t0, DELAY_BEFORE_TRIG, camD);
    end

    fprintf('Holding at %.1f V...\n', FAN_V_END);
    for i = 1:nHoldSteps
        pause(FAN_HOLD_T / nHoldSteps);
        sentT(end+1) = toc(t0); sentV(end+1) = FAN_V_END; %#ok<SAGROW>
        [camerasStarted, camStartElapsed] = maybeStartCameras(camerasStarted, camStartElapsed, t0, DELAY_BEFORE_TRIG, camD);
    end

    disp('Ramping down...');
    for v = linspace(FAN_V_END, FAN_V_START, nStepsDown)
        write(fanD, v);
        fprintf('%.2f V\n', v);
        pause(FAN_DT);
        sentT(end+1) = toc(t0); sentV(end+1) = v; %#ok<SAGROW>
        % [camerasStarted, camStartElapsed] = maybeStartCameras(camerasStarted, camStartElapsed, t0, DELAY_BEFORE_TRIG, camD);
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
[camerasStarted, camStartElapsed] = maybeStartCameras(camerasStarted, camStartElapsed, t0, DELAY_BEFORE_TRIG, camD);

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
disp('Hot-wire acquisition stopped (fan ramp ended).');

%% Read back the buffered hot-wire data
data = read(hw, "all");
hwT = seconds(data.Time);   % elapsed time (s) for each hot-wire sample, from t0
disp('Hot-wire acquisition complete.');

%% Convert raw voltages to velocity (same pipeline as convertE2U.m)
E1 = data{:, 1};   % Probe 1 (ai0)
E2 = data{:, 2};   % Probe 2 (ai1)
E3 = data{:, 3};   % Probe 3 (ai2)

[U, V, W] = convert_E2U_fn(E1, E2, E3, CAL_FILE);
disp('Hot-wire voltages converted to velocity (U,V,W).');

% Save raw voltages, time, and converted velocity together
save(sprintf('hotwire_%s.mat', datestr(now,'yyyymmdd_HHMMSS')), 'data', 'hwT', 'U', 'V', 'W');

disp('Experiment complete.');

%% Plot actual experiment signals
% Trigger window reconstructed from the real recorded camStartElapsed /
% camStopElapsed (not a guess) -- drawn as a representative 50 Hz square
% wave using the first enabled camera's frequency. Per plot_experiment_signals.m's
% own caveat, this does NOT capture per-camera pulse-width differences
% (e.g. Color's exposure vs IR) -- it's a visual "cameras were active here"
% window, not a literal per-channel waveform.
if ~isempty(camConfig) && ~isnan(camStartElapsed)
    camFreq = camConfig(1).freq;
    dtTrig = 1 / (camFreq * 200);
    trigT = 0 : dtTrig : max(sentT(end), hwT(end));
    inWindow = trigT >= camStartElapsed & trigT <= camStopElapsed;
    trigV = double(mod(trigT, 1/camFreq) < (1/camFreq/2)) .* inWindow;
else
    trigT = [0, max(sentT(end), hwT(end))];
    trigV = [0, 0];
end

plot_experiment_signals('actual', ...
    'FanT', sentT, 'FanV', sentV, ...
    'HotwireT', hwT, 'HotwireV', E1, ...
    'TrigT', trigT, 'TrigV', trigV);

figure;
plot(hwT, U, hwT, V, hwT, W, 'k', 'LineWidth', 1);
xlabel('Time (s)');
ylabel('Velocity (m/s)');
legend('U', 'V', 'W');
title('Converted velocity (probe coordinates)');
grid on;


%% Local functions
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
