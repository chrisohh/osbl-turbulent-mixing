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
%   - PCI-6621 (camera counters) is "Dev1" in NI MAX -- CHANGE if differentnn
%   - Fan on Dev4/ao0, hot-wire on Dev4/ai0:ai2
%   - Fan ramp: 1.5->8V over 60s, hold 30s, ramp down 5s -- adjust as needed
%   - Hot-wire: 4 kHz sample rate
%   - Camera counters, all 50 Hz: IR/ctr2 (duty arbitrary -- IR sets its own
%     exposure), Color/ctr1 and Mono/ctr3 (0.035 duty = 700us exposure, since
%     these trigger on HIGH level so pulse width = exposure time)
%   - ctr3 feeds BOTH the Mono side camera and the Mono under-tank camera
%
% NOT verified -- confirm before relying on this:
%   - Whether start(camD, "continuous") or plain start(camD) is correct for
%     this DAQ Toolbox version (see test_camera_triggers.m notes)

clear; clc;

DEV_ID = "Dev4";   % <-- confirm this matches NI MAX for the USB-6451

% V->U: a straight line fits to <1% RMS of range over 2.1-9.5V, but the
% residuals are NOT noise -- a real, if mild, S-curve (see build_fan_transfer_
% combined.m's merged data). So PCHIP (via fan_transfer.mat) is used, not a
% line -- it tracks that shape instead of averaging over it.
FAN_V_START = 2.1;   % fan cycles on/off at 2.0-2.01 V -- start clear of that.
                      % Flow jumps 0 -> U(2.1)~1.85 m/s at wind start; the ramp
                      % below only covers that -> FAN_V_END's U, not literally 0.
FAN_V_END   = 9.1;    % U(9.1V) ~= 10.0 m/s (interpolated between 9V->9.91
                      % and 9.5V->10.35, hotwire_cal_20260923_034304).
                      % KNOWN RISK: water came out of the fan at 9V (user,
                      % 2026-09-23) -- this holds AT/ABOVE that voltage for
                      % part of FAN_HOLD_T (10s at V_END). User chose to
                      % accept this risk to hit 10 m/s (2026-09-23) rather
                      % than cap at 8.5V (~9.4 m/s, no water seen). Watch
                      % for water during FAN_HOLD_T; if it recurs, drop back
                      % to 8.5V and treat 10 m/s as not achievable with this
                      % hardware until the water source is fixed.
FAN_RAMP_LINEAR_U = true;
FAN_TF_FILE = 'D:\Chris\osbl-turbulent-mixing\experiment\processing_code\Run_exp\fan_transfer.mat';
% Held at FAN_V_START before the ramp begins climbing, so the 0->U(V_START)
% jump has time to settle (flow settles in ~18s after a step, same measured
% value used for STEP_SETTLE_T throughout calibration) before the ramp piles
% more voltage on top of it. Without this the first ~18-20s of "linear in U"
% would really be jump-settling dynamics plus ramp command, not a clean climb.
% DELAY_BEFORE_TRIG (10s) is measured from wind start as before (user,
% 2026-09-23) -- with the hold below, the cameras now start DURING it,
% before the ramp itself begins climbing. That's the chosen behavior, not a Ybug.
% NOTE: 10s is BELOW the ~18-20s settle time this rig actually measured
% (hotwire_cal_20260923_013607: plateau not reached until ~18-20s). Left as
% typed rather than silently reverted to 20 -- but the climb below will start
% from flow that likely hasn't fully settled yet at 10s.
FAN_V_START_HOLD_T = 2;   % s
% FAN_RAMPUP_T is the TOTAL time from wind start to reaching FAN_V_END,
% hold included (user, 2026-09-23) -- so the actual climb only gets
% FAN_RAMPUP_T - FAN_V_START_HOLD_T seconds. The first FAN_RAMP_EASE_T of
% that climb eases in from zero slope (matching the flat hold) up to the
% constant rate used for the rest, instead of jumping straight to full rate
% -- avoids a kink in acceleration at the hold->ramp transition.
FAN_RAMP_EASE_T = 5;   % s, see block above the ramp-up loop for the closed form
FAN_RAMPUP_T  = 70;   % seconds, TOTAL wind-start -> FAN_V_END, hold included
FAN_RAMPDOWN_T  = 5;   % seconds, ramp up duration
FAN_HOLD_T  = 0;    % seconds, hold at V_END
FAN_DT      = 0.5;  % seconds, step interval

HOTWIRE_FS       = 1000;   % Hz
DELAY_BEFORE_TRIG = 10;    % seconds after wind start before starting camera counters
PRE_WIND_BASELINE_T = 0;   % seconds of zero-flow hot-wire data before the fan startsn (0 to skip)

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
ENABLED_AI = {'Hotwire'};%'Hotwire','Hotwire', 'RefProbe
aiConfig   = allAiConfig(ismember({allAiConfig.group}, ENABLED_AI));
aiNames    = {aiConfig.name};
aiGroups   = {aiConfig.group};
hasHotwire = any(strcmp(aiGroups, 'Hotwire'));
hasRef     = any(strcmp(aiGroups, 'RefProbe'));

if isempty(aiConfig)
    error('ENABLED_AI selected no channels -- valid groups are ''Hotwire'' and ''RefProbe''.');
end

DEV_CAMERAS = "Dev3";   % <-- confirm this matches NI MAX for the PCI-6621
% Pulse width IS the exposure: these cameras are set to trigger on HIGH level,
% so they expose for as long as the line is high. At 50 Hz (20000us period),
% dutyCycle = exposure_us / 20000 -- so 0.035 = 700us.
% IR ignores pulse width (its exposure is set in its own software), so its
% 0.5 is arbitrary.
%
% ctr3 DRIVES TWO CAMERAS: the Mono side camera and the Mono under-tank
% cross-view camera are wired to the same counter output, so they necessarily
% share this signal -- same 50 Hz, same 700us exposure. Changing 'Mono' below
% changes both. To give the under-tank camera its own exposure it would need
% a separate counter (ctr0 is free) and its own entry here.
% ctr1 (Color) and ctr2 (IR) each drive a single camera.
allCamConfig = struct( ...
    'name',      {'IR',   'Color',  'Mono'}, ...
    'ctr',       {'ctr2', 'ctr1', 'ctr3'}, ...
    'freq',      {50,     50,   50}, ...
    'dutyCycle', {0.5,    0.035, 0.035}, ...
    'delay',     {0,      0,    0});

% 'Mono' here means the ctr3 line, i.e. Mono side + Mono under-tank.
ENABLED_CAMERAS = {'IR','Color','Mono'};   % <-- edit to trigger only certain cameras, e.g. {'IR'} or {} for none
camConfig = allCamConfig(ismember({allCamConfig.name}, ENABLED_CAMERAS));

% probe4.txt holds BOTH probes used below: the tri-axial hot-wire (55P95,
% read by parse_calibration.m/convert_E2U_fn.m, which take the file's FIRST
% probe block) and the reference probe (T29, read by parse_probe_section.m
% by name further down).
CAL_FILE = 'D:\Chris\osbl-turbulent-mixing\experiment\data\260909\probe4.txt'; % <-- this probe's cal/header
CTA_DIR  = 'D:\Chris\osbl-turbulent-mixing\experiment\processing_code\CTA';     % parse_calibration + convert_E2U_fn + parse_probe_section
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
% Force a known starting state -- ao0 holds its LAST commanded voltage
% between MATLAB sessions, so without this an interrupted or otherwise
% incomplete previous run could leave the fan spinning at whatever voltage
% it last had, and this run would silently jump from there instead of 0.
write(fanD, 0);

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
    % FAN_RAMPUP_T is TOTAL wind-start -> FAN_V_END (hold included), so the
    % climb itself only gets what's left over.
    T_climb = FAN_RAMPUP_T - FAN_V_START_HOLD_T;
    if T_climb <= FAN_RAMP_EASE_T
        error(['FAN_RAMPUP_T (%.1fs) - FAN_V_START_HOLD_T (%.1fs) = %.1fs of climb, ' ...
               'which is not longer than FAN_RAMP_EASE_T (%.1fs). Raise FAN_RAMPUP_T ' ...
               'or shorten the hold/ease.'], FAN_RAMPUP_T, FAN_V_START_HOLD_T, T_climb, FAN_RAMP_EASE_T);
    end
    nStepsUp = round(T_climb / FAN_DT) + 1;
    nStepsDown = round(FAN_RAMPDOWN_T / FAN_DT) + 1;
    nHoldSteps = max(round(FAN_HOLD_T / FAN_DT), 1);

    % Every step is scheduled against an ABSOLUTE time from t0 (rather than
    % pausing a relative duration each iteration), so a slow step is absorbed
    % by the next one instead of pushing the whole ramp later. Without this
    % the ~10ms per-step overshoot accumulated to ~2s over the full ramp.
    holdDt    = FAN_HOLD_T / nHoldSteps;
    holdEndT  = FAN_RAMPUP_T + FAN_HOLD_T;   % FAN_RAMPUP_T already includes the hold

    disp('Ramping up...');
    if FAN_RAMP_LINEAR_U
        % Equal velocity increments per step, mapped back to voltage through
        % the measured steady-state V->U curve. No extrapolation: both ends
        % must sit inside the calibrated voltage range.
        tf = load(FAN_TF_FILE, 'fanTF'); tf = tf.fanTF;
        if FAN_V_START < tf.V(1) || FAN_V_END > tf.V(end)
            error('Ramp %.2f..%.2f V is outside the fan transfer range %.2f..%.2f V (%s).', ...
                FAN_V_START, FAN_V_END, tf.V(1), tf.V(end), FAN_TF_FILE);
        end
        uEnds = interp1(tf.V, tf.U, [FAN_V_START FAN_V_END], 'pchip');
        U0 = uEnds(1); U1 = uEnds(2);
        % Ease-in (0 slope at t=0, matching the flat hold) blended into a
        % constant rate for the rest of the climb, reaching U1 exactly at
        % T_climb. Quadratic on [0,Tease]: U=U0 + r*t^2/(2*Tease), value+slope
        % match the linear piece on [Tease,T_climb]: U=U0+r*Tease/2+r*(t-Tease).
        % Solving U(T_climb)=U1 for r gives the denominator below.
        Te = FAN_RAMP_EASE_T;
        r = (U1 - U0) / (T_climb - Te/2);   % m/s per s, the plateau rate
        tSteps = linspace(0, T_climb, nStepsUp);
        rampUpU = nan(size(tSteps));
        easeMask = tSteps <= Te;
        rampUpU(easeMask)  = U0 + r * tSteps(easeMask).^2 / (2*Te);
        rampUpU(~easeMask) = U0 + r*Te/2 + r*(tSteps(~easeMask) - Te);
        rampUpV = interp1(tf.U, tf.V, rampUpU, 'pchip');
        fprintf('Ramp linear in U: %.2f -> %.2f m/s over %.1fs climb (%.1fs ease-in), rate %.4f m/s/s (%s)\n', ...
            U0, U1, T_climb, Te, r, tf.created);
    else
        rampUpV = linspace(FAN_V_START, FAN_V_END, nStepsUp);
    end

    % The FIRST fan write is the true "wind start" -- it lands ~0.19 s after
    % t0 because start(hw,...) and setup run first. Anchoring both the ramp
    % schedule and the camera delay here (instead of t0) makes
    % DELAY_BEFORE_TRIG a true delay from wind start, removing a ~143 ms
    % systematic offset.
    write(fanD, rampUpV(1));
    fanStartElapsed = toc(t0);
    sentT(end+1) = fanStartElapsed; sentV(end+1) = rampUpV(1); %#ok<SAGROW>
    camTargetElapsed = fanStartElapsed + DELAY_BEFORE_TRIG;   % from wind start, unaffected by the hold below
    fprintf('%.2f V\n', rampUpV(1));

    if FAN_V_START_HOLD_T > 0
        fprintf('Holding at %.2f V for %.1f s (letting the jump settle before ramping)...\n', ...
            FAN_V_START, FAN_V_START_HOLD_T);
        [camerasStarted, camStartElapsed] = waitUntilWithCameraCheck(fanStartElapsed + FAN_V_START_HOLD_T, camerasStarted, camStartElapsed, t0, camTargetElapsed, camD);
    end

    for k = 2:nStepsUp
        [camerasStarted, camStartElapsed] = waitUntilWithCameraCheck(fanStartElapsed + FAN_V_START_HOLD_T + (k-1)*FAN_DT, camerasStarted, camStartElapsed, t0, camTargetElapsed, camD);
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
    % Down to 0V, not just back to FAN_V_START -- the fan should start AND
    % end the run at 0, not be left spinning at FAN_V_START afterward. Below
    % ~2.0V nothing sustains flow anyway (dead band), so this last stretch of
    % the ramp is just as physically discontinuous on the way down as the
    % 0->FAN_V_START jump was on the way up.
    rampDownV = linspace(FAN_V_END, 0, nStepsDown);
    for j = 1:nStepsDown
        [camerasStarted, camStartElapsed] = waitUntilWithCameraCheck(fanStartElapsed + holdEndT + (j-1)*FAN_DT, camerasStarted, camStartElapsed, t0, camTargetElapsed, camD);
        v = rampDownV(j);
        write(fanD, v);
        sentT(end+1) = toc(t0); sentV(end+1) = v; %#ok<SAGROW> logged at the write
        fprintf('%.2f V\n', v);
    end

    disp('Fan ramp done -- fan at 0V.');
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
U_ref = [];
if hasRef
    E_ref = data{:, strcmp(aiNames, 'RefProbe')};

    % Reference-probe velocity from CAL_FILE's T29 block (parse_probe_section,
    % in CTA_DIR) instead of the certificate table hardcoded in
    % convert_Eref2Uref.m -- same cal file as the hot-wire's 55P95 block
    % above, different probe section. Below the calibrated floor the
    % polynomial diverges, so those samples fall back to the origin-to-floor
    % line (extrapolation, not calibration); above the ceiling U_ref is
    % clamped. Same handling as convert_E2U_fn.m uses for the hot-wire.
    calRef = parse_probe_section(CAL_FILE, 'T29');
    U_ref  = polyval(fliplr(calRef.C(1,:)), E_ref);
    belowRef = E_ref < calRef.E_floor(1);
    U_ref(belowRef) = (calRef.U_floor(1) / calRef.E_floor(1)) * E_ref(belowRef);
    aboveRef = E_ref > calRef.E_ceil(1);
    U_ref(aboveRef) = calRef.U_ceil(1);
    fprintf('Reference probe (T29): %.3f..%.3f m/s (%.1f%% below floor, %.1f%% above ceiling)\n', ...
        min(U_ref), max(U_ref), 100*nnz(belowRef)/numel(E_ref), 100*nnz(aboveRef)/numel(E_ref));
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

% What gets written if you choose to save (prompt is at the end, after the
% plots, so you can inspect the run before deciding).
% hwT_wind is the one to use for analysis: 0 = wind start. aiConfig/aiNames
% record which channels were actually logged in this run.
saveName = sprintf('hotwire_%s.mat', datestr(now,'yyyymmdd_HHMMSS'));
saveVars = {'data', 'hwT', 'hwT_t0', 'hwT_wind', 'CAL_FILE', ...
            'E1', 'E2', 'E3', 'U', 'V', 'W', 'E_ref', 'U_ref', ...
            'aiConfig', 'aiNames', 'aiGroups', 'sentT', 'sentV', ...
            'fanStartElapsed', 'hwStartElapsed', 'hwStopElapsed', ...
            'camStartElapsed', 'camStopElapsed', 'camDelayActual', 'camStartJitter'};

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
nSub = 2*double(hasHotwire) + 2*double(hasRef);
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

        iSub = iSub + 1;
        subplot(nSub,1,iSub)
        plot(hwT_wind, U_ref, 'LineWidth', 1);
        xlabel('Time since wind start (s)');
        ylabel('Velocity (m/s)');
        legend('U', 'V', 'W');
        title('Converted velocity (probe coordinates)');
        grid on;
    end
end

%% Save (prompted) -- after the plots so the run can be inspected first.
% Enter defaults to YES: an accidental keypress should not discard a run.
drawnow;   % make sure the figures are rendered before the prompt blocks
resp = input(sprintf('\nSave this run to %s? [Y/n]: ', saveName), 's');
if isempty(resp) || strncmpi(strtrim(resp), 'y', 1)
    save(saveName, saveVars{:});
    fprintf('Saved %s\n', saveName);
else
    disp('Not saved. All variables are still in the workspace -- to save later, run:');
    fprintf('  save(saveName, saveVars{:})\n');
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
