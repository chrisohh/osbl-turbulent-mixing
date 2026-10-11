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
%   - Fan on Dev4/ao0, hot-wire on Dev4/ai0:ai2, reference probe on Dev4/ai4,
%     wave gauge on Dev4/ai7 (SingleEnded) -- all selected via ENABLED_AI
%   - Fan ramp: 1.5->8V over 60s, hold 30s, ramp down 5s -- adjust as needed
%   - Hot-wire: 4 kHz sample rate
%   - Camera counters all share CAM_FS (default 50 Hz): IR/ctr2 (exposure
%     arbitrary -- IR sets its own), Color/ctr1 (500us) and Mono/ctr3 (700us),
%     set via CAM_EXPOSURE_US, since these trigger on HIGH level so pulse
%     width = exposure time
%   - ctr3 feeds BOTH the Mono side camera and the Mono under-tank camera
%
% NOT verified -- confirm before relying on this:
%   - Whether start(camD, "continuous") or plain start(camD) is correct for
%     this DAQ Toolbox version (see test_camera_triggers.m notes)

% Batch mode: run_experiment_batch.m calls this script with a BATCH_RUN struct
% (.name, .overrides) in the workspace -- keep it through the clear. Run by
% hand, BATCH_RUN is empty and everything behaves as before (prompts etc.).
if exist('BATCH_RUN', 'var') && isstruct(BATCH_RUN)
    clearvars -except BATCH_RUN
else
    clear; clc;
    BATCH_RUN = [];
end

DEV_ID = "Dev4";   % <-- confirm this matches NI MAX for the USB-6451
RUN_TAG = '';   % <-- optional short label for the saved filename, e.g. 'ease5_hold10'

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
% Cut the ramp early (s since wind start) WITHOUT changing its shape or rate:
% the run is identical to a full one up to this time, then holds/ramps down
% from the voltage reached. Inf = run the full ramp. Use this (not a lower
% FAN_V_END) to stop before waves reach the hot-wire -- lowering FAN_V_END
% alone would slow the whole ramp and change the forcing history.
FAN_RAMP_STOP_T = Inf;
FAN_RAMPDOWN_T  = 5;   % seconds, ramp up duration
FAN_HOLD_T  = 0;    % seconds, hold at V_END
FAN_DT      = 0.5;  % seconds, step interval

HOTWIRE_FS       = 4000;   % Hz
DELAY_BEFORE_TRIG = 10;    % seconds after wind start before starting camera counters
PRE_WIND_BASELINE_T = 0;   % seconds of zero-flow hot-wire data before the fan startsn (0 to skip)
% With the wave gauge enabled this is also the still-water window eta is zeroed
% against, which is what removes the gauge's offset drift (measured at
% 0.007 cm/min, and 0.2 cm between calibrations 1.5 h apart -- both far larger
% than anything else in the eta error budget). A wave gauge run therefore wants
% a non-zero baseline, so if ENABLED_AI includes WaveGauge and this is still 0,
% it is raised to WG_BASELINE_MIN_T below. An explicit batch override is never
% overridden, even to 0.
WG_BASELINE_MIN_T   = 10;  % s, minimum still-water window for zeroing the gauge
WG_ZERO_FROM_BASELINE = true;   % false = keep eta on the calibration's own datum

% ---- Channel + camera selection ----------------------------------------
% These live ABOVE the batch-override block on purpose: a wrapper script
% (run_surface_clean.m) needs to change the camera rate and the logged
% channels, and only settings defined above that block are overridable.
% What each one may contain is documented at the config it feeds, below.
ENABLED_AI      = {'Hotwire'};   % 'Hotwire', 'RefProbe', 'WaveGauge'
DEV_CAMERAS     = "Dev3";        % <-- confirm this matches NI MAX for the PCI-6621
ENABLED_CAMERAS = {'IR','Mono'}; % {'IR','Color','Mono'}, or {} for none
CAM_FS          = 50;            % Hz, trigger rate shared by every enabled camera

% Exposure per camera in MICROSECONDS -- not duty cycle, deliberately.
% setup_camera_triggers.m sets ch.DutyCycle (the PulseWidth path its own header
% describes is NOT implemented), and these cameras expose for as long as the
% line is HIGH, so exposure = dutyCycle/freq. A fixed duty cycle would
% therefore stretch every exposure in proportion as CAM_FS drops -- at 5 Hz the
% old 0.025 duty would be 5ms, not 500us. Giving exposure here and deriving
% duty where allCamConfig is built keeps exposure fixed when CAM_FS changes,
% which is what a low-rate run needs.
% IR ignores pulse width entirely (exposure set in its own software).
% Mono is 700us, matching the 0.035 duty this file's header always documented;
% the 0.95 duty (19ms) previously hard-coded here was for laser work and is not
% wanted for camera runs (user, 2026-10-05).
CAM_EXPOSURE_US = struct('IR', 10000, 'Color', 500, 'Mono', 19000);

% Batch overrides (e.g. FAN_RAMP_STOP_T, PRE_WIND_BASELINE_T) replace the
% values above for this run only. Only settings defined ABOVE this point can
% be overridden -- a misspelt name errors instead of being silently ignored.
% batchOverrideNames is kept (not cleared) so later defaults can tell an
% explicitly overridden setting from one still at its scripted value.
batchOverrideNames = {};
if ~isempty(BATCH_RUN) && isfield(BATCH_RUN, 'overrides')
    ovNames = fieldnames(BATCH_RUN.overrides);
    for iOv = 1:numel(ovNames)
        if ~exist(ovNames{iOv}, 'var')
            error('Batch override "%s" is not a setting defined above this point in run_experiment.m.', ovNames{iOv});
        end
        eval([ovNames{iOv} ' = BATCH_RUN.overrides.(ovNames{iOv});']);
        fprintf('Batch override: %s = %s\n', ovNames{iOv}, fmt_override(BATCH_RUN.overrides.(ovNames{iOv})));
    end
    batchOverrideNames = ovNames;
    clear ovNames iOv
end

% Analog input channels on the USB-6451. Hot-wire probes are differential
% (ai0/ai1/ai2, hardware-paired with ai8/ai9/ai10). RefProbe is the velocity
% reference transducer.
% CONFIRM: RefProbe is SingleEnded here because only one pin (ai4) was given.
% If it is actually wired differentially it uses ai4/ai12 and 'term' must be
% "Differential" instead.
% WaveGauge is the resistance wave gauge on ai7, SingleEnded -- same channel
% and terminal config as WaveGauge_receiver.m and Wave_Gauge/
% run_wavegauge_fan_calibration.m use. It rides along in this same AI task, so
% eta shares one clock with the hot-wire and the camera triggers; no second
% session, no cross-session drift.
% NOT VERIFIED: enabling WaveGauge alongside Hotwire/RefProbe makes this a
% MIXED-terminal task (SingleEnded + Differential) for the first time --
% setup_hotwire_daq.m's header already flags that as unconfirmed on this
% device. If addinput errors, either run the gauge on its own or rewire it
% differentially (ai7 pairs with ai15) and change 'term' below.
% Channels are selected by GROUP, not individually: Probe1-3 are the three
% sensors of one tri-axial probe and are always logged together (the U,V,W
% decomposition needs all three), so they cannot be split apart here.
allAiConfig = struct( ...
    'group', {'Hotwire',      'Hotwire',      'Hotwire',      'RefProbe',     'WaveGauge'}, ...
    'name',  {'Probe1',       'Probe2',       'Probe3',       'RefProbe',     'WaveGauge'}, ...
    'chan',  {"ai0",          "ai1",          "ai2",          "ai4",          "ai7"}, ...
    'term',  {"Differential", "Differential", "Differential", "Differential", "SingleEnded"}, ...
    'range', {[-10 10],       [-10 10],       [-10 10],       [-10 10],       [-10 10]});

% Valid groups: 'Hotwire', 'RefProbe', 'WaveGauge'. Any combination -- e.g.
% {'RefProbe'} for reference only, {'Hotwire','WaveGauge'} to log waves
% alongside the probe, {'WaveGauge'} for the gauge alone.
% NOTE: every channel in the task shares HOTWIRE_FS. At 4 kHz the gauge is
% heavily oversampled (waves are a few Hz) -- harmless, just a bigger file.
% ENABLED_AI itself is set above the batch-override block (so a wrapper can
% override it); the groups it may contain are the ones documented just above.
aiConfig   = allAiConfig(ismember({allAiConfig.group}, ENABLED_AI));
aiNames    = {aiConfig.name};
aiGroups   = {aiConfig.group};
hasHotwire = any(strcmp(aiGroups, 'Hotwire'));
hasRef     = any(strcmp(aiGroups, 'RefProbe'));
hasWave    = any(strcmp(aiGroups, 'WaveGauge'));

if isempty(aiConfig)
    error('ENABLED_AI selected no channels -- valid groups are ''Hotwire'', ''RefProbe'' and ''WaveGauge''.');
end
badGroups = setdiff(ENABLED_AI, {allAiConfig.group});
if ~isempty(badGroups)
    error('ENABLED_AI contains unknown group(s): %s', strjoin(badGroups, ', '));
end

% A wave gauge run with no still-water window has nothing to zero eta against,
% so the gauge's offset drift stays in the data. Raise the baseline rather than
% let that happen silently -- but never argue with an explicit batch override.
if hasWave && PRE_WIND_BASELINE_T < WG_BASELINE_MIN_T
    if ismember('PRE_WIND_BASELINE_T', batchOverrideNames)
        warning(['WaveGauge is enabled but PRE_WIND_BASELINE_T was overridden to %g s ' ...
                 '(< %g s): eta will keep the calibration''s own zero, including its drift.'], ...
                 PRE_WIND_BASELINE_T, WG_BASELINE_MIN_T);
    else
        fprintf('WaveGauge enabled: raising PRE_WIND_BASELINE_T %g -> %g s for still-water zeroing.\n', ...
            PRE_WIND_BASELINE_T, WG_BASELINE_MIN_T);
        PRE_WIND_BASELINE_T = WG_BASELINE_MIN_T;
    end
end

% DEV_CAMERAS, ENABLED_CAMERAS, CAM_FS and CAM_EXPOSURE_US are all set above
% the batch-override block, so a wrapper script can change them.
% Pulse width IS the exposure: these cameras trigger on HIGH level, so they
% expose for as long as the line is high -- so exposure = dutyCycle/freq, and
% dutyCycle is DERIVED from CAM_EXPOSURE_US below rather than hard-coded, so
% that changing CAM_FS does not silently change every exposure with it.
% IR ignores pulse width (its exposure is set in its own software), so its
% 10000us entry is arbitrary.
%
% ctr3 DRIVES TWO CAMERAS: the Mono side camera and the Mono under-tank
% cross-view camera are wired to the same counter output, so they necessarily
% share this signal -- same rate, same 700us exposure. Changing 'Mono' below
% changes both. To give the under-tank camera its own exposure it would need
% a separate counter (ctr0 is free) and its own entry here.
% ctr1 (Color) and ctr2 (IR) each drive a single camera.
allCamConfig = struct( ...
    'name',      {'IR',   'Color',  'Mono'}, ...
    'ctr',       {'ctr2', 'ctr1', 'ctr3'}, ...
    'freq',      {CAM_FS, CAM_FS, CAM_FS}, ...
    'dutyCycle', {NaN,    NaN,   NaN}, ...   % derived from CAM_EXPOSURE_US below
    'delay',     {0,      0,    0});

% 'Mono' here means the ctr3 line, i.e. Mono side + Mono under-tank.
% duty = exposure * rate. This errors rather than clipping: duty >= 1 means the
% requested exposure is longer than the frame period, which would leave the
% trigger line permanently high -- no usable edge at all.
for iCam = 1:numel(allCamConfig)
    camName = allCamConfig(iCam).name;
    if ~isfield(CAM_EXPOSURE_US, camName)
        error('CAM_EXPOSURE_US has no field "%s" -- every camera in allCamConfig needs one.', camName);
    end
    expUs = CAM_EXPOSURE_US.(camName);
    duty  = expUs * CAM_FS / 1e6;
    if ~(duty > 0 && duty < 1)
        error(['%s: exposure %g us at %g Hz gives duty %.4f, which is not in (0,1). ' ...
               'The frame period is %g us -- lower CAM_EXPOSURE_US.%s or CAM_FS.'], ...
               camName, expUs, CAM_FS, duty, 1e6/CAM_FS, camName);
    end
    allCamConfig(iCam).dutyCycle = duty;
end
clear iCam camName expUs duty

% A typo in ENABLED_CAMERAS used to silently trigger nothing. Now that a
% wrapper can override it, catch it -- same treatment as ENABLED_AI gets above.
badCams = setdiff(ENABLED_CAMERAS, {allCamConfig.name});
if ~isempty(badCams)
    error('ENABLED_CAMERAS contains unknown camera(s): %s', strjoin(badCams, ', '));
end

camConfig = allCamConfig(ismember({allCamConfig.name}, ENABLED_CAMERAS));

% probe4.txt holds BOTH probes used below: the tri-axial hot-wire (55P95,
% read by parse_calibration.m/convert_E2U_fn.m, which take the file's FIRST
% probe block) and the reference probe (T29, read by parse_probe_section.m
% by name further down).
CAL_FILE = 'D:\Chris\osbl-turbulent-mixing\experiment\data\260909\probe4.txt'; % <-- this probe's cal/header
CTA_DIR  = 'D:\Chris\osbl-turbulent-mixing\experiment\processing_code\CTA';     % parse_calibration + convert_E2U_fn + parse_probe_section
addpath(CTA_DIR);

% inferno lives with the CSIG code -- same lookup as
% Wave_Gauge/plotting/plot_wave_ramp_run.m, used by the eta spectrogram below.
SLOPE_UTIL = 'D:\Chris\osbl-turbulent-mixing\experiment\processing_code\Slope_Gauge\util';
if ~exist('inferno', 'file') && exist(SLOPE_UTIL, 'dir'), addpath(SLOPE_UTIL); end

% Wave gauge volts -> cm, from a STATIC calibration (still water, fan off)
% made by Wave_Gauge/run_wavegauge_calibration.m. Empty = use the newest
% wave_gauge_calibration_*.mat found in the Wave_Gauge folder, here, or the
% current folder. Without one, eta stays in VOLTS -- the run is still valid,
% it just isn't in physical units. Same lookup as
% Wave_Gauge/run_wavegauge_fan_calibration.m.
WG_DIR      = 'D:\Chris\osbl-turbulent-mixing\experiment\processing_code\Wave_Gauge';
WG_CAL_FILE = '';
wgCal = [];
if hasWave
    if isempty(WG_CAL_FILE)
        wgSearch = unique({WG_DIR, fileparts(mfilename('fullpath')), pwd});
        fWg = [];
        for iWg = 1:numel(wgSearch)
            fWg = [fWg; dir(fullfile(wgSearch{iWg}, 'wave_gauge_calibration_*.mat'))]; %#ok<AGROW>
        end
        if ~isempty(fWg)
            [~, iNewest] = max([fWg.datenum]);
            WG_CAL_FILE = fullfile(fWg(iNewest).folder, fWg(iNewest).name);
        end
        clear wgSearch fWg iWg iNewest
    end
    if ~isempty(WG_CAL_FILE)
        Swg = load(WG_CAL_FILE, 'calibration');
        wgCal = Swg.calibration;
        fprintf('Wave gauge calibration: %s\n', WG_CAL_FILE);
        clear Swg
    else
        warning('No wave gauge calibration found -- eta will be in VOLTS.');
    end
end

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
    holdDt    = FAN_HOLD_T / nHoldSteps;   % holdEndT is set after the ramp is built (it depends on FAN_RAMP_STOP_T)

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

    % Early stop: cut the SAME ramp at FAN_RAMP_STOP_T instead of reshaping it,
    % so everything up to the cut matches a full-length run exactly (same
    % rate, same forcing history) -- then hold/ramp down from wherever it got.
    tUp = FAN_V_START_HOLD_T + (0:nStepsUp-1) * FAN_DT;   % wind clock of each step
    keepUp = tUp <= FAN_RAMP_STOP_T + 1e-9;
    rampUpV = rampUpV(keepUp);
    nStepsUp = numel(rampUpV);
    rampEndT = FAN_V_START_HOLD_T + (nStepsUp-1) * FAN_DT;   % = FAN_RAMPUP_T if not cut
    V_top = rampUpV(end);
    holdEndT = rampEndT + FAN_HOLD_T;
    if rampEndT < FAN_RAMPUP_T
        if FAN_RAMP_LINEAR_U
            fprintf('Ramp CUT at t=%.1fs (of %.1fs): stops at %.2f V (~%.2f m/s), same rate as the full ramp.\n', ...
                rampEndT, FAN_RAMPUP_T, V_top, interp1(tf.V, tf.U, V_top, 'pchip'));
        else
            fprintf('Ramp CUT at t=%.1fs (of %.1fs): stops at %.2f V.\n', rampEndT, FAN_RAMPUP_T, V_top);
        end
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

    fprintf('Holding at %.2f V...\n', V_top);
    for i = 1:nHoldSteps
        [camerasStarted, camStartElapsed] = waitUntilWithCameraCheck(fanStartElapsed + rampEndT + i*holdDt, camerasStarted, camStartElapsed, t0, camTargetElapsed, camD);
        sentT(end+1) = toc(t0); sentV(end+1) = V_top; %#ok<SAGROW>
    end

    disp('Ramping down...');
    % Down to 0V, not just back to FAN_V_START -- the fan should start AND
    % end the run at 0, not be left spinning at FAN_V_START afterward. Below
    % ~2.0V nothing sustains flow anyway (dead band), so this last stretch of
    % the ramp is just as physically discontinuous on the way down as the
    % 0->FAN_V_START jump was on the way up.
    rampDownV = linspace(V_top, 0, nStepsDown);
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

%% Wave gauge: volts -> surface elevation (same pipeline as the Wave_Gauge scripts)
E_wg = [];
eta = [];
etaUnit = '';
etaZeroOffset = 0;
if hasWave
    E_wg = data{:, strcmp(aiNames, 'WaveGauge')};
    if ~isempty(wgCal)
        % Cubic fit from the static calibration; outside the calibrated
        % voltage span this is extrapolation off the end of that cubic, which
        % diverges fast -- so the share is reported, not silently accepted.
        eta = polyval(wgCal.polynomial_coefficients, E_wg);
        etaUnit = 'cm';
        vc = wgCal.voltage_calib;
        if isfield(wgCal, 'clipped'), vc = vc(~wgCal.clipped); end   % older cal files lack it
        outOfCal = E_wg < min(vc) | E_wg > max(vc);
        fprintf('Wave gauge: eta %.2f..%.2f cm (%.1f%% of samples outside the calibrated voltage range -- extrapolated).\n', ...
            min(eta), max(eta), 100*nnz(outOfCal)/numel(E_wg));
    else
        eta = E_wg;
        etaUnit = 'V';
        fprintf('Wave gauge: %.3f..%.3f V (no calibration -- eta is in volts).\n', min(E_wg), max(E_wg));
    end
    wgRange = aiConfig(strcmp(aiNames, 'WaveGauge')).range;
    if any(E_wg <= wgRange(1) + 0.05 | E_wg >= wgRange(2) - 0.05)
        warning('%.2f%% of wave gauge samples are at the +-%g V input limit.', ...
            100*nnz(E_wg <= wgRange(1) + 0.05 | E_wg >= wgRange(2) - 0.05)/numel(E_wg), wgRange(2));
    end

    % Still-water zeroing. The gauge's offset is the one error term big enough
    % to matter (0.2 cm between calibrations 1.5 h apart, drifting at
    % 0.007 cm/min) and it cancels exactly against a pre-wind mean. The shift
    % is SAVED, not just applied, so absolute elevation on the calibration's
    % own datum is recoverable as eta + etaZeroOffset.
    etaZeroOffset = 0;
    baseIdx = hwT_wind < 0;
    if WG_ZERO_FROM_BASELINE && nnz(baseIdx) > 0
        etaZeroOffset = mean(eta(baseIdx));
        eta = eta - etaZeroOffset;
        fprintf('Wave gauge zeroed on %.1f s of still water (%d samples): offset %+.4f %s removed.\n', ...
            hwT_wind(find(baseIdx, 1, 'last')) - hwT_wind(1), nnz(baseIdx), etaZeroOffset, etaUnit);
        % Baseline scatter is the floor on how well that zero is known.
        fprintf('  baseline std %.4f %s (that is the uncertainty on the new zero).\n', ...
            std(eta(baseIdx)), etaUnit);
    elseif WG_ZERO_FROM_BASELINE
        warning(['No pre-wind samples, so eta keeps the calibration''s own zero ' ...
                 'and carries its offset drift. Set PRE_WIND_BASELINE_T > 0.']);
    end
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
% RUN_TAG makes runs distinguishable by filename alone (e.g. "ease5_hold10")
% instead of only by timestamp -- set it above with the other config, or
% leave '' for the old timestamp-only behavior.
if isempty(RUN_TAG)
    saveName = sprintf('hotwire_%s.mat', datestr(now,'yyyymmdd_HHMMSS'));
else
    saveName = sprintf('hotwire_%s_%s.mat', RUN_TAG, datestr(now,'yyyymmdd_HHMMSS'));
end
% Every setting that shapes the run, so a saved file says how it was made
% (earlier files didn't, and the ramp had to be reverse-engineered from
% sentV). runConfig.fanTF is a copy of the transfer function actually used,
% since fan_transfer.mat may be rebuilt later.
runConfig = struct( ...
    'RUN_TAG', RUN_TAG, ...
    'FAN_V_START', FAN_V_START, 'FAN_V_END', FAN_V_END, ...
    'FAN_RAMP_LINEAR_U', FAN_RAMP_LINEAR_U, 'FAN_TF_FILE', FAN_TF_FILE, ...
    'FAN_V_START_HOLD_T', FAN_V_START_HOLD_T, 'FAN_RAMP_EASE_T', FAN_RAMP_EASE_T, ...
    'FAN_RAMPUP_T', FAN_RAMPUP_T, 'FAN_RAMP_STOP_T', FAN_RAMP_STOP_T, ...
    'FAN_RAMPDOWN_T', FAN_RAMPDOWN_T, 'FAN_HOLD_T', FAN_HOLD_T, 'FAN_DT', FAN_DT, ...
    'rampEndT', rampEndT, 'V_top', V_top, ...           % where the ramp actually stopped
    'HOTWIRE_FS', HOTWIRE_FS, 'DELAY_BEFORE_TRIG', DELAY_BEFORE_TRIG, ...
    'PRE_WIND_BASELINE_T', PRE_WIND_BASELINE_T, ...
    'ENABLED_AI', {ENABLED_AI}, 'ENABLED_CAMERAS', {ENABLED_CAMERAS}, ...
    'CAM_FS', CAM_FS, 'CAM_EXPOSURE_US', CAM_EXPOSURE_US, ...
    'WG_CAL_FILE', WG_CAL_FILE, 'WG_ZERO_FROM_BASELINE', WG_ZERO_FROM_BASELINE, ...
    'fanTF', []);
if FAN_RAMP_LINEAR_U, runConfig.fanTF = tf; end
saveVars = {'data', 'hwT', 'hwT_t0', 'hwT_wind', 'CAL_FILE', ...
            'E1', 'E2', 'E3', 'U', 'V', 'W', 'E_ref', 'U_ref', ...
            'E_wg', 'eta', 'etaUnit', 'etaZeroOffset', 'wgCal', 'WG_CAL_FILE', ...
            'aiConfig', 'aiNames', 'aiGroups', 'sentT', 'sentV', ...
            'fanStartElapsed', 'hwStartElapsed', 'hwStopElapsed', ...
            'camStartElapsed', 'camStopElapsed', 'camDelayActual', 'camStartJitter', ...
            'runConfig'};

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
% (No hwT_wind=hwT here: that overwrote the wind clock before the save, so
% saved files had hwT_wind ~0.1 s off, and it would hide any pre-wind
% baseline. hwT_wind is already 0 = wind start from the read-back above.)
% The wave gauge contributes a spectrogram row plus an eta row, and a raw-volts
% row when calibrated (uncalibrated, eta IS the voltage, so that row would just
% repeat the eta row). Raw volts get their own row rather than a second y axis
% on the eta row: the calibration is a cubic, so volts and cm are not a linear
% rescale of each other and a twin axis would imply they are.
nSub = 2*double(hasHotwire) + 2*double(hasRef) + double(hasWave)*(2 + double(~isempty(wgCal)));
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

    if hasWave
        if ~isempty(wgCal)
            iSub = iSub + 1;
            subplot(nSub,1,iSub)
            plot(hwT_wind, E_wg, 'LineWidth', 0.5);
            % Dashed = ends of the calibrated voltage span; beyond them eta is
            % extrapolated off the cubic (percentage printed above).
            yline(min(vc), 'r--'); yline(max(vc), 'r--');
            xlabel('Time since wind start (s)');
            ylabel('Voltage (V)');
            legend('E_{wg}');
            title('Wave gauge raw voltage (ai7; red dashed = calibrated range)');
            grid on;
        end

        iSub = iSub + 1;
        subplot(nSub,1,iSub)
        plot(hwT_wind, eta, 'LineWidth', 0.5);
        xlabel('Time since wind start (s)');
        ylabel(sprintf('\\eta (%s)', etaUnit));
        legend('\eta');
        title(sprintf('Wave gauge surface elevation (%s)', etaUnit));
        grid on;

        % omega-vs-t spectrogram, same form as the tile in
        % Wave_Gauge/plotting/plot_wave_ramp_run.m: angular frequency on y,
        % log10 of the one-sided PSD, inferno. Shows the peak marching down in
        % frequency as the waves grow through the ramp.
        iSub = iSub + 1;
        axSpec = subplot(nSub,1,iSub);
        WG_SPEC_WIN    = 4;    % s, window length
        WG_SPEC_F_MAX  = 10;   % Hz, top of the plotted band
        % Colour range, in DECADES of S_eta below the loudest bin in the band --
        % same units and default as Wave_Gauge/plotting/plot_wave_ramp_run.m.
        % Without it the scale stretches to the pre-wind noise floor ~11 decades
        % down, which compresses the wave field into the top few colours and
        % promotes the constant ~7.2 Hz instrumental line into a feature.
        WG_SPEC_DYN_DECADES = 6;
        nWin = round(WG_SPEC_WIN * HOTWIRE_FS);
        % 4th output is the one-sided PSD (etaUnit^2/Hz); the 1st output |s|^2
        % would be raw FFT power with no physical units.
        [~, fS, tS, pS] = spectrogram(detrend(eta), hann(nWin), round(nWin/2), [], HOTWIRE_FS);
        % S(omega) = S(f)/(2 pi) keeps the variance under the curve unchanged
        % when the abscissa is stretched by 2 pi.
        specLog = log10(pS/(2*pi) + eps);
        imagesc(tS + hwT_wind(1), 2*pi*fS, specLog);
        set(gca, 'YDir', 'normal');
        colormap(axSpec, pick_wg_cmap());
        % Max over the plotted band only, so an out-of-view high-frequency bin
        % cannot set the scale for the band you are looking at.
        topLog = max(max(specLog(fS <= WG_SPEC_F_MAX, :)));
        clim([topLog - WG_SPEC_DYN_DECADES, topLog]);
        cb = colorbar;
        cb.Label.String = sprintf('log_{10} S_\\eta(\\omega)  (%s^2 s/rad)', etaUnit);
        ylim([0 2*pi*WG_SPEC_F_MAX]);
        xlabel('Time since wind start (s)');
        ylabel('\omega (rad/s)');
        title(sprintf('S_\\eta(\\omega, t) (%g s windows)', WG_SPEC_WIN));
    end

    % The spectrogram's colorbar steals width from its axes, so without this
    % its time axis sits narrower than the line plots above it and the columns
    % stop lining up. Same fix as plot_wave_ramp_run.m: give every other axes
    % an invisible colorbar, then link the time axes so zooming one zooms all.
    if hasWave
        axAll = flipud(findobj(gcf, 'Type', 'axes'));
        for iA = 1:numel(axAll)
            if isempty(axAll(iA).Colorbar)
                cbHidden = colorbar(axAll(iA)); cbHidden.Visible = 'off';
            end
        end
        linkaxes(axAll, 'x');
        xlim(axAll(1), [hwT_wind(1), hwT_wind(end)]);
        clear axAll iA cbHidden
    end
end

%% Save (prompted) -- after the plots so the run can be inspected first.
% Naming happens HERE, after you've seen the plots -- usually you know what
% a run was worth calling only once you've looked at it, not before it ran.
% Blank keeps whatever RUN_TAG produced (timestamp-only if RUN_TAG was '').
% Enter defaults to YES on the save itself: an accidental keypress should
% not discard a run.
drawnow;   % make sure the figures are rendered before the prompt blocks
if ~isempty(BATCH_RUN)
    % Batch mode (run_experiment_batch.m): no prompts, name comes from the batch.
    saveName = next_repeat_name(regexprep(strtrim(BATCH_RUN.name), '[^\w-]', '_'));
    resp = 'y';
else
    tagResp = input(sprintf('\nName for this run (blank = keep "%s"): ', saveName), 's');
    tagResp = strtrim(tagResp);
    if ~isempty(tagResp)
        % Repeats get the same typed name -- never overwrite an earlier run
        % (save() does so silently); they're numbered _r1, _r2, ... instead.
        saveName = next_repeat_name(regexprep(tagResp, '[^\w-]', '_'));
    end
    resp = input(sprintf('Save this run to %s? [Y/n]: ', saveName), 's');
end
if isempty(resp) || strncmpi(strtrim(resp), 'y', 1)
    save(saveName, saveVars{:});
    fprintf('Saved %s\n', saveName);
else
    disp('Not saved. All variables are still in the workspace -- to save later, run:');
    fprintf('  save(saveName, saveVars{:})\n');
end


%% Local functions
function cm = pick_wg_cmap()
% inferno when Slope_Gauge/util is on the path (added near CAL_FILE above),
% parula otherwise -- same fallback as the Wave_Gauge plotting functions, so a
% missing path costs you the colormap, not the figure.
    if exist('inferno', 'file'), cm = inferno; else, cm = parula; end
end

function name = next_repeat_name(base)
% Next free repeat file for a run name: base_r1.mat, base_r2.mat, ... in the
% current folder (where save() writes). A plain base.mat from before this
% numbering counts as repeat 1, so nothing is ever overwritten.
% If the name already ends in _r<n> (typed by hand), that number is used as
% given -- unless that file exists, in which case the next free number is.
    tok = regexp(base, '^(.*)_r(\d+)$', 'tokens', 'once');
    if ~isempty(tok)
        if ~isfile([base '.mat'])
            name = [base '.mat'];
            return
        end
        fprintf('%s.mat already exists -- using the next free repeat number instead.\n', base);
        base = tok{1};
    end
    used = [];
    if isfile([base '.mat']), used(end+1) = 1; end
    d = dir([base '_r*.mat']);
    for k = 1:numel(d)
        tok = regexp(d(k).name, ['^' regexptranslate('escape', base) '_r(\d+)\.mat$'], 'tokens', 'once');
        if ~isempty(tok), used(end+1) = str2double(tok{1}); end %#ok<AGROW>
    end
    if isempty(used), n = 1; else, n = max(used) + 1; end
    name = sprintf('%s_r%d.mat', base, n);
end

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

function s = fmt_override(v)
% Render an override value for the "Batch override:" log line. mat2str covers
% numeric, logical and char, but ERRORS on the cell arrays (ENABLED_AI,
% ENABLED_CAMERAS) and structs (CAM_EXPOSURE_US) that became overridable when
% they moved above the override block -- and an error here would kill the run
% while merely trying to log what it was about to do.
if iscell(v)
    if isempty(v)
        s = '{}';
    elseif all(cellfun(@(x) ischar(x) || isstring(x), v))
        quoted = cellfun(@(x) ['''' char(x) ''''], v, 'UniformOutput', false);
        s = ['{' strjoin(quoted, ', ') '}'];
    else
        s = sprintf('<%dx%d cell>', size(v,1), size(v,2));
    end
elseif isstruct(v) && isscalar(v)
    f = fieldnames(v);
    parts = cell(1, numel(f));
    for i = 1:numel(f)
        parts{i} = sprintf('%s=%s', f{i}, fmt_override(v.(f{i})));
    end
    s = ['struct(' strjoin(parts, ', ') ')'];
elseif isstring(v)
    s = char(strjoin(v, ', '));
else
    s = mat2str(v);
end
end
