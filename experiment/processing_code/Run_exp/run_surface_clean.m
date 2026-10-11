function run_surface_clean()
%% run_surface_clean.m
% Surface-clean ("dust-off") run: a long, LOW-wind blow that sweeps the
% surfactant film / floating dust downwind so it can be skimmed, while every
% instrument logs and the cameras run at a reduced frame rate so the run costs
% little disk space. Thin wrapper -- it only overrides settings in
% run_experiment.m, so all the DAQ logic has exactly one home.
%
% WHY LOW WIND: the point is to move the surface film, not to force the water.
% At ~1.85 m/s the interfacial Richardson number for a 10 cm fresh layer is
% O(10^3), i.e. two orders of magnitude above the ~10-30 where entrainment
% starts, so the blow costs ~0.2 mm of layer depth over 3 minutes. The erosion
% estimate scales as U^3 (see WHAT THIS COSTS below), which is why this script
% refuses to run at a voltage high enough to matter.
%
% WHAT THIS DOES NOT DO: wind alone does not REMOVE anything -- it sweeps the
% film to the downwind end, where a surfactant film re-spreads by Marangoni
% within seconds to a minute of the fan stopping. The blow is step 1 of 3:
%   1. this script (sweep the film downwind, instruments logging)
%   2. fan off, then IMMEDIATELY skim the downwind accumulation line
%      (syringe/pipette at the surface, or lay and lift clean lens tissue)
%   3. CT profile, then your actual wind run
% Steps 2 and 3 are manual. Skipping step 2 largely wastes step 1.
%
% BONUS USE: with the fan running and (at this speed) essentially no waves,
% the wave gauge record here is your noise floor WITH THE FAN ON -- fan
% vibration, rig coupling and electrical pickup on ai7, none of which appear
% in a still-water PRE_WIND_BASELINE. You need that floor to tell "no waves"
% from "no waves plus 0.03 cm of fan noise" when you later judge whether the
% surface came clean.
%
% SAFETY: Ctrl+C at any point sets the fan to 0 V -- this is a function, so
% its onCleanup fires on Ctrl+C, unlike a plain script. If anything still
% looks stuck after an interrupt, run cleanup_daq.m anyway.
%
% Edit the settings below, then press Run (F5).

%% ---- Settings ------------------------------------------------------------
SAVE_NAME  = 'surfclean';   % goes into the filename as hotwire_<SAVE_NAME>_<timestamp>.mat

WIND_V     = 0;    % V. ~1.85 m/s (2.1 V -> 1.82-1.85 m/s across the
                     % 2026-09-22/23 calibrations). Do NOT use 2.0: the fan
                     % cycles on/off at 2.0-2.01 V and the flow is unsteady
                     % there. The V->U curve is nearly flat coming out of the
                     % dead band, so 2.1 V buys steady flow for ~2% more wind.
BLOW_T     = 240;    % s held at WIND_V. Surface drift is ~3% of the wind, so
                     % ~5.5 cm/s here -- one pass down a ~15 m tank takes
                     % ~270 s, which is why this is 300 and not 180: below one
                     % full pass, film upwind of the measurement region never
                     % reaches the downwind end to be skimmed. Costs ~0.9 mm of
                     % a 10 cm layer; 600 s (two passes) costs ~1.8 mm, also
                     % fine. The U^3 scaling is what makes a long blow cheap
                     % HERE and expensive at any higher voltage.
CAM_FS     = 10;     % Hz. 50 -> 10 cuts frames 5x (~2970 instead of ~14850 per
                     % camera over this run). Enough to SEE the surface state;
                     % it aliases capillary waves badly, so let the wave gauge
                     % carry anything quantitative.
HOTWIRE_FS = 1000;   % Hz. The .mat size is set by the AI task, not the
                     % cameras. 4 kHz is for turbulence spectra; a sweep run
                     % does not need it, and 1 kHz is still ~100x the wave
                     % frequencies the gauge sees.
AI_GROUPS  = {'Hotwire'};   % log everything
CAMERAS    = {'IR','Mono'};   % as run_experiment.m's default; {} for none

% Guard: this script exists for a gentle sweep. The U^3 scaling means a
% "quick" blow at a higher voltage quietly turns into a mixing run, so refuse
% rather than let that happen by a typo. Raise deliberately if you mean it.
WIND_V_MAX = 3.0;    % V (~2.97 m/s)

DEV_ID = "Dev4";     % fan device, only used to zero the fan on interrupt

%% ---- Checks --------------------------------------------------------------
if WIND_V < 2.05
    error(['WIND_V = %.2f V is in or below the fan''s on/off cycling band ' ...
           '(2.0-2.01 V) -- flow there is unsteady. Use 2.1 V or higher.'], WIND_V);
end
if WIND_V > WIND_V_MAX
    error(['WIND_V = %.2f V exceeds WIND_V_MAX = %.2f V. This script is for a ' ...
           'gentle sweep; erosion scales as U^3, so a higher blow is a mixing ' ...
           'run and belongs in run_experiment.m with a stated FAN_RAMP_STOP_T.'], ...
           WIND_V, WIND_V_MAX);
end

% FAN_RAMP_STOP_T is a TIME, not a voltage. Setting it equal to
% FAN_V_START_HOLD_T keeps exactly the first ramp step, so the fan holds at
% FAN_V_START for FAN_HOLD_T and the climb never begins. It must not go BELOW
% the hold: that empties rampUpV and run_experiment.m errors on rampUpV(1).
V_START_HOLD_T = 2;    % s, as in run_experiment.m
RAMP_STOP_T    = V_START_HOLD_T;

%% ---- Overrides -----------------------------------------------------------
% Only settings defined ABOVE run_experiment.m's batch-override block can be
% overridden; a misspelt name errors there rather than being ignored.
OVERRIDES = struct( ...
    'RUN_TAG',             SAVE_NAME, ...
    'FAN_V_START',         WIND_V, ...
    'FAN_V_START_HOLD_T',  V_START_HOLD_T, ...
    'FAN_RAMP_STOP_T',     RAMP_STOP_T, ...   % = hold -> flat at FAN_V_START
    'FAN_HOLD_T',          BLOW_T, ...
    'HOTWIRE_FS',          HOTWIRE_FS, ...
    'PRE_WIND_BASELINE_T', 10, ...            % still-water window to zero eta
    'ENABLED_AI',          {AI_GROUPS}, ...
    'ENABLED_CAMERAS',     {CAMERAS}, ...
    'CAM_FS',              CAM_FS);

%% ---- Run -----------------------------------------------------------------
% A manual run of run_experiment.m leaves hw/fanD/camD in the BASE workspace
% still holding Dev4, so this would fail with "Hardware is reserved". Zero the
% fan through them and release them first.
release_base_daq();

fanGuard = onCleanup(@() zero_fan(DEV_ID)); %#ok<NASGU>

fprintf('\n===== Surface-clean run =====\n');
fprintf('  %.2f V held for %g s, cameras %g Hz, AI %g Hz\n', WIND_V, BLOW_T, CAM_FS, HOTWIRE_FS);
fprintf('  Total run ~%.0f s (10 s baseline + %g s hold + 5 s ramp-down)\n', ...
    10 + RAMP_STOP_T + BLOW_T + 5, BLOW_T);
report_cost(WIND_V, BLOW_T);
fprintf('  AFTER THIS: fan off -> skim the downwind line -> CT profile.\n\n');

do_one_run(struct('name', SAVE_NAME, 'overrides', OVERRIDES));

fprintf('\n===== Done. Skim the downwind accumulation line NOW, before the film re-spreads. =====\n');
end


function report_cost(windV, blowT)
% Rough erosion estimate so the run states its own cost. Kato & Phillips
% entrainment, u_e/u* = 2.5/Ri*, integrated over the hold:
%   dh = 2.5 * u*^3 * t / (g' * h)
% with u*_water ~= 1.249e-3 * U (rho_a/rho_w, C_d = 1.3e-3).
% g' is NOT measured here -- 1% density contrast is assumed, and the answer
% scales as 1/g', so halve it at 2%. Treat as order-of-magnitude, not truth:
% the 2.5 coefficient comes from an annulus with no breaking waves.
try
    tf = load('fan_transfer.mat', 'fanTF'); tf = tf.fanTF;
    U = interp1(tf.V, tf.U, windV, 'pchip');
catch
    U = NaN;
end
if isnan(U)
    fprintf('  (could not read fan_transfer.mat -- no erosion estimate)\n');
    return
end
gPrime = 0.1;   % m/s^2, i.e. drho/rho = 1% -- ASSUMED
h0     = 0.10;  % m, fresh layer depth -- ASSUMED
uStar  = 1.249e-3 * U;
dh     = 2.5 * uStar^3 * blowT / (gPrime * h0);
fprintf('  ~%.2f m/s; estimated deepening ~%.2f mm of a %.0f cm layer (assumes drho/rho=1%%)\n', ...
    U, dh*1000, h0*100);
end


function do_one_run(BATCH_RUN) %#ok<INUSD> -- read by run_experiment.m
% Runs the script in THIS function's workspace, so its clearvars can't touch
% the caller's variables. Same trick run_experiment_batch.m uses.
run_experiment;
end


function zero_fan(DEV_ID)
% Two tries: a fresh session; if the output is still reserved by a leftover
% object, release everything (daqreset) and try again.
for attempt = 1:2
    try
        d = daq("ni");
        addoutput(d, DEV_ID, "ao0", "Voltage");
        write(d, 0);
        disp('Fan set to 0 V.');
        return
    catch ME
        if attempt == 1
            release_base_daq();
            daqreset;
        else
            fprintf(2, 'Could not zero the fan (%s). CHECK THE FAN PHYSICALLY.\n', ME.message);
        end
    end
end
end


function release_base_daq()
% Stop and clear DAQ objects a manual run left in the base workspace, writing
% 0 V through any fan object first (it's the one that can still drive ao0).
for nm = ["fanD", "fanCleanup", "hw", "camD"]
    if ~evalin('base', sprintf('exist(''%s'', ''var'')', nm)), continue; end
    obj = evalin('base', nm);
    if startsWith(nm, "fan")
        try, write(obj, 0); catch, end
    end
    try, stop(obj); catch, end
    evalin('base', sprintf('clear %s', nm));
    fprintf('Released leftover %s from the base workspace.\n', nm);
end
end
