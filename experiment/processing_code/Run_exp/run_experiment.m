%% run_experiment.m
% Top-level script: starts fan ramp + hot-wire at t=0 on the USB-6451,
% waits 30 s, then fires a single start pulse on DIO0 to arm the
% PCI-6621 camera counters (IR/color/mono, run separately in LabVIEW).
%
% ASSUMPTIONS (check these before running):
%   - USB-6451 is "Dev1" in NI MAX -- CHANGE if different
%   - Fan on Dev1/ao0, hot-wire on Dev1/ai0, trigger pulse on Dev1/port0/line0 (DIO0)
%   - Fan ramp: 0->7V over 30s, hold 5s, ramp down -- adjust as needed
%   - Hot-wire: 4 kHz sample rate
%   - Trigger: single TTL pulse, active high, ~10 ms wide (adjust WIDTH if 6621 needs longer)
%
% NOT verified -- confirm before relying on this:
%   - Whether DIO0 can be written while the AI task is running on the same device
%     (background AI acquisition + separate DIO write should be fine on the 6451,
%     but if you see conflicts, move the trigger to a counter output instead)

clear; clc;

DEV_ID = "Dev4";   % <-- confirm this matches NI MAX for the USB-6451

FAN_V_START = 1.5;
FAN_V_END   = 8;
FAN_RAMPUP_T  = 60;   % seconds, ramp up duration
FAN_RAMPDOWN_T  = 5;   % seconds, ramp up duration
FAN_HOLD_T  = 30;    % seconds, hold at V_END
FAN_DT      = 0.5;  % seconds, step interval

HOTWIRE_FS       = 4000;   % Hz
DELAY_BEFORE_TRIG = 30;    % seconds after t=0 before firing camera trigger

TRIGGER_PULSE_WIDTH = 0.01; % seconds (10 ms) -- confirm 6621 min pulse width requirement

%% Set up hot-wire AI task (background acquisition, starts immediately)
HOTWIRE_TOTAL_T   = FAN_RAMPUP_T + FAN_HOLD_T + FAN_RAMPDOWN_T + 5; % rough guess, adjust to your real run length

hw = setup_hotwire_daq(DEV_ID, HOTWIRE_FS, HOTWIRE_TOTAL_T);
% start(hw);   % begins acquiring in the background
disp('Hot-wire acquisition started.');

%% Set up fan AO channel
fanD = daq("ni");
addoutput(fanD, DEV_ID, "ao0", "Voltage");

%% Run fan ramp (blocking, foreground) while hot-wire logs in background
disp('Starting fan ramp...');
run_fan_ramp(fanD, FAN_V_START, FAN_V_END, FAN_RAMPUP_T, FAN_HOLD_T, FAN_RAMPDOWN_T, FAN_DT);
disp('Fan ramp complete (holding/ramped down per profile).');

%% Wait until 30 s mark (measured from t=0, i.e. fan/hot-wire start), then trigger cameras
elapsed = FAN_RAMP_T + FAN_HOLD_T + FAN_RAMP_T; % time already spent in ramp above
remaining_wait = DELAY_BEFORE_TRIG - elapsed;
if remaining_wait > 0
    fprintf('Waiting additional %.1f s before camera trigger...\n', remaining_wait);
    pause(remaining_wait);
end

send_trigger_pulse(DEV_ID, "port0/line0", TRIGGER_PULSE_WIDTH);
disp('Camera trigger pulse sent to PCI-6621 (CTR7 via DIO0).');

%% Let hot-wire finish logging
disp('Waiting for hot-wire acquisition to finish...');
wait(hw);
data = read(hw, "all");
disp('Hot-wire acquisition complete.');

% Optionally save
save(sprintf('hotwire_%s.mat', datestr(now,'yyyymmdd_HHMMSS')), 'data');

disp('Experiment complete.');
