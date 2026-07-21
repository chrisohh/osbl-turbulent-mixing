%% test_camera_triggers.m
% Standalone scope-validation test for two of three planned camera-trigger
% counters (IR/ctr2, Mono/ctr3) on the PCI-6621, run in isolation before
% wiring counter-based triggering into a combined experiment script.
%
% ASSUMPTIONS (check these before running):
%   - PCI-6621 is "Dev3" in NI MAX -- CHANGE if different (not yet confirmed
%     anywhere else in this project)
%   - ctr2 = IR, ctr3 = Mono, per prior LabVIEW wiring
%   - scope probes are connected to the correct PFI/counter output lines on Dev3
%
% NOT verified -- confirm before relying on this:
%   - Whether start(d, "continuous") or plain start(d) is the correct call
%     for continuous PulseGeneration counter output in your DAQ Toolbox
%     version -- if "continuous" errors, try start(camD) alone.
%   - The daqlist("ni") table's exact column name (DeviceID assumed below) --
%     check the disp(devices) output on first run.
%   - DutyCycle/InitialDelay property behavior -- see setup_camera_triggers.m

clear; clc;

DEV_CAMERAS   = "Dev3";   % <-- confirm this matches NI MAX for the PCI-6621
TEST_DURATION = 5;        % seconds, scope-capture run length

%% Check that the camera DAQ device is present
disp('Checking for connected NI devices...');
devices = daqlist("ni");
disp(devices);

if ~any(strcmp(devices.DeviceID, DEV_CAMERAS))
    error('Device %s not found in daqlist -- check NI MAX / cabling before proceeding.', DEV_CAMERAS);
end
disp(['Found ' char(DEV_CAMERAS) ' -- proceeding.']);

%% Camera config -- IR and Mono only
% Color/ctr1 intentionally excluded from this test: its ~700us duty-cycle
% math has not been confirmed yet, and IR+Mono are being validated on the
% scope first before adding a third channel.
camConfig = struct( ...
    'name',      {'IR',   'Mono'}, ...
    'ctr',       {'ctr2', 'ctr3'}, ...
    'freq',      {50,     50}, ...
    'dutyCycle', {0.5,    0.5}, ...
    'delay',     {0,      0});

%% Set up camera counters (not started yet)
disp('Setting up camera trigger counters...');
camD = setup_camera_triggers(DEV_CAMERAS, camConfig);
disp('Camera trigger channels configured (IR/ctr2, Mono/ctr3).');

%% Start, hold for scope capture, stop
try
    disp('Starting counter outputs...');
    start(camD, "continuous");   % <-- NOT VERIFIED, see header

    fprintf('Running for %.1f s -- check scope now on %s ctr2 (IR) and ctr3 (Mono)...\n', ...
        TEST_DURATION, DEV_CAMERAS);
    pause(TEST_DURATION);

    stop(camD);
    disp('Counter outputs stopped cleanly.');
catch ME
    disp('Interrupted -- stopping counter outputs.');
    stop(camD);
    rethrow(ME);
end

%% Wrap-up
disp('Test complete: IR (ctr2) and Mono (ctr3) trigger signals exercised for scoping.');
disp('Color/ctr1 was NOT included in this test -- pulse-width math still unconfirmed.');
disp('Next: verify both channels on scope (~50 Hz, ~50% duty), then add Color/ctr1,');
disp('and eventually integrate counter-based triggering into a full combined experiment script.');
