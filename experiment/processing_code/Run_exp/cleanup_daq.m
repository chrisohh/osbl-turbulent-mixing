%% cleanup_daq.m
% Run this after any interrupted experiment (Ctrl-C, error, closed figure
% mid-run). MATLAB's Ctrl-C does NOT execute catch blocks, so run_experiment.m's
% cleanup is skipped on a manual interrupt -- leaving the fan at its last
% commanded voltage, the camera counters still pulsing, and the hot-wire task
% still acquiring.
%
% This script zeroes the fan FIRST (the only physically-consequential state),
% then releases every DAQ task on all devices.
%
% ASSUMPTIONS (check these before running):
%   - USB-6451 (fan) is "Dev4" in NI MAX -- must match run_experiment.m
%   - Fan is on Dev4/ao0
%
% NOTE: daqreset invalidates any existing daq objects (camD, hw, fanD), so
% re-run run_experiment.m from the top afterwards rather than resuming.

DEV_ID = "Dev4";   % <-- confirm this matches NI MAX for the USB-6451

%% Zero the fan first -- do this before daqreset, which would leave the AO
%% at its last value rather than driving it low.
try
    fanCleanup = daq("ni");
    addoutput(fanCleanup, DEV_ID, "ao0", "Voltage");
    write(fanCleanup, 0);
    disp('Fan output set to 0 V.');
catch ME
    fprintf(2, 'Could not zero the fan (%s): %s\n', ME.identifier, ME.message);
    fprintf(2, 'CHECK THE FAN PHYSICALLY -- it may still be running.\n');
end

%% Stop any tasks left running in the base workspace (counters, hot-wire).
for name = ["camD", "hw", "fanD", "fanCleanup"]
    if evalin('base', sprintf('exist(''%s'', ''var'')', name))
        try
            evalin('base', sprintf('stop(%s);', name));
            fprintf('Stopped %s.\n', name);
        catch
            % Already stopped, or not a running task -- nothing to do.
        end
    end
end

%% Release all DAQ hardware
daqreset;
disp('daqreset complete -- all DAQ tasks released.');
disp('Camera counters are no longer pulsing; re-run run_experiment.m from the top.');
