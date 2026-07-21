function d = setup_hotwire_daq(devID, fs, totalTime)
% SETUP_HOTWIRE_DAQ  Configure background AI acquisition for hot-wire.
%
%   d = setup_hotwire_daq(devID, fs, totalTime)
%
%   devID     - string, e.g. "Dev1"
%   fs        - sample rate in Hz (e.g. 4000)
%   totalTime - total duration to acquire, in seconds
%
% NOT VERIFIED: the exact background-acquisition call (start/read pattern)
% differs across MATLAB Data Acquisition Toolbox versions. Two common
% patterns are shown below in run_experiment.m comments -- test with a
% short duration (e.g. 2 seconds) first to confirm which one your MATLAB
% version accepts before running the full experiment.
%
% Also NOT VERIFIED: hot-wire CTA voltage output range -- confirm against
% your Dantec system spec before trusting the Range setting below.

    d = daq("ni");
    ch = addinput(d, devID, "ai0", "Voltage");
    ch.Range = [-10 10];   % <-- CONFIRM against your CTA's actual output range

    d.Rate = fs;
end
