function d = setup_hotwire_daq(devID, fs, totalTime)
% SETUP_HOTWIRE_DAQ  Configure background AI acquisition for 3 hot-wire probes.
%
%   d = setup_hotwire_daq(devID, fs, totalTime)
%
%   devID     - string, e.g. "Dev1"
%   fs        - sample rate in Hz (e.g. 4000)
%   totalTime - total duration to acquire, in seconds
%
% Three probes, each differential per the USB-6451's fixed AI pairing.
% Only the positive input is named in software -- the negative pin is
% hardwired by the device and must NOT be specified separately:
%   Probe 1 - ai0 (paired with ai8)
%   Probe 2 - ai1 (paired with ai9)
%   Probe 3 - ai2 (paired with ai10)
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

    probeChannels = ["ai0", "ai1", "ai2"];
    for i = 1:numel(probeChannels)
        ch = addinput(d, devID, probeChannels(i), "Voltage");
        ch.TerminalConfig = "Differential";
        ch.Range = [-10 10];  
    end

    d.Rate = fs;
end
