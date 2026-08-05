function d = setup_hotwire_daq(devID, fs, chanConfig)
% SETUP_HOTWIRE_DAQ  Configure background AI acquisition for hot-wire probes
% and/or a reference probe.
%
%   d = setup_hotwire_daq(devID, fs)              % default: 3 hot-wire probes
%   d = setup_hotwire_daq(devID, fs, chanConfig)  % explicit channel list
%
%   devID      - string, e.g. "Dev4"
%   fs         - sample rate in Hz (e.g. 4000)
%   chanConfig - struct array selecting which AI channels to log, fields:
%                  name  - label, used only for logging
%                  chan  - channel id, e.g. "ai0"
%                  term  - "Differential" or "SingleEnded"
%                  range - [min max] volts
%                If omitted, empty, or numeric (the old unused totalTime
%                argument), defaults to the three differential hot-wire probes.
%
% Channels are added in the order given, which is also the column order of
% the timetable returned by read() -- callers should index by that order
% rather than assuming fixed columns.
%
% Differential channels use the USB-6451's fixed AI pairing. Only the positive
% input is named in software; the negative pin is hardwired by the device and
% must NOT be listed separately:
%   ai0<->ai8, ai1<->ai9, ai2<->ai10, ai3<->ai11, ai4<->ai12
%
% NOT VERIFIED: the exact background-acquisition call (start/read pattern)
% differs across MATLAB Data Acquisition Toolbox versions. Two common
% patterns are shown below in run_experiment.m comments -- test with a
% short duration (e.g. 2 seconds) first to confirm which one your MATLAB
% version accepts before running the full experiment.
%
% Also NOT VERIFIED: whether this device accepts mixed terminal configurations
% (Differential hot-wire channels alongside a SingleEnded reference channel)
% in one task. If addinput errors on the mixed task, put all channels on the
% same terminal config.
%
% Also NOT VERIFIED: hot-wire CTA voltage output range -- confirm against
% your Dantec system spec before trusting the Range settings passed in.

    if nargin < 3 || isempty(chanConfig) || isnumeric(chanConfig)
        chanConfig = struct( ...
            'name',  {'Probe1',       'Probe2',       'Probe3'}, ...
            'chan',  {"ai0",          "ai1",          "ai2"}, ...
            'term',  {"Differential", "Differential", "Differential"}, ...
            'range', {[-10 10],       [-10 10],       [-10 10]});
    end

    d = daq("ni");

    for i = 1:numel(chanConfig)
        cfg = chanConfig(i);
        ch = addinput(d, devID, cfg.chan, "Voltage");
        ch.TerminalConfig = cfg.term;
        ch.Range = cfg.range;
        fprintf('  AI col %d: %s on %s/%s (%s, %g..%g V)\n', ...
            i, cfg.name, devID, cfg.chan, cfg.term, cfg.range(1), cfg.range(2));
    end

    d.Rate = fs;
end
