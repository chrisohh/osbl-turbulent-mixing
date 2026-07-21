function d = setup_camera_triggers(devID, camConfig)
% SETUP_CAMERA_TRIGGERS  Configure continuous counter-output pulse trains
% for camera triggering (PCI-6621).
%
%   d = setup_camera_triggers(devID, camConfig)
%
%   devID     - string, e.g. "Dev3"
%   camConfig - struct array with fields:
%                 name      - string, camera label (e.g. "IR"), used only for logging
%                 ctr       - string, counter channel (e.g. "ctr2")
%                 freq      - pulse frequency in Hz (e.g. 50)
%                 dutyCycle - fraction 0-1 (e.g. 0.5 for 50%)
%                 delay     - initial delay in seconds before first pulse (0 for none)
%
% ASSUMPTIONS (check these before running):
%   - devID is the PCI-6621's actual NI MAX device string
%   - each camConfig.ctr names a valid, unique counter on that device (ctr0-ctr3)
%
% NOT VERIFIED -- confirm before relying on this:
%   - Whether DutyCycle is directly settable on a PulseGeneration channel in
%     your DAQ Toolbox version, or whether PulseWidth must be set instead --
%     after the first addoutput call, run properties(ch) at the command line
%     to see what's actually exposed, and adjust the property-set block below
%     if needed.
%   - Whether InitialDelay is the correct property name for a per-channel
%     start delay in your toolbox version.

    d = daq("ni");

    for i = 1:numel(camConfig)
        cfg = camConfig(i);

        ch = addoutput(d, devID, cfg.ctr, "PulseGeneration");
        ch.Frequency = cfg.freq;
        ch.DutyCycle = cfg.dutyCycle;   % <-- NOT VERIFIED, see header

        if cfg.delay ~= 0
            ch.InitialDelay = cfg.delay;   % <-- NOT VERIFIED, see header
        end

        fprintf('%s configured on %s: %.1f Hz, %.1f%% duty\n', ...
            cfg.name, cfg.ctr, cfg.freq, cfg.dutyCycle * 100);
    end
end
