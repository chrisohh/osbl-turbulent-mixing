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
%   - Whether InitialDelay is the correct property name for a per-channel
%     start delay in your toolbox version.
%
% Pulse width vs. duty cycle: if the channel exposes a PulseWidth property
% (seconds), that's used directly -- more precise for cameras like Color
% that need an exact pulse width (e.g. 700us) regardless of frequency.
% Falls back to DutyCycle (fraction of period) if PulseWidth isn't available
% in this DAQ Toolbox version. Prints which one was actually used per channel.

    d = daq("ni");

    for i = 1:numel(camConfig)
        cfg = camConfig(i);

        ch = addoutput(d, devID, cfg.ctr, "PulseGeneration");
        ch.Frequency = cfg.freq;

        ch.DutyCycle = cfg.dutyCycle;

        fprintf('%s configured on %s: %.1f Hz, exposure %.0f us (duty %.4f), idle %s\n', ...
            cfg.name, cfg.ctr, cfg.freq, (cfg.dutyCycle / cfg.freq) * 1e6, ...
            cfg.dutyCycle, string(ch.IdleState));

        if cfg.delay ~= 0
            ch.InitialDelay = cfg.delay;   % <-- NOT VERIFIED, see header
        end
    end
end
