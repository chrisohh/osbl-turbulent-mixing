function nudge_camera_trigger(devID, ctr, freq, dutyCycle, duration)
% NUDGE_CAMERA_TRIGGER  Emit a short burst of trigger pulses on one counter.
%
%   nudge_camera_trigger(devID, ctr)
%   nudge_camera_trigger(devID, ctr, freq, dutyCycle, duration)
%
%   devID     - string, e.g. "Dev3" (the PCI-6621)
%   ctr       - counter channel, e.g. "ctr2" (IR), "ctr1" (Color), "ctr3" (Mono)
%   freq      - pulse frequency in Hz (default 50)
%   dutyCycle - fraction 0-1 (default 0.5; irrelevant for just registering an edge)
%   duration  - burst length in seconds (default 0.5, i.e. ~25 pulses at 50 Hz)
%
% WHY THIS EXISTS: the camera software treats Record / Stop Record as PENDING
% states that are only processed when the next trigger arrives. If you arm or
% stop a camera while the counters are idle, it sits waiting forever. This
% sends a short burst so the pending command actually takes effect, without
% running a full experiment.
%
% NOTE: if a previous task still holds this counter (e.g. run_experiment.m was
% interrupted with Ctrl-C, so stop(camD) never ran), this will fail with a
% resource-busy error. Run cleanup_daq first in that case.

    if nargin < 3 || isempty(freq),      freq = 50;        end
    if nargin < 4 || isempty(dutyCycle), dutyCycle = 0.5;  end
    if nargin < 5 || isempty(duration),  duration = 0.5;   end

    d = daq("ni");
    ch = addoutput(d, devID, ctr, "PulseGeneration");
    ch.Frequency = freq;
    ch.DutyCycle = dutyCycle;

    try
        start(d, "continuous");
        fprintf('Pulsing %s/%s at %.0f Hz for %.2f s (~%.0f pulses)...\n', ...
            devID, ctr, freq, duration, freq * duration);
        pause(duration);
        stop(d);
        disp('Burst complete -- the camera should have processed its pending command.');
    catch ME
        stop(d);
        rethrow(ME);
    end
end
