function send_trigger_pulse(devID, lineID, pulseWidth)
% SEND_TRIGGER_PULSE  Write a single TTL high pulse on a digital line.
%
%   send_trigger_pulse(devID, lineID, pulseWidth)
%
%   devID      - string, e.g. "Dev1"
%   lineID     - string, e.g. "port0/line0" (DIO0)
%   pulseWidth - pulse duration in seconds (e.g. 0.01 for 10 ms)
%
% NOT VERIFIED: minimum pulse width the PCI-6621 needs to reliably
% register a start trigger on its PFI input -- 10 ms is a safe guess
% for a software-timed pulse, but check your 6621 trigger spec or just
% test on the scope to confirm it's being seen cleanly.
%
% NOTE: this uses software timing (write True, pause, write False), so
% the pulse width will have some jitter (typically a few ms) due to
% MATLAB/OS scheduling. If the 6621 needs a very precise or very short
% pulse, this software approach may not be reliable -- a hardware-timed
% digital pulse (via a counter output) would be more robust.

    d = daq(devID);
    ch = addoutput(d, devID, lineID, "Digital");

    write(d, 1);
    pause(pulseWidth);
    write(d, 0);
end
