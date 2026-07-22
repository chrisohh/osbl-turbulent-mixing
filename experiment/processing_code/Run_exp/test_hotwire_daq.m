%% test_hotwire_daq.m
% Simple standalone test: acquire from the 3 hot-wire probes for a fixed
% duration and plot the raw voltage traces. No fan, no camera trigger, no
% StreamWare pulse -- manually blow air on the probes during the
% acquisition window and confirm each channel's voltage responds.
%
% ASSUMPTIONS (check before running):
%   - USB-6451 is "Dev4" in NI MAX -- CHANGE if different
%   - Probes wired per setup_hotwire_daq.m: ai0/ai8, ai1/ai9, ai2/ai10 (differential)
%
% NOT verified -- confirm before relying on this:
%   - Whether ch.TerminalConfig = "Differential" is accepted by your DAQ
%     Toolbox version -- this will error immediately inside
%     setup_hotwire_daq(), before any acquisition starts, if not.
%   - This uses a blocking foreground read(hw, seconds(DURATION)) instead
%     of the background start/read pattern noted in setup_hotwire_daq.m,
%     specifically to sidestep that pattern's toolbox-version uncertainty
%     for this quick test.

clear; clc;

DEV_ID   = "Dev4";   % <-- confirm this matches NI MAX for the USB-6451
FS       = 4000;     % Hz
DURATION = 15;        % seconds -- have air ready to blow on the probes during this window

hw = setup_hotwire_daq(DEV_ID, FS, DURATION);

fprintf('Acquiring for %d s -- blow air on the probes now...\n', DURATION);
data = read(hw, seconds(DURATION));
disp('Acquisition complete.');

%% Plot raw voltage traces
t = seconds(data.Time);

figure;
plot(t, data.Variables);
xlabel('Time (s)');
ylabel('Voltage (V)');
legend('Probe 1 (ai0)', 'Probe 2 (ai1)', 'Probe 3 (ai2)');
title('Hot-wire probe response');
grid on;
