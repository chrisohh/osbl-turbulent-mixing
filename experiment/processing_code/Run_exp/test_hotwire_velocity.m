%% test_hotwire_velocity.m
% Acquire from the 3 hot-wire probes and convert to velocity (U,V,W) live,
% using the same calibration pipeline as the offline convertE2U.m. No fan,
% no triggers -- manually blow air on the probes during the window and
% confirm the velocity responds.
%
% ASSUMPTIONS (check before running):
%   - USB-6451 is "Dev4" in NI MAX -- CHANGE if different
%   - Probes wired per setup_hotwire_daq.m: ai0/ai8, ai1/ai9, ai2/ai10 (differential)
%   - CAL_FILE points to the StreamWare calibration/header .txt for THIS probe
%
% NOT verified: no fluid-temperature channel here, so temperature correction
% is skipped (Ta omitted from convert_E2U_fn). Absolute velocities below the
% probe's calibration floor (~0.5 m/s) rely on the origin-to-floor branch.

clear; clc;

DEV_ID   = "Dev4";     % <-- confirm this matches NI MAX for the USB-6451
FS       = 4000;       % Hz
DURATION = 15;         % seconds -- have air ready to blow during this window

CAL_FILE = 'D:\Chris\osbl-turbulent-mixing\experiment\data\260721\probe3.txt'; % <-- this probe's cal/header
CTA_DIR  = 'D:\Chris\osbl-turbulent-mixing\experiment\processing_code\CTA';     % parse_calibration + convert_E2U_fn

addpath(CTA_DIR);   % so convert_E2U_fn / parse_calibration are on the path

%% Acquire
hw = setup_hotwire_daq(DEV_ID, FS, DURATION);

fprintf('Acquiring for %d s -- blow air on the probes now...\n', DURATION);
data = read(hw, seconds(DURATION));
disp('Acquisition complete.');

E1 = data{:, 1};   % Probe 1 (ai0)
E2 = data{:, 2};   % Probe 2 (ai1)
E3 = data{:, 3};   % Probe 3 (ai2)
t  = seconds(data.Time);

%% Convert voltages -> velocities (probe coordinates)
[U, V, W] = convert_E2U_fn(E1, E2, E3, CAL_FILE);
U_mag = sqrt(U.^2 + V.^2 + W.^2);

fprintf('Mean |U| = %.3f m/s  (U=%.3f, V=%.3f, W=%.3f)\n', ...
    mean(U_mag), mean(U), mean(V), mean(W));

%% Plot raw voltages and converted velocity
figure;
subplot(2,1,1);
plot(t, [E1 E2 E3]);
ylabel('Voltage (V)');
legend('Probe 1 (ai0)', 'Probe 2 (ai1)', 'Probe 3 (ai2)');
title('Raw hot-wire voltages');
grid on;

subplot(2,1,2);
plot(t, U, t, V, t, W, 'k', 'LineWidth', 1);
xlabel('Time (s)');
ylabel('Velocity (m/s)');
legend('U', 'V', 'W');
title('Converted velocity (probe coordinates)');
grid on;
