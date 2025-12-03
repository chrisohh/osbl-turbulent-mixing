% TRIGGER GENERATOR - Using Dev1/port0/line0
% Physical: P0.0 → USER 1 → PFI0 (on BNC-2110)

clear; clc;

fprintf('=== Trigger Generator (Dev1) ===\n');
fprintf('Using: Dev1/port0/line7 (P0.7)\n');
fprintf('Path: P0.0 → USER 1 → PFI0\n\n');

% Create digital output
dq = daq("ni");
ch = addoutput(dq, "Dev1", "port0/line7", "Digital");

% Initialize to LOW
fprintf('Initializing to LOW...\n');
write(dq, 0);
pause(0.5);

% Wait
fprintf('\n=== READY TO TRIGGER ===\n');
fprintf('Make sure wave gauge receiver is waiting!\n\n');
input('Press ENTER to send trigger...', 's');

% Send pulse
fprintf('Sending trigger pulse...');
write(dq, 1);     % HIGH
pause(0.01);      % 10ms
write(dq, 0);     % LOW
fprintf(' DONE!\n');

fprintf('\nTrigger sent! Check receiver window.\n');

clear dq ch;
