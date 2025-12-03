% Simple Wave Gauge Acquisition with External Trigger
% Chris - SOARS Lab
% Start simple: just wave gauge + external trigger

clear; clc; close all;

%% ==================== SIMPLE CONFIGURATION ====================

% Wave Gauge Settings
SAMPLE_RATE = 1000;              % Hz (1 kHz is plenty for waves)
DURATION = 10;                   % seconds
WAVEGAUGE_CHANNEL = 'ai0';       % Wave gauge connected to AI0
WAVEGAUGE_CALIBRATION = 0.2;     % V/cm (from amplifier)

% Trigger Settings
USE_TRIGGER = true;              % Set to false to test without trigger
TRIGGER_CHANNEL = 'PFI0'; % Change based on what WWT team provides
                                 % Options: 'port0/line2', 'port0/line7', 'PFI0', etc.

%% ==================== SETUP ====================

fprintf('=== Simple Wave Gauge Acquisition ===\n\n');
fprintf('Sample Rate: %d Hz\n', SAMPLE_RATE);
fprintf('Duration: %.1f sec\n', DURATION);
fprintf('Trigger: %s\n', TRIGGER_CHANNEL);

% Create DAQ
dq = daq("ni");
dq.Rate = SAMPLE_RATE;

% Add wave gauge channel
ch = addinput(dq, "Dev1", WAVEGAUGE_CHANNEL, "Voltage");
ch.TerminalConfig = 'SingleEnded';
fprintf('Wave gauge: Dev1/%s\n', WAVEGAUGE_CHANNEL);

%% ==================== TRIGGER SETUP ====================

if USE_TRIGGER
    fprintf('\n=== EXTERNAL TRIGGER MODE ===\n');
    fprintf('Trigger input: Dev1/%s\n', TRIGGER_CHANNEL);
    
    % Remove any existing triggers
    try
        removetrigger(dq);
    catch
    end
    
    % Add digital trigger
    trig = addtrigger(dq, "Digital", "StartTrigger", ...
        "External", sprintf("Dev1/%s", TRIGGER_CHANNEL));
    trig.Condition = "RisingEdge";
    
    fprintf('\n>>> ARMED AND WAITING FOR TRIGGER <<<\n');
    fprintf('Trigger the WWT system to start acquisition\n\n');
else
    fprintf('\n=== NO TRIGGER - Starting immediately ===\n');
    input('Press ENTER to start...', 's');
end

%% ==================== ACQUIRE DATA ====================

if USE_TRIGGER
    fprintf('Waiting for trigger on %s...\n', TRIGGER_CHANNEL);
    fprintf('Send trigger from other MATLAB window (you have ~10 seconds)!\n\n');
end

tic;

% Read data - will wait for trigger if enabled
try
    fprintf('Aquiring data (%.2f sec)\n', DURATION);
    data = read(dq, seconds(DURATION));
    elapsed = toc;
    fprintf('Acquisition complete! (%.2f sec)\n', elapsed);
catch ME
    if contains(ME.message, 'Timeout') || contains(ME.message, 'timed out')
        fprintf('\n=== TRIGGER TIMEOUT ===\n');
        fprintf('No trigger received after waiting.\n\n');
        fprintf('Troubleshooting checklist:\n');
        fprintf('  1. Is BNC cable connected? Dev2 → Dev1/PFI0\n');
        fprintf('  2. Is trigger generator running in another MATLAB window?\n');
        fprintf('  3. Did you press ENTER in the trigger generator?\n');
        fprintf('  4. Try running trigger generator first, then receiver\n\n');
        error('Trigger timeout - see checklist above');
    else
        rethrow(ME);
    end
end


%% ==================== PROCESS DATA ====================

% Extract time and voltage
time = (0:height(data)-1)' / SAMPLE_RATE;
voltage = data.Variables;

% Convert to wave height (cm)
still_water_voltage = mean(voltage(1:min(100, length(voltage))));
wave_height_cm = (voltage - still_water_voltage) / WAVEGAUGE_CALIBRATION;


%% ==================== PLOT ====================

figure('Position', [100 100 1200 600]);

% Plot 1: Voltage
subplot(2,1,1);
plot(time, voltage, 'b-', 'LineWidth', 1);
hold on;
yline(still_water_voltage, 'r--', 'Still Water', 'LineWidth', 1.5);
grid on;
xlabel('Time (s)');
ylabel('Voltage (V)');
title('Wave Gauge Voltage');

% Plot 2: Wave Height
subplot(2,1,2);
plot(time, wave_height_cm, 'b-', 'LineWidth', 1);
hold on;
yline(0, 'r--', 'Still Water Level', 'LineWidth', 1.5);
grid on;
xlabel('Time (s)');
ylabel('Wave Height (cm)');
title('Wave Height (relative to still water)');

%% ==================== SAVE DATA ====================

% Save with timestamp
timestamp = datestr(now, 'yyyymmdd_HHMMSS');
filename = sprintf('wave_gauge_%s.mat', timestamp);

save_data = struct();
save_data.time = time;
save_data.voltage = voltage;
save_data.wave_height_cm = wave_height_cm;
save_data.sample_rate = SAMPLE_RATE;
save_data.calibration = WAVEGAUGE_CALIBRATION;
save_data.still_water_voltage = still_water_voltage;
save_data.trigger_used = USE_TRIGGER;
save_data.trigger_channel = TRIGGER_CHANNEL;
save_data.timestamp = timestamp;

save(filename, 'save_data');
fprintf('\nData saved: %s\n', filename);

% Export to CSV
csv_filename = sprintf('wave_gauge_%s.csv', timestamp);
csv_table = table(time, voltage, wave_height_cm, ...
    'VariableNames', {'Time_s', 'Voltage_V', 'WaveHeight_cm'});
writetable(csv_table, csv_filename);
fprintf('CSV exported: %s\n', csv_filename);

fprintf('\n=== Acquisition Complete ===\n');