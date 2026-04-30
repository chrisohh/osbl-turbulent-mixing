% Script to extract and plot wind velocity (Sm) from Vaisala data
% File: viasala_wind_flowStraightner_glasschannel_10cmML_20260116.txt
% Data sampled at 1 Hz

% File path
filename = '\\Airseaserver28\D\HLAB_2026\Vaisala\viasala_wind_flowStraightner_glasschannel_10cmML_20260116.txt';

% Read the file
fid = fopen(filename, 'r');
if fid == -1
    error('Cannot open file: %s', filename);
end

% Initialize arrays
Sm = [];
line_count = 0;

% Read line by line
while ~feof(fid)
    line = fgetl(fid);
    if ischar(line)
        line_count = line_count + 1;

        % Extract Sm value using regular expression
        % Pattern: Sm=X.XM where X.X is the value
        pattern = 'Sm=(\d+\.\d+)M';
        tokens = regexp(line, pattern, 'tokens');

        if ~isempty(tokens)
            Sm(end+1) = str2double(tokens{1}{1});
        end
    end
end
fclose(fid);

% Create time vector (1 Hz sampling)
time = (0:length(Sm)-1)';  % Time in seconds

% Plot the wind velocity time series
figure('Position', [100, 100, 1000, 500]);
plot(time, Sm, 'b-', 'LineWidth', 1.5);
xlabel('Time (s)');
ylabel('Wind Velocity, S_m (m/s)');
title('Vaisala Wind Velocity Time Series with Flow Straightener 10cm Fresh Water');
grid on;
set(gca, 'FontSize', 11);

time_honeycomb=time;
velocity_honeycomb=Sm;
% File path
filename = '\\Airseaserver28\D\HLAB_2026\Vaisala\viasala_wind_open_lid_glasschannel_10cmML_20260116.txt';

% Read the file
fid = fopen(filename, 'r');
if fid == -1
    error('Cannot open file: %s', filename);
end

% Initialize arrays
Sm = [];
line_count = 0;

% Read line by line
while ~feof(fid)
    line = fgetl(fid);
    if ischar(line)
        line_count = line_count + 1;

        % Extract Sm value using regular expression
        % Pattern: Sm=X.XM where X.X is the value
        pattern = 'Sm=(\d+\.\d+)M';
        tokens = regexp(line, pattern, 'tokens');

        if ~isempty(tokens)
            Sm(end+1) = str2double(tokens{1}{1});
        end
    end
end
fclose(fid);

% Create time vector (1 Hz sampling)
time = (0:length(Sm)-1)';  % Time in seconds
hold on
plot(time, Sm, 'k-', 'LineWidth', 1.5);
xlabel('Time (s)');
ylabel('Wind Velocity (m/s)');
title('Vaisala Ultrasonic Wind Velocity 10cm Fresh + 49cm Salt Water');
set(gca, 'FontSize', 12);
legend('with flow straightner','without flow straightner')

grid off

time_openlid=time;
velocity_openlid=Sm;

%% Save the extracted data
save('vaisala_wind_10cmFresh_49cmSalt.mat', 'time_openlid', 'time_honeycomb', 'velocity_openlid','velocity_honeycomb');