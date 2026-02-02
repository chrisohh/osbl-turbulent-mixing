% wind_speed_profile.m
% Generates a wind speed profile table with:
% - Constant acceleration ramp from 1 m/s to 10 m/s over 60 seconds
% - Dwell at 10 m/s for 2 minutes (120 seconds)

% Define parameters
v_start = 1;        % Starting wind speed (m/s)
v_end = 10;         % Ending wind speed (m/s)
t_ramp = 60;        % Total ramp duration (s)
t_dwell_final = 120;      % Final dwell duration at max speed (s)

% Dwell time at each step during ramp (adjust based on system response)
% Recommended: 1-2 seconds for typical wind-wave flume PID controllers
dt_ramp = 10;        % Dwell time per step during ramp (s)

% Calculate wind speed step size to achieve constant acceleration
% For constant acceleration: delta_v = acceleration * dt_ramp
% Total ramp time: t_ramp = n_steps * dt_ramp
% Total velocity change: v_end - v_start = n_steps * dv
% Therefore: dv = (v_end - v_start) * dt_ramp / t_ramp
dv = (v_end - v_start) * dt_ramp / t_ramp;  % Wind speed increment (m/s)
n_steps = round(t_ramp / dt_ramp);          % Number of steps

% Generate ramp-up phase with constant acceleration
wind_speed_ramp = (v_start:dv:(v_start + (n_steps)*dv))';
dwell_ramp = dt_ramp * ones(size(wind_speed_ramp));

% Adjust last dwell time of ramp if needed
% (since we go from 0 to t_ramp inclusive)

% Add dwell phase at constant 10 m/s
wind_speed_dwell = v_end;
dwell_constant = t_dwell_final;

% Combine ramp and dwell phases
wind_speed = [wind_speed_ramp; wind_speed_dwell];
dwell = [dwell_ramp; dwell_constant];

% Create table
wind_profile = table(wind_speed, dwell, ...
    'VariableNames', {'Airspeed_m_s_', 'Dwell_sec_'});

% Display summary
fprintf('\n=== Wind Speed Profile Summary ===\n');
fprintf('Ramp: %.1f to %.1f m/s over %.0f seconds\n', v_start, v_end, t_ramp);
fprintf('Dwell per step during ramp: %.2f s\n', dt_ramp);
fprintf('Wind speed increment per step: %.4f m/s\n', dv);
fprintf('Number of ramp steps: %d\n', n_steps);
fprintf('Acceleration: %.4f m/s²\n', dv/dt_ramp);
fprintf('Final dwell at %.1f m/s: %.0f s\n', v_end, t_dwell_final);
fprintf('Total duration: %.0f s\n', sum(dwell));
fprintf('==================================\n\n');

% Display table
disp('Wind Speed Profile:');
disp(wind_profile);

%% Optional: Plot the profile
figure;
time_cumulative = [0; cumsum(dwell)];
wind_speed_plot = [wind_speed(1); wind_speed];

% subplot(2,1,1);
stairs(time_cumulative, wind_speed_plot, 'LineWidth', 2);
xlabel('Time (s)');
ylabel('Wind Speed (m/s)');
title('Wind Speed Profile');
grid on;

% subplot(2,1,2);
% bar(1:length(dwell), dwell);
% xlabel('Step Number');
% ylabel('Dwell Time (s)');
% title('Dwell Time at Each Step');
% grid on;

%% Save to CSV file with custom headers
csv_filename = 'wind_speed_profile.csv';
[fid, msg] = fopen(csv_filename, 'w');
if fid == -1
    warning('Could not create CSV file: %s', msg);
    fprintf('Make sure the file is not open in another program.\n');
else
    fprintf(fid, 'Airspeed (m/s),Dwell (sec)\n');
    for i = 1:length(wind_speed)
        fprintf(fid, '%.4f,%.2f\n', wind_speed(i), dwell(i));
    end
    fclose(fid);
    fprintf('Profile saved to: %s\n', csv_filename);
end
