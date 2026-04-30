% Tri-axial Hot-Wire Probe Data Processing for Vertical Traverses
% Modular version with functions for single or multiple heights

%% Main Script
% Load calibration data
cal_file = 'G:\My Drive\OSBL\CalibratorTest1\probe2\probe2_50-200rpm_global_export.txt';
cal = parse_calibration(cal_file);

% Display extracted calibration
display_calibration(cal);

%% Setup reference velocity calibration
cal_ref = setup_reference_calibration();

%% Example 1: Process single height measurement
base_path = 'G:\My Drive\OSBL\CalibratorTest1\probe2\probe2_voltage\';
single_file = 'traverse_single_height.txt';

% For single height, just provide the height value (or leave empty)
profile_data = process_traverse_data(fullfile(base_path, single_file), ...
    cal, cal_ref, 0); % Single height at z=0mm

fprintf('\n=== Single Height Results ===\n');
fprintf('Mean: U=%.3f, V=%.3f, W=%.3f m/s, T=%.2f°C\n', ...
    profile_data.U_mean, profile_data.V_mean, profile_data.W_mean, profile_data.T_mean);
fprintf('RMS: U_rms=%.3f, V_rms=%.3f, W_rms=%.3f m/s\n', ...
    profile_data.U_rms, profile_data.V_rms, profile_data.W_rms);

%% Example 2: Process multiple heights (vertical profile)
traverse_file = 'traverse_400rpm_400mm.txt';
heights = 400;% [300,200,150,100,75,50]; % Specify your heights in mm
time_per_height = 60; % seconds at each height
sampling_rate = 4000; % Hz

% For multiple heights, provide all parameters
profile_data = process_traverse_data(fullfile(base_path, traverse_file), ...
    cal, cal_ref, heights, time_per_height, sampling_rate);

% Get data_collection for full time series at each height
[profile_data, data_collection] = process_traverse_data(fullfile(base_path, traverse_file), ...
    cal, cal_ref, heights, time_per_height, sampling_rate);

% Plot vertical profiles (only if multiple heights)
if length(heights) > 1
    plot_vertical_profiles(profile_data, heights);
    save_profile_results(profile_data, heights, 'vertical_profile_results.csv');
    
    % Optional: Plot time series for each height
    plot_time_series_by_height(data_collection);
end

%% ========================================================================
%  FUNCTION DEFINITIONS
%% ========================================================================

%% Function: Display calibration info
function display_calibration(cal)
    fprintf('Calibration Coefficients Extracted:\n');
    fprintf('Sensor 1: C0=%.4f, C1=%.4f, C2=%.4f, C3=%.4f, C4=%.4f\n', ...
        cal.C0_1, cal.C1_1, cal.C2_1, cal.C3_1, cal.C4_1);
    fprintf('Sensor 2: C0=%.4f, C1=%.4f, C2=%.4f, C3=%.4f, C4=%.4f\n', ...
        cal.C0_2, cal.C1_2, cal.C2_2, cal.C3_2, cal.C4_2);
    fprintf('Sensor 3: C0=%.4f, C1=%.4f, C2=%.4f, C3=%.4f, C4=%.4f\n', ...
        cal.C0_3, cal.C1_3, cal.C2_3, cal.C3_3, cal.C4_3);
    fprintf('\nYaw factors squared: k1²=%.6f, k2²=%.6f, k3²=%.6f\n', ...
        cal.k1_sq, cal.k2_sq, cal.k3_sq);
    fprintf('Pitch factors squared: h1²=%.6f, h2²=%.6f, h3²=%.6f\n', ...
        cal.h1_sq, cal.h2_sq, cal.h3_sq);
    fprintf('Reference Temperature: %.2f °C\n', cal.T_ref);
end

%% Function: Setup reference velocity calibration
function cal_ref = setup_reference_calibration()
    cal_ref.TR = 24.68;
    cal_ref.PR = 100.2;
    cal_ref.G1 = 60;
    cal_ref.G2 = 20;
    cal_ref.U0 = 0.05;
    cal_ref.U1 = 0.5;
    cal_ref.U2 = 8;
    cal_ref.U3 = 33;
end

%% Function: Process traverse data (single or multiple heights)
function [profile_data, data_collection] = process_traverse_data(filepath, cal, cal_ref, ...
    heights, time_per_height, sampling_rate, apply_T_correction)
    
    % Handle optional arguments
    if nargin < 4 || isempty(heights)
        heights = 0; % Default single height
    end
    if nargin < 5
        time_per_height = [];
    end
    if nargin < 6
        sampling_rate = [];
    end
    if nargin < 7
        apply_T_correction = false;
    end
    
    % Read data
    data = readmatrix(filepath);
    total_samples = size(data, 1);
    num_heights = length(heights);
    
    % Determine if single or multiple heights
    if num_heights == 1
        % Single height - process all data (remove NaN rows)
        fprintf('\nProcessing single height measurement:\n');
        valid_rows = ~isnan(data(:,1));
        data = data(valid_rows, :);
        fprintf('Total valid samples: %d\n', size(data,1));
        samples_per_height = size(data,1);
    else
        % Multiple heights - find segments separated by NaN rows
        if isempty(time_per_height) || isempty(sampling_rate)
            error('For multiple heights, must specify time_per_height and sampling_rate');
        end
        
        fprintf('\nProcessing vertical profile:\n');
        fprintf('Total rows in file: %d\n', total_samples);
        
        % Find NaN separator rows
        nan_rows = find(isnan(data(:,1)));
        fprintf('Found %d NaN separator rows\n', length(nan_rows));
        
        % Define segments based on NaN separators
        segment_starts = [1; nan_rows + 1];
        segment_ends = [nan_rows - 1; total_samples];
        
        % Remove any invalid segments
        valid_segments = segment_starts <= total_samples & segment_ends >= segment_starts;
        segment_starts = segment_starts(valid_segments);
        segment_ends = segment_ends(valid_segments);
        
        num_segments = length(segment_starts);
        fprintf('Number of data segments: %d\n', num_segments);
        
        if num_segments ~= num_heights
            warning('Number of segments (%d) does not match number of heights (%d)', ...
                num_segments, num_heights);
            num_heights = min(num_segments, num_heights);
        end
    end
    
    % Initialize profile arrays
    profile_data.U_mean = zeros(num_heights, 1);
    profile_data.V_mean = zeros(num_heights, 1);
    profile_data.W_mean = zeros(num_heights, 1);
    profile_data.T_mean = zeros(num_heights, 1);
    profile_data.U_ref_mean = zeros(num_heights, 1);
    profile_data.U_mag_mean = zeros(num_heights, 1);
    profile_data.U_rms = zeros(num_heights, 1);
    profile_data.V_rms = zeros(num_heights, 1);
    profile_data.W_rms = zeros(num_heights, 1);
    
    % Initialize data_collection structure (like original code)
    data_collection = struct();
    
    % Process each height
    for h = 1:num_heights
        if num_heights == 1
            % Single height - already cleaned data
            start_idx = 1;
            end_idx = size(data, 1);
        else
            % Multiple heights - use segment boundaries
            start_idx = segment_starts(h);
            end_idx = segment_ends(h);
        end
        
        fprintf('Height %d (%.1f mm): Processing rows %d to %d (%d samples)\n', ...
            h, heights(h), start_idx, end_idx, end_idx-start_idx+1);
        
        % Extract data for this height
        time = data(start_idx:end_idx, 1);
        E1 = data(start_idx:end_idx, 2);
        E2 = data(start_idx:end_idx, 3);
        E3 = data(start_idx:end_idx, 4);
        E_T = data(start_idx:end_idx, 5);
        E_ref = data(start_idx:end_idx, 6);
        
        % Remove any remaining NaN values within the segment
        valid_idx = ~isnan(E1) & ~isnan(E2) & ~isnan(E3) & ~isnan(E_T) & ~isnan(E_ref);
        time = time(valid_idx);
        E1 = E1(valid_idx);
        E2 = E2(valid_idx);
        E3 = E3(valid_idx);
        E_T = E_T(valid_idx);
        E_ref = E_ref(valid_idx);
        
        if isempty(E1)
            warning('No valid data for height %d (%.1f mm)', h, heights(h));
            continue;
        end
        
        % Process velocities
        [U, V, W, T, U_ref] = process_velocities(E1, E2, E3, E_T, E_ref, ...
            cal, cal_ref, apply_T_correction);
        
        % Calculate statistics
        U_mag = sqrt(U.^2 + V.^2 + W.^2);
        
        profile_data.U_mean(h) = mean(U);
        profile_data.V_mean(h) = mean(V);
        profile_data.W_mean(h) = mean(W);
        profile_data.T_mean(h) = mean(T);
        profile_data.U_ref_mean(h) = mean(U_ref);
        profile_data.U_mag_mean(h) = mean(U_mag);
        profile_data.U_rms(h) = std(U);
        profile_data.V_rms(h) = std(V);
        profile_data.W_rms(h) = std(W);
        
        % Store full time series in data_collection structure
        data_collection(h).height = heights(h);
        data_collection(h).time = time;
        data_collection(h).U = U;
        data_collection(h).V = V;
        data_collection(h).W = W;
        data_collection(h).T = T;
        data_collection(h).U_ref = U_ref;
        data_collection(h).U_mag = U_mag;
        data_collection(h).E_ref = E_ref;
        
        fprintf('  Results: U=%.3f, V=%.3f, W=%.3f m/s, T=%.2f°C\n', ...
            profile_data.U_mean(h), profile_data.V_mean(h), ...
            profile_data.W_mean(h), profile_data.T_mean(h));
    end
end

%% Function: Process velocities from raw voltages
function [U, V, W, T, U_ref] = process_velocities(E1, E2, E3, E_T, E_ref, ...
    cal, cal_ref, apply_T_correction)
    
    if nargin < 7
        apply_T_correction = false;
    end
    
    % Temperature calibration constants
    C0_T = 5.3065;
    C1_T = -25.309;
    C2_T = -0.7779;
    C3_T = -0.5073;
    
    % Temperature conversion
    T = C0_T + C1_T*log(E_T) + C2_T*log(E_T).^2 + C3_T*log(E_T).^3;
    
    % Reference velocity conversion
    E0 = 0;
    E1_ref = E0 + ((cal_ref.U1 - cal_ref.U0) / cal_ref.G1)^2;
    E2_ref = E1_ref + ((cal_ref.U2 - cal_ref.U1) / cal_ref.G2)^2;
    
    U_ref = zeros(size(E_ref));
    mask1 = (E_ref >= E0) & (E_ref < E1_ref);
    U_ref(mask1) = cal_ref.U0 + cal_ref.G1 * sqrt(E_ref(mask1) - E0);
    mask2 = (E_ref >= E1_ref) & (E_ref < E2_ref);
    U_ref(mask2) = cal_ref.U1 + cal_ref.G2 * sqrt(E_ref(mask2) - E1_ref);
    mask3 = (E_ref >= E2_ref);
    U_ref(mask3) = cal_ref.U1 + cal_ref.G2 * sqrt(E_ref(mask3) - E1_ref);
    U_ref(E_ref < E0) = cal_ref.U0;
    
    % Temperature correction (optional)
    if apply_T_correction
        T_w = 200 + cal.T_ref;
        E1_corr = E1 .* sqrt((T_w - cal.T_ref) ./ (T_w - T));
        E2_corr = E2 .* sqrt((T_w - cal.T_ref) ./ (T_w - T));
        E3_corr = E3 .* sqrt((T_w - cal.T_ref) ./ (T_w - T));
    else
        E1_corr = E1;
        E2_corr = E2;
        E3_corr = E3;
    end
    
    % Linearization
    Ucal1 = cal.C0_1 + cal.C1_1*E1_corr + cal.C2_1*E1_corr.^2 + ...
            cal.C3_1*E1_corr.^3 + cal.C4_1*E1_corr.^4;
    Ucal2 = cal.C0_2 + cal.C1_2*E2_corr + cal.C2_2*E2_corr.^2 + ...
            cal.C3_2*E2_corr.^3 + cal.C4_2*E2_corr.^4;
    Ucal3 = cal.C0_3 + cal.C1_3*E3_corr + cal.C2_3*E3_corr.^2 + ...
            cal.C3_3*E3_corr.^3 + cal.C4_3*E3_corr.^4;
    
    % Decomposition into wire velocities
    cos_angle = cosd(35.3);
    A = [cal.k1_sq, 1, cal.h1_sq;
         cal.h2_sq, cal.k2_sq, 1;
         1, cal.h3_sq, cal.k3_sq];
    
    U1_sq = zeros(size(Ucal1));
    U2_sq = zeros(size(Ucal2));
    U3_sq = zeros(size(Ucal3));
    
    for i = 1:length(Ucal1)
        B = [Ucal1(i)^2 * (1 + cal.k1_sq + cal.h1_sq).^2 * cos_angle^2;
             Ucal2(i)^2 * (1 + cal.k2_sq + cal.h2_sq).^2 * cos_angle^2;
             Ucal3(i)^2 * (1 + cal.k3_sq + cal.h3_sq).^2 * cos_angle^2];
        
        U_sq = A \ B;
        U1_sq(i) = max(U_sq(1), 0);
        U2_sq(i) = max(U_sq(2), 0);
        U3_sq(i) = max(U_sq(3), 0);
    end
    
    U1 = sqrt(U1_sq);
    U2 = sqrt(U2_sq);
    U3 = sqrt(U3_sq);
    
    % Transform to probe coordinates
    U = U1 * cosd(54.74) + U2 * cosd(54.74) + U3 * cosd(54.74);
    V = -U1 * cosd(45) - U2 * cosd(135) + U3 * cosd(90);
    W = -U1 * cosd(114.09) - U2 * cosd(114.09) - U3 * cosd(35.26);
end



%% Function: Plot vertical profiles
function plot_vertical_profiles(profile_data, heights)
    figure('Position', [100, 100, 1200, 800]);
    
    % Streamwise velocity profile
    subplot(2,3,1);
    plot(profile_data.U_mean, heights, 'b-o', 'LineWidth', 2, 'MarkerSize', 8);
    xlabel('U (m/s)', 'FontSize', 12);
    ylabel('Height (mm)', 'FontSize', 12);
    title('Streamwise Velocity Profile', 'FontSize', 14);
    grid on;
    
    % Lateral velocity profile
    subplot(2,3,2);
    plot(profile_data.V_mean, heights, 'r-o', 'LineWidth', 2, 'MarkerSize', 8);
    xlabel('V (m/s)', 'FontSize', 12);
    ylabel('Height (mm)', 'FontSize', 12);
    title('Lateral Velocity Profile', 'FontSize', 14);
    grid on;
    
    % Vertical velocity profile
    subplot(2,3,3);
    plot(profile_data.W_mean, heights, 'g-o', 'LineWidth', 2, 'MarkerSize', 8);
    xlabel('W (m/s)', 'FontSize', 12);
    ylabel('Height (mm)', 'FontSize', 12);
    title('Vertical Velocity Profile', 'FontSize', 14);
    grid on;
    
    % Temperature profile
    subplot(2,3,4);
    plot(profile_data.T_mean, heights, 'k-o', 'LineWidth', 2, 'MarkerSize', 8);
    xlabel('Temperature (°C)', 'FontSize', 12);
    ylabel('Height (mm)', 'FontSize', 12);
    title('Temperature Profile', 'FontSize', 14);
    grid on;
    
    % Velocity magnitude profile
    subplot(2,3,5);
    plot(profile_data.U_mag_mean, heights, 'm-o', 'LineWidth', 2, 'MarkerSize', 8);
    xlabel('|U| (m/s)', 'FontSize', 12);
    ylabel('Height (mm)', 'FontSize', 12);
    title('Velocity Magnitude Profile', 'FontSize', 14);
    grid on;
    
    % RMS profiles
    subplot(2,3,6);
    hold on;
    plot(profile_data.U_rms, heights, 'b-o', 'LineWidth', 2, 'MarkerSize', 6, 'DisplayName', 'U_{rms}');
    plot(profile_data.V_rms, heights, 'r-o', 'LineWidth', 2, 'MarkerSize', 6, 'DisplayName', 'V_{rms}');
    plot(profile_data.W_rms, heights, 'g-o', 'LineWidth', 2, 'MarkerSize', 6, 'DisplayName', 'W_{rms}');
    xlabel('RMS Velocity (m/s)', 'FontSize', 12);
    ylabel('Height (mm)', 'FontSize', 12);
    title('RMS Velocity Profiles', 'FontSize', 14);
    legend('Location', 'best');
    grid on;
    hold off;
end

%% Function: Save profile results
function save_profile_results(profile_data, heights, filename)
    results = [heights', profile_data.U_mean, profile_data.V_mean, ...
               profile_data.W_mean, profile_data.U_mag_mean, ...
               profile_data.T_mean, profile_data.U_rms, ...
               profile_data.V_rms, profile_data.W_rms];
    
    writematrix(results, filename);
    fprintf('\nResults saved to %s\n', filename);
    
    % Display summary
    fprintf('\n=== Vertical Profile Summary ===\n');
    fprintf('Height range: %.1f - %.1f mm\n', min(heights), max(heights));
    fprintf('U range: %.3f - %.3f m/s\n', min(profile_data.U_mean), max(profile_data.U_mean));
    fprintf('Temperature range: %.2f - %.2f °C\n', min(profile_data.T_mean), max(profile_data.T_mean));
end

%% Function: Plot time series for each height
function plot_time_series_by_height(data_collection)
    num_heights = length(data_collection);
    
    figure('Position', [100, 100, 1400, 900]);
    
    % U velocity time series
    subplot(3,2,1);
    hold on;
    for h = 1:num_heights
        plot(data_collection(h).time, data_collection(h).U, 'DisplayName', ...
            sprintf('z=%.0f mm', data_collection(h).height));
    end
    ylabel('U (m/s)');
    title('Streamwise Velocity Time Series');
    legend('Location', 'best');
    grid on;
    hold off;
    
    % V velocity time series
    subplot(3,2,3);
    hold on;
    for h = 1:num_heights
        plot(data_collection(h).time, data_collection(h).V, 'DisplayName', ...
            sprintf('z=%.0f mm', data_collection(h).height));
    end
    ylabel('V (m/s)');
    title('Lateral Velocity Time Series');
    legend('Location', 'best');
    grid on;
    hold off;
    
    % W velocity time series
    subplot(3,2,5);
    hold on;
    for h = 1:num_heights
        plot(data_collection(h).time, data_collection(h).W, 'DisplayName', ...
            sprintf('z=%.0f mm', data_collection(h).height));
    end
    xlabel('Time (s)');
    ylabel('W (m/s)');
    title('Vertical Velocity Time Series');
    legend('Location', 'best');
    grid on;
    hold off;
    
    % Temperature time series
    subplot(3,2,2);
    hold on;
    for h = 1:num_heights
        plot(data_collection(h).time, data_collection(h).T, 'DisplayName', ...
            sprintf('z=%.0f mm', data_collection(h).height));
    end
    ylabel('Temperature (°C)');
    title('Temperature Time Series');
    legend('Location', 'best');
    grid on;
    hold off;
    
    % Reference velocity time series
    subplot(3,2,4);
    hold on;
    for h = 1:num_heights
        plot(data_collection(h).time, data_collection(h).U_ref, 'DisplayName', ...
            sprintf('z=%.0f mm', data_collection(h).height));
    end
    ylabel('U_{ref} (m/s)');
    title('Reference Velocity Time Series');
    legend('Location', 'best');
    grid on;
    hold off;
    
    % Velocity magnitude time series
    subplot(3,2,6);
    hold on;
    for h = 1:num_heights
        plot(data_collection(h).time, data_collection(h).U_mag, 'DisplayName', ...
            sprintf('z=%.0f mm', data_collection(h).height));
    end
    xlabel('Time (s)');
    ylabel('|U| (m/s)');
    title('Velocity Magnitude Time Series');
    legend('Location', 'best');
    grid on;
    hold off;
end