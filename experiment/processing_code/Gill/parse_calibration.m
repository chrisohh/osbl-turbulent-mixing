%% Parse Calibration File
function cal_data = parse_calibration(filename)
    % Read the entire file
    fid = fopen(filename, 'r');
    text = fread(fid, '*char')';
    fclose(fid);
    
    % Initialize structure
    cal_data = struct();
    
    % Extract calibration coefficients for each sensor
    for sensor = 1:3
        % Corrected regex pattern that handles newlines
        pattern = sprintf(['Probe sensor no\\.:\\s*%d[\\s\\S]*?' ...
                  'C0:\\s*([\\-\\d\\.]+)[\\s\\S]*?' ...
                  'C1:\\s*([\\-\\d\\.]+)[\\s\\S]*?' ...
                  'C2:\\s*([\\-\\d\\.]+)[\\s\\S]*?' ...
                  'C3:\\s*([\\-\\d\\.]+)[\\s\\S]*?' ...
                  'C4:\\s*([\\-\\d\\.]+)'], sensor);
        tokens = regexp(text, pattern, 'tokens', 'once', 'dotexceptnewline');
        
        if ~isempty(tokens)
            cal_data.(sprintf('C0_%d', sensor)) = str2double(tokens{1});
            cal_data.(sprintf('C1_%d', sensor)) = str2double(tokens{2});
            cal_data.(sprintf('C2_%d', sensor)) = str2double(tokens{3});
            cal_data.(sprintf('C3_%d', sensor)) = str2double(tokens{4});
            cal_data.(sprintf('C4_%d', sensor)) = str2double(tokens{5});
        end
    end
    
    % Extract yaw (k) and pitch (h) constants
    k_line = regexp(text, 'k\s*=\s*(.+)', 'tokens', 'once', 'lineanchors');
    k_numbers = regexp(k_line{1}, '[\d\.]+', 'match');
    cal_data.k1_sq = str2double(k_numbers{1})^2;
    cal_data.k2_sq = str2double(k_numbers{2})^2;
    cal_data.k3_sq = str2double(k_numbers{3})^2;
    fprintf('Found k values: %.6f, %.6f, %.6f\n', ...
        str2double(k_numbers{1}), str2double(k_numbers{2}), str2double(k_numbers{3}));

    h_line = regexp(text, 'h\s*=\s*(.+)', 'tokens', 'once', 'lineanchors');
    % Extract all numbers from the line
    h_numbers = regexp(h_line{1}, '[\d\.]+', 'match');
    cal_data.h1_sq = str2double(h_numbers{1})^2;
    cal_data.h2_sq = str2double(h_numbers{2})^2;
    cal_data.h3_sq = str2double(h_numbers{3})^2;
    fprintf('Found h values: %.6f, %.6f, %.6f\n', ...
        str2double(h_numbers{1}), str2double(h_numbers{2}), str2double(h_numbers{3}));
    
    % Extract reference temperature
    temp_pattern = 'Cal\\. ref\\. temp\\..*?([\\d\\.]+)';
    temp_tokens = regexp(text, temp_pattern, 'tokens', 'once');
    if ~isempty(temp_tokens)
        cal_data.T_ref = str2double(temp_tokens{1});
    else
        cal_data.T_ref = 20; % Default
    end
    
    % Extract probe transformation matrix
    matrix_pattern = 'Probe transformation matrix Mp\\s*([\\-\\d\\.eE]+)\\s*([\\-\\d\\.eE]+)\\s*([\\-\\d\\.eE]+)\\s*([\\-\\d\\.eE]+)\\s*([\\-\\d\\.eE]+)\\s*([\\-\\d\\.eE]+)\\s*([\\-\\d\\.eE]+)\\s*([\\-\\d\\.eE]+)\\s*([\\-\\d\\.eE]+)';
    matrix_tokens = regexp(text, matrix_pattern, 'tokens', 'once');
    
    if ~isempty(matrix_tokens)
        cal_data.Mp = [str2double(matrix_tokens{1}), str2double(matrix_tokens{2}), str2double(matrix_tokens{3});
                       str2double(matrix_tokens{4}), str2double(matrix_tokens{5}), str2double(matrix_tokens{6});
                       str2double(matrix_tokens{7}), str2double(matrix_tokens{8}), str2double(matrix_tokens{9})];
    end
end