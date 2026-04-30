%%
% Method 1: Check file size
file_info = dir(strcat(obs_folderName,obs_fileName));
file_bytes = file_info.bytes;

% For 12MP camera (4000 x 3000):
expected_bayer_8bit = 4000 * 3000 * 1;      % 12 MB
expected_rgb_8bit = 4000 * 3000 * 3;        % 36 MB
expected_bayer_16bit = 4000 * 3000 * 2;     % 24 MB
expected_rgb_16bit = 4000 * 3000 * 3 * 2;   % 72 MB

fprintf('File size: %d bytes (%.1f MB)\n', file_bytes, file_bytes/1e6);
fprintf('Expected Bayer 8-bit: %.1f MB\n', expected_bayer_8bit/1e6);
fprintf('Expected RGB 8-bit: %.1f MB\n', expected_rgb_8bit/1e6);

% Method 2: Try loading as binary and check
fid = fopen(strcat(obs_folderName,obs_fileName), 'r');
data = fread(fid, 100, 'uint8');  % Read first 100 bytes
fclose(fid);
fprintf('First 20 values: '); disp(data(1:20)');

% Method 3: Try loading and reshaping
fid = fopen(strcat(obs_folderName,obs_fileName), 'r');
all_data = fread(fid, '*uint8');
fclose(fid);
fprintf('Total bytes read: %d\n', length(all_data));

% Try as RGB
if mod(length(all_data), 3) == 0
    fprintf('Could be RGB (divisible by 3)\n');
    % Try reshape
    try
        rgb_test = reshape(all_data, 3, 4000, 3000);
        rgb_test = permute(rgb_test, [3 2 1]);  % [height, width, 3]
        fprintf('Successfully reshaped as RGB!\n');
        figure; imshow(rgb_test);
    catch
        fprintf('Failed to reshape as RGB\n');
    end
end
%%
%% Inspect Raw File Header
% Examine the header bytes to understand the file format
filename = strcat(ref_folderName,ref_fileName);

fprintf('===== HEADER INSPECTOR =====\n');
fprintf('File: %s\n\n', filename);

fid = fopen(filename, 'r');
if fid == -1
    error('Cannot open file: %s', filename);
end

%% Method 1: View as ASCII text (if header is text-based)
fprintf('===== FIRST 2000 BYTES AS TEXT =====\n');
header_text = fread(fid, 2000, '*uint8');
ascii_text = char(header_text');

% Replace non-printable characters with dots for readability
ascii_clean = ascii_text;
ascii_clean(header_text < 32 | header_text > 126) = '.';

fprintf('%s\n\n', ascii_clean);

% Check if it looks like text
printable_ratio = sum(header_text >= 32 & header_text <= 126) / length(header_text);
fprintf('Printable characters: %.1f%%\n', printable_ratio * 100);
if printable_ratio > 0.5
    fprintf('→ Header appears to be TEXT-based\n\n');
else
    fprintf('→ Header appears to be BINARY\n\n');
end

%% Method 2: View as hex dump (like a hex editor)
fprintf('===== FIRST 512 BYTES AS HEX =====\n');
fseek(fid, 0, 'bof');
hex_bytes = fread(fid, 512, '*uint8');

% Print in hex dump format (16 bytes per line)
for i = 1:16:min(256, length(hex_bytes))
    % Address
    fprintf('%08X: ', i-1);
    
    % Hex values
    for j = i:min(i+15, length(hex_bytes))
        fprintf('%02X ', hex_bytes(j));
    end
    
    % Padding
    if i+15 > length(hex_bytes)
        fprintf(repmat('   ', 1, i+15-length(hex_bytes)));
    end
    
    % ASCII representation
    fprintf(' | ');
    for j = i:min(i+15, length(hex_bytes))
        if hex_bytes(j) >= 32 && hex_bytes(j) <= 126
            fprintf('%c', hex_bytes(j));
        else
            fprintf('.');
        end
    end
    fprintf('\n');
end
fprintf('\n');

%% Method 3: Look for common patterns/markers
fprintf('===== SEARCHING FOR PATTERNS =====\n');
fseek(fid, 0, 'bof');
search_bytes = fread(fid, 10000, '*uint8');

% Look for common strings
search_strings = {'width', 'height', 'size', 'format', 'RGB', 'Bayer', ...
                  'bits', 'depth', 'data', 'offset', 'header', 'camera', ...
                  'image', 'resolution', 'pixel', 'END', '---', 'START'};

fprintf('Looking for keywords...\n');
for k = 1:length(search_strings)
    str = search_strings{k};
    % Search case-insensitive
    idx = strfind(lower(char(search_bytes')), lower(str));
    if ~isempty(idx)
        fprintf('  Found "%s" at byte %d\n', str, idx(1));
        % Show context (50 chars before and after)
        context_start = max(1, idx(1)-50);
        context_end = min(length(search_bytes), idx(1)+50);
        context = char(search_bytes(context_start:context_end)');
        context(search_bytes(context_start:context_end) < 32) = '.';
        fprintf('    Context: ...%s...\n', context);
    end
end
fprintf('\n');

%% Method 4: Check for common binary patterns
fprintf('===== BINARY STRUCTURE ANALYSIS =====\n');

% Check if first bytes are a known file signature
fseek(fid, 0, 'bof');
sig = fread(fid, 4, '*uint8');
fprintf('First 4 bytes: %02X %02X %02X %02X\n', sig(1), sig(2), sig(3), sig(4));

% Common signatures:
if all(sig == [255 216 255 224]') || all(sig == [255 216 255 225]')
    fprintf('→ This is a JPEG file!\n');
elseif all(sig == [137 80 78 71]')
    fprintf('→ This is a PNG file!\n');
elseif all(sig == [73 73 42 0]') || all(sig == [77 77 0 42]')
    fprintf('→ This is a TIFF file!\n');
else
    fprintf('→ Not a standard image format (custom raw)\n');
end

% Look for null-terminated strings
fseek(fid, 0, 'bof');
all_header = fread(fid, min(5000, 1748764), '*uint8');
null_positions = find(all_header == 0);
if ~isempty(null_positions)
    fprintf('\nFound %d null bytes (0x00) in first 5000 bytes\n', length(null_positions));
    fprintf('First null at byte %d\n', null_positions(1));
end

% Check for repeated patterns (might indicate structure)
fseek(fid, 0, 'bof');
pattern_bytes = fread(fid, 1000, '*uint8');
[unique_vals, ~, idx] = unique(pattern_bytes);
fprintf('\nFirst 1000 bytes contain %d unique values\n', length(unique_vals));
most_common = mode(pattern_bytes);
fprintf('Most common byte value: %d (0x%02X) appears %d times\n', ...
    most_common, most_common, sum(pattern_bytes == most_common));

%% Method 5: Try to parse as structured data
fprintf('\n===== ATTEMPTING STRUCTURED PARSING =====\n');
fseek(fid, 0, 'bof');

% Try reading as 32-bit integers
int32_values = fread(fid, 20, 'int32');
fprintf('First 20 bytes as int32: ');
fprintf('%d ', int32_values(1:5));
fprintf('...\n');

% Try reading as floats
fseek(fid, 0, 'bof');
float_values = fread(fid, 20, 'float32');
fprintf('First 20 bytes as float32: ');
fprintf('%.2f ', float_values(1:5));
fprintf('...\n');

fclose(fid);

fprintf('\n===== RECOMMENDATIONS =====\n');
fprintf('1. Look at the TEXT output above - can you read any metadata?\n');
fprintf('2. Check HEX dump for recognizable patterns\n');
fprintf('3. If you see strings like "width=4000", the format is probably text-based\n');
fprintf('4. The header size is probably 1748764 bytes\n');
fprintf('5. Try opening in ImageJ: File -> Import -> Raw, skip first 1748764 bytes\n');