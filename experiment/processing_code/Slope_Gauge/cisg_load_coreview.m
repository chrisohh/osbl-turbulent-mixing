function [rgb_image, metadata] = cisg_load_coreview(filename)
% CISG_LOAD_COREVIEW Load CoreView/Flare camera raw files
%
% This format has a 20-byte binary header:
%   Bytes 0-3:   Width (int32)
%   Bytes 4-7:   Height (int32)
%   Bytes 8-11:  Bits per pixel (int32)
%   Bytes 12-15: Unknown (int32)
%   Bytes 16-19: Unknown (int32)
%   Then: RGB image data (width × height × 3 bytes)
%
% Input:
%   filename - path to .raw file
%
% Output:
%   rgb_image - RGB image [height x width x 3], uint8
%   metadata - struct with header information

    % Open file
    fid = fopen(filename, 'r');
    if fid == -1
        error('Cannot open file: %s', filename);
    end
    
    % Read binary header
    width = fread(fid, 1, 'int32');
    height = fread(fid, 1, 'int32');
    bits_per_pixel = fread(fid, 1, 'int32');
    skipframe = fread(fid, 1, 'int32');
    skipsomething = fread(fid, 1, 'int32');% check the export settings
    
    % Store metadata
    metadata = struct();
    metadata.width = width;
    metadata.height = height;
    metadata.bits_per_pixel = bits_per_pixel;
    metadata.header_bytes = 20;
    
    % Calculate expected data size
    bytes_per_pixel = bits_per_pixel / 8;
    expected_data = width * height * bytes_per_pixel;
    
    % Determine data type
    if bits_per_pixel <= 8
        data_type = 'uint8';
        metadata.bit_depth = 8;
    elseif bits_per_pixel <= 16
        data_type = 'uint8';  % Still read as uint8 since it's 24-bit (3×8)
        metadata.bit_depth = 8;
    else
        data_type = 'uint16';
        metadata.bit_depth = 16;
    end
    
    % Read RGB data (starts at byte 20)
    % The data after the header appears to be image data based on pattern
    % But there's more header data between byte 20 and actual image
    
    % Get file size to calculate actual image data start
    fseek(fid, 0, 'eof');
    file_size = ftell(fid);
    actual_data_size = width * height * 3;  % RGB
    header_size = file_size - actual_data_size;
    
    % Seek to start of image data
    fseek(fid, header_size, 'bof');
    
    % Read RGB data
    rgb_data = fread(fid, actual_data_size, '*uint8');
    fclose(fid);
    
    if length(rgb_data) ~= actual_data_size
        warning('Read %d bytes but expected %d', length(rgb_data), actual_data_size);
    end
    
    % Reshape to RGB image
    % Data format: interleaved, but channel order may vary
    rgb_image = reshape(rgb_data, 3, width, height);
    rgb_image = permute(rgb_image, [3 2 1]);  % [height x width x 3]
    
    % CHANNEL ORDER FIX
    % If colors are wrong, the camera might store in BGR or GRB order
    % Based on your description: red→green, blue→red, green→blue
    % This suggests GRB order (channel 1=Green, 2=Red, 3=Blue)
    % Reorder to RGB: [R G B] = [channel2 channel1 channel3]
    rgb_image = rgb_image(:,:,[2 3 1]);  % Swap first two channels (GBR → RGB)
    
    fprintf('  RGB image: %d x %d x 3\n', size(rgb_image,1), size(rgb_image,2));
    fprintf('  Loaded successfully\n');
    
    % Quick visualization check
    % figure; imshow(rgb_image); title(filename);
end