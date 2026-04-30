function [rgb_image, metadata] = cisg_load_raw_simple(filename, width, height, bit_depth, bayer_pattern, header_bytes)
% CISG_LOAD_RAW_SIMPLE Simple raw loader with manual parameters
%
% Use this if your raw format doesn't have parseable headers
%
% Inputs:
%   filename - path to raw file
%   width - image width in pixels (default: 4000)
%   height - image height in pixels (default: 3000)
%   bit_depth - bits per pixel: 8, 10, 12, or 16 (default: 8)
%   bayer_pattern - 'rggb', 'bggr', 'grbg', or 'gbrg' (default: 'rggb')
%   header_bytes - number of header bytes to skip (default: 0)
%
% Output:
%   rgb_image - demosaiced RGB image
%   metadata - struct with parameters used

% Set defaults
if nargin < 2 || isempty(width), width = 4000; end
if nargin < 3 || isempty(height), height = 3000; end
if nargin < 4 || isempty(bit_depth), bit_depth = 8; end
if nargin < 5 || isempty(bayer_pattern), bayer_pattern = 'rggb'; end
if nargin < 6 || isempty(header_bytes), header_bytes = 0; end

% Store metadata
metadata = struct();
metadata.width = width;
metadata.height = height;
metadata.bit_depth = bit_depth;
metadata.bayer_pattern = lower(bayer_pattern);
metadata.header_bytes = header_bytes;

% Determine data type
if bit_depth <= 8
    data_type = 'uint8';
else
    data_type = 'uint16';
end

% Open and read file
fid = fopen(filename, 'r');
if fid == -1
    error('Cannot open file: %s', filename);
end

% Skip header
if header_bytes > 0
    fseek(fid, header_bytes, 'bof');
end

% Read raw data
raw_data = fread(fid, width * height, ['*' data_type]);
fclose(fid);

% Check if we got expected amount of data
if length(raw_data) ~= width * height
    warning('Read %d pixels but expected %d', length(raw_data), width*height);
end

% Reshape to 2D Bayer pattern
raw_bayer = reshape(raw_data, width, height)';

fprintf('Loaded raw: %dx%d, %d-bit, %s\n', ...
    width, height, bit_depth, upper(bayer_pattern));

% Demosaic
rgb_image = demosaic(raw_bayer, metadata.bayer_pattern);

fprintf('Demosaiced to RGB\n');
end