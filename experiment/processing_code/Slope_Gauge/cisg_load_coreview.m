function [rgb_image, metadata] = cisg_load_coreview(filename, expected_size)
% CISG_LOAD_COREVIEW Load CoreView/Flare camera raw files
%
% This format has a binary header (the file is header + width*height*3 bytes):
%   Bytes 0-3:   Width (int32)
%   Bytes 4-7:   Height (int32)
%   Bytes 8-11:  Bits per pixel (int32)
%   Bytes 12-15: Unknown (int32)
%   Bytes 16-19: Unknown (int32)
%   Then: RGB image data (width × height × 3 bytes), read from the END of the file
%
% Input:
%   filename      - path to .raw file
%   expected_size - optional [height width] (or struct with fields height, width) of a
%                   good frame of this run.  Enables the two checks below.
%
% Output:
%   rgb_image - RGB image [height x width x 3], uint8
%   metadata  - struct with header information (header_missing = true if the file had
%               no header)
%
% Bad files (seen in run 95): 20-byte files (header only), and files of exactly
% height*width*3 bytes (pixels without the header).  With expected_size the second
% kind is read from offset 0; everything else invalid raises 'cisg:badFrame'.

    if nargin < 2, expected_size = []; end
    if isstruct(expected_size), expected_size = [expected_size.height expected_size.width]; end

    d = dir(filename);
    if isempty(d)
        error('cisg:badFrame', 'Cannot open file: %s', filename);
    end
    file_size = d.bytes;

    metadata = struct();
    metadata.header_missing = false;
    metadata.header_bytes = 20;

    % Pixels without header: size is exactly h*w*3 of a good frame
    if ~isempty(expected_size) && file_size == prod(expected_size) * 3
        height = expected_size(1);  width = expected_size(2);
        metadata.width = width;  metadata.height = height;
        metadata.bits_per_pixel = 24;  metadata.bit_depth = 8;
        metadata.header_missing = true;  metadata.header_bytes = 0;
        rgb_data = read_bytes(filename, 0, width * height * 3);
    else
        fid = fopen(filename, 'r');
        if fid == -1
            error('cisg:badFrame', 'Cannot open file: %s', filename);
        end
        width = fread(fid, 1, 'int32');
        height = fread(fid, 1, 'int32');
        bits_per_pixel = fread(fid, 1, 'int32');
        fread(fid, 2, 'int32');                              % unknown header words (export settings)
        fclose(fid);
        ok = ~isempty(width) && ~isempty(height) && width > 0 && height > 0 && width <= 20000 && height <= 20000 ...
             && file_size >= width * height * 3;
        if ok && ~isempty(expected_size) && (height ~= expected_size(1) || width ~= expected_size(2))
            ok = false;
        end
        if ~ok
            error('cisg:badFrame', 'Bad CoreView frame (%d bytes, header w=%s h=%s): %s', file_size, ...
                  num2str(width), num2str(height), filename);
        end

        % Store metadata
        metadata.width = width;
        metadata.height = height;
        metadata.bits_per_pixel = bits_per_pixel;
        if bits_per_pixel <= 16
            metadata.bit_depth = 8;     % 24-bit (3×8) data read as uint8
        else
            metadata.bit_depth = 16;
        end

        % Image data sits at the end of the file
        rgb_data = read_bytes(filename, file_size - width * height * 3, width * height * 3);
    end

    if length(rgb_data) ~= width * height * 3
        error('cisg:badFrame', 'Read %d bytes but expected %d: %s', length(rgb_data), width * height * 3, filename);
    end

    % Reshape to RGB image (interleaved)
    rgb_image = reshape(rgb_data, 3, width, height);
    rgb_image = permute(rgb_image, [3 2 1]);  % [height x width x 3]

    % CHANNEL ORDER FIX
    % The camera stores GBR; reorder to RGB
    if ~metadata.header_missing
        rgb_image = rgb_image(:,:,[2 3 1]);
    end                                  % headerless files already store R,G,B (checked against neighbouring frames of run 95)
end

function data = read_bytes(filename, offset, n)
    fid = fopen(filename, 'r');
    if fid == -1
        error('cisg:badFrame', 'Cannot open file: %s', filename);
    end
    fseek(fid, offset, 'bof');
    data = fread(fid, n, '*uint8');
    fclose(fid);
end
