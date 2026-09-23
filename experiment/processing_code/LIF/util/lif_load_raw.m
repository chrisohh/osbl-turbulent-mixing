function [img, meta] = lif_load_raw(filename, width, height)
%LIF_LOAD_RAW Load a monochrome Flare/CoreView .raw frame (MCL cameras).
%
%   [img, meta] = LIF_LOAD_RAW(filename)               uses 4096 x 3072
%   [img, meta] = LIF_LOAD_RAW(filename, width, height)
%
% Unlike the colour (CCL) CoreView files handled by cisg_load_coreview, the
% MCL files carry NO header: the file is exactly width*height uint16 samples
% in row-major (x fastest) order.  The sensor is 10-bit, so values run 0-1023
% and are returned as double (raw counts, no scaling).
%
% Output:
%   img  - [height x width] double, raw counts
%   meta - struct with width, height, bytes, bit_depth estimate

    if nargin < 2 || isempty(width),  width  = 4096; end
    if nargin < 3 || isempty(height), height = 3072; end

    d = dir(filename);
    if isempty(d)
        error('lif_load_raw:missing', 'File not found: %s', filename);
    end

    expected = width * height * 2;
    if d.bytes ~= expected
        error('lif_load_raw:size', ...
              '%s is %d bytes; expected %d for %dx%d uint16.', ...
              filename, d.bytes, expected, width, height);
    end

    fid = fopen(filename, 'r');
    if fid == -1
        error('lif_load_raw:open', 'Cannot open file: %s', filename);
    end
    raw = fread(fid, width*height, '*uint16');
    fclose(fid);

    % x fastest -> reshape as [width height] then transpose to [height width]
    img = double(reshape(raw, width, height).');

    meta = struct('width', width, 'height', height, 'bytes', d.bytes, ...
                  'bit_depth', 10, 'file', filename);
end
