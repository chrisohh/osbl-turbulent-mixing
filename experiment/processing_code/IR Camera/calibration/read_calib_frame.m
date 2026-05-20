function frame = read_calib_frame(path, idx)
% READ_CALIB_FRAME  Load one calibration frame from .ats or synthetic .mat.
%
%   frame = read_calib_frame(path, idx)
%
% Returns a single-precision frame. For .ats inputs the FlirMovieReader is
% cached per-path in a persistent map so repeated calls within one
% calibration run don't pay the reset-and-step cost from readFrame.m.
% For .mat inputs, the file must contain a variable `frames` of size
% [H W N] (single).

    [~,~,ext] = fileparts(path);
    switch lower(ext)
        case '.ats'
            persistent readers
            if isempty(readers)
                readers = containers.Map('KeyType','char','ValueType','any');
            end
            if ~isKey(readers, path)
                addpath('C:\Program Files\FLIR Systems\sdks\file\bin\Release');
                v = FlirMovieReader(path);
                v.unit = 'temperatureFactory';
                readers(path) = v;
            end
            v = readers(path);
            reset(v);
            for i = 1:(idx-1)
                step(v);
            end
            frame = single(step(v));
        case '.mat'
            s = matfile(path);
            frame = single(s.frames(:,:,idx));
        otherwise
            error('read_calib_frame:badExt', ...
                'Unsupported extension "%s" (expected .ats or .mat)', ext);
    end
end
