function [Sx, Sy] = load_slope_frame(cache_dir, n)
% LOAD_SLOPE_FRAME  Read one cached slope frame from disk.

frame_path = fullfile(cache_dir, sprintf('frame_%04d.mat', n));
S = load(frame_path, 'Sx', 'Sy', 't');
Sx = S.Sx;
Sy = S.Sy;
end
