function [Sx, Sy] = load_slope_cube(cache_dir, frame_numbers, roi)
% LOAD_SLOPE_CUBE  Materialize a frame range into 3D arrays [Ny x Nx x Nt].
%
% frame_numbers : vector of raw frame indices to load (e.g. 1:2291).
% roi (optional): struct(row_min,row_max,col_min,col_max) in CACHED-frame
%                 coords. If supplied, each frame is cropped to this ROI
%                 before being stacked — big RAM win when only the inner
%                 region is needed downstream.

if isempty(frame_numbers)
    error('load_slope_cube:empty', 'frame_numbers is empty.');
end
crop = (nargin >= 3) && ~isempty(roi);

[Sx1, Sy1] = load_slope_frame(cache_dir, frame_numbers(1));
if crop
    Sx1 = Sx1(roi.row_min:roi.row_max, roi.col_min:roi.col_max);
    Sy1 = Sy1(roi.row_min:roi.row_max, roi.col_min:roi.col_max);
end
[ny, nx] = size(Sx1);
Nt = numel(frame_numbers);

Sx = zeros(ny, nx, Nt, 'single');
Sy = zeros(ny, nx, Nt, 'single');

Sx(:,:,1) = Sx1;
Sy(:,:,1) = Sy1;

tic;
for k = 2:Nt
    [Sxk, Syk] = load_slope_frame(cache_dir, frame_numbers(k));
    if crop
        Sxk = Sxk(roi.row_min:roi.row_max, roi.col_min:roi.col_max);
        Syk = Syk(roi.row_min:roi.row_max, roi.col_min:roi.col_max);
    end
    Sx(:,:,k) = Sxk;
    Sy(:,:,k) = Syk;
    if mod(k, 200) == 0
        fprintf('  loaded %d/%d frames (%.1f s elapsed)\n', k, Nt, toc);
    end
end
end
