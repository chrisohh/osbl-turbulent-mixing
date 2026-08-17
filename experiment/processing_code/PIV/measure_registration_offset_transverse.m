function offset = measure_registration_offset_transverse(surfImg, pivImg, nPoints)
%MEASURE_REGISTRATION_OFFSET_TRANSVERSE Interactively measure the pixel
%   offset between a transverse surf-camera image and the corresponding
%   transverse PIV image, using matching wave features (e.g. a crest
%   apex) as tie points instead of a still-water calibration shot.
%
%   offset = measure_registration_offset_transverse(surfImg, pivImg)
%   offset = measure_registration_offset_transverse(surfImg, pivImg, nPoints)
%
%   Inputs:
%     surfImg  transverse surf-camera image, already rectified/rescaled
%              to PIV resolution (e.g. SurfImgToPIVDimsTransverse output,
%              or Surf.ImgScaledToPIVSmallCrop from
%              FindSurfaceCapillaryTransverse)
%     pivImg   the corresponding transverse PIV image (e.g. imgPiv,
%              rectified with RectifyPIVFrameTransverse if you're
%              checking against FindSurfaceCapillaryTransverse's mask)
%     nPoints  number of tie points to click (default 3). Use more than
%              1 so you can sanity-check that the offset is roughly
%              constant across the frame -- if it drifts with x, the
%              horizontal crop/resize needs remeasuring too, not just
%              the vertical offset.
%
%   Usage: pick a frame where a distinct wave crest (or trough, or any
%   sharp recognisable feature) is visible in both camera views. You'll
%   be shown the surf image first -- click the SAME feature nPoints
%   times if you have several usable ones, or the same feature is fine
%   if that's all you have. Then you'll be shown the PIV image -- click
%   the matching feature(s) in the same order.
%
%   Output offset.dy is the vertical registration constant: it directly
%   replaces the "-797+433" style constant in CropSurfToPIVDimsTransverse
%   (i.e. surfCropped = surf - offset.dy). offset.dx is informational --
%   it should be close to 0 if the horizontal crop start index is
%   correct; a large or inconsistent dx means the crop/resize factor
%   needs rechecking, not just the vertical offset.

if nargin < 3, nPoints = 3; end

figure('Name', sprintf('Surf image -- click %d matching wave feature(s) (e.g. a crest apex)', nPoints), ...
       'Color', 'w');
imagesc(surfImg); colormap gray; axis image;
title('Click the feature(s), one per pair, in any order you will repeat on the next image');
[surfX, surfY] = ginput(nPoints);
close(gcf);

figure('Name', 'PIV image -- click the SAME feature(s), in the SAME order', 'Color', 'w');
imagesc(pivImg); colormap gray; axis image;
title('Click the matching feature(s), same order as the surf image');
[pivX, pivY] = ginput(nPoints);
close(gcf);

dy = surfY - pivY;
dx = surfX - pivX;

offset.dy     = mean(dy);
offset.dx     = mean(dx);
offset.dy_std = std(dy);
offset.dx_std = std(dx);

fprintf('Per-point offsets (surf - piv):\n');
fprintf('  dx: '); fprintf('%.1f  ', dx); fprintf('\n');
fprintf('  dy: '); fprintf('%.1f  ', dy); fprintf('\n');
fprintf('Mean dy = %.1f px (std %.1f over %d points)\n', offset.dy, offset.dy_std, nPoints);
fprintf('  -> use in CropSurfToPIVDimsTransverse as: surfCropped = surf - %.0f;\n', offset.dy);
fprintf('Mean dx = %.1f px (std %.1f over %d points)\n', offset.dx, offset.dx_std, nPoints);
fprintf('  -> should be near 0; a large/consistent dx means the crop start index is off,\n');
fprintf('     a large std across points means the horizontal resize factor is off.\n');

if offset.dy_std > 5 || offset.dx_std > 5
    warning(['Offsets vary by more than 5 px across points -- check that the same physical ' ...
             'feature was clicked in both images, or that the resize/dewarp calibration is still correct.']);
end

end
