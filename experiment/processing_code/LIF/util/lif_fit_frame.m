function rgb = lif_fit_frame(rgb, vsz)
% LIF_FIT_FRAME  Force an RGB frame to the size a VideoWriter was opened with.
%   rgb = lif_fit_frame(rgb, vsz) crops or white-pads rgb to vsz = [rows cols].
%   print('-RGBImage') can return a figure a pixel or two larger/smaller from
%   one call to the next (axis labels changing width, DPI rounding), and
%   VideoWriter locks the frame size on the first writeVideo, so every later
%   frame has to be nudged back to it.

    [ny, nx, ~] = size(rgb);

    if ny > vsz(1), rgb = rgb(1:vsz(1), :, :); end
    if nx > vsz(2), rgb = rgb(:, 1:vsz(2), :); end

    [ny, nx, ~] = size(rgb);
    if ny < vsz(1) || nx < vsz(2)
        padded = repmat(uint8(255), [vsz(1) vsz(2) size(rgb, 3)]);   % figure is white
        padded(1:ny, 1:nx, :) = rgb;
        rgb = padded;
    end
end
