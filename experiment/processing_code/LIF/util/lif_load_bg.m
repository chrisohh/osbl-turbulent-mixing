function bg = lif_load_bg(fname, frames, downsample)
%LIF_LOAD_BG  Mean of several dye-free frames: the theta reference image.
%
%   bg = LIF_LOAD_BG(fname, frames, downsample)
%
%   fname      - handle, fname(n) = path of frame n
%   frames     - frame numbers to average.  EVERY one must be dye-free, so
%                keep well clear of the frame where the dye first appears.
%   downsample - row/column stride, applied after loading (as in the scripts)
%
% theta = ln(I_bg/I) inherits the shot noise of I_bg at full strength, and
% that noise is what makes the 99% integral wander when there is little
% dye: a single reference frame leaves sigma_theta ~ 1/sqrt(counts) at
% every depth.  Averaging N dye-free frames divides it by sqrt(N).  The dye
% normally arrives long after the recording starts, so there are plenty of
% frames to average.
%
% Stays double (not uint16): the mean is fractional counts.

    if isempty(frames), bg = []; return; end
    frames = frames(:).';
    bg = 0;
    for n = frames
        I  = lif_load_raw(fname(n));
        bg = bg + double(I(1:downsample:end, 1:downsample:end));   % double: uint16 sums round and saturate
    end
    bg = bg / numel(frames);
end
