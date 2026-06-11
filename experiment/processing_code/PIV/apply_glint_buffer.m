function mask_out = apply_glint_buffer(imSurf, glint_buffer_px)
% apply_glint_buffer  Extend surface mask downward to exclude glint/foam band.
%
%   mask_out = apply_glint_buffer(imSurf, glint_buffer_px)
%
%   Takes a FindSurfaceCapillary output struct and returns a mask in the
%   same NaN/1 convention used by ComputeVelocities_Quick_Filt_Deform_Water_dcorFilt:
%       1   = valid water pixel
%       NaN = excluded (air above surface OR glint band below surface)
%
%   Inputs:
%     imSurf          struct from FindSurfaceCapillary (needs .mask and .surfacePIVImg)
%     glint_buffer_px number of rows below the detected surface to also exclude
%
%   Example (in run_decomposition_linearWaveTheory.m):
%     GLINT_BUFFER = 40;
%     fprintf('Glint buffer: %d px = %.2f mm\n', GLINT_BUFFER, GLINT_BUFFER * DX * 1e3);
%     Mask1 = apply_glint_buffer(imSurfa, GLINT_BUFFER);
%     Mask2 = apply_glint_buffer(imSurfb, GLINT_BUFFER);
%     compVel = ComputeVelocities_Quick_Filt_Deform_Water_dcorFilt(...
%         IM_a, IM_b, Mask1, Mask2, IntrWndw, GrdSpc, dcorGate, iuod);

% Start from the existing mask (NaN/1 or 0/1 — normalise to NaN/1)
mask_out = double(imSurf.mask);
mask_out(mask_out == 0) = NaN;   % handle 0/1 input convention

[h, w] = size(mask_out);
rows    = (1:h)';                          % h×1
surf    = round(imSurf.surfacePIVImg(:)'); % 1×w

% Set glint band (rows just below the surface) to NaN
glint_zone = (rows > surf) & (rows <= min(surf + glint_buffer_px, h));
mask_out(glint_zone) = NaN;

