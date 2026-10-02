function [z0, z_cm, dy] = lif_surface_datum(Ibg, z_cm, dx_mm, dy)
%LIF_SURFACE_DATUM Put z = 0 at the water level of a reference frame.
%
%   [z0, z_cm, dy] = LIF_SURFACE_DATUM(Ibg, z_cm, dx_mm, dy)
%
% Ibg   - reference frame (the first, dye-free frame), same grid as z_cm
% z_cm  - [1 x ny] row heights (cm) in the old datum (board centre)
% dx_mm - pixel size (mm)
% dy    - lif_dye_front options; surf_window_cm is in the old datum
%
% The surface is found exactly as in lif_dye_front (darkest smoothed row
% per column) and averaged over x -> z0 (old datum).  Returns z_cm - z0
% and dy with surf_window_cm shifted by -z0, so everything downstream is
% measured from the reference water level.

    [~, z_surf] = lif_dye_front(Ibg, z_cm, dx_mm, dy);
    z0   = mean(z_surf, 'omitnan');
    z_cm = z_cm - z0;
    dy.surf_window_cm = dy.surf_window_cm - z0;
    fprintf('Datum: reference water level at %.2f cm above the board centre (x-std %.2f cm) -> z = 0\n', ...
            z0, std(z_surf, 'omitnan'));
end
