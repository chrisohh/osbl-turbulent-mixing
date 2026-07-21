function [v_zeta, w_zeta, zeta_grid] = apply_vertical_shift_transverse(v_meas, w_meas, z_piv, eta_t, zeta_grid)
%APPLY_VERTICAL_SHIFT_TRANSVERSE  Surface-following coords for transverse PIV.
%
%   In the transverse (y-z) plane, the wave surface is uniform across y:
%   eta(y,t) ~ eta(t).  The full curvilinear mapping degenerates to a
%   simple vertical shift at each time step:
%
%       zeta = z - eta(t)     (depth below instantaneous surface)
%
%   Then interpolate onto a uniform zeta grid so that zeta=0 always
%   corresponds to the surface.  (See LaTeX doc Section 4.2.)
%
%   Inputs:
%     v_meas   [Nt x Ny x Nz]  cross-wind velocity in Cartesian z
%     w_meas   [Nt x Ny x Nz]  vertical velocity in Cartesian z
%     z_piv    [Nz x 1]        Cartesian z, positive downward (depth)
%     eta_t    [Nt x 1]        surface elevation at transverse station
%                               (positive = crest, metres)
%     zeta_grid [Nzeta x 1]   (optional) target uniform zeta grid.
%                               Default: same spacing as z_piv, starting at 0.
%
%   Outputs:
%     v_zeta   [Nt x Ny x Nzeta]  v on uniform zeta grid
%     w_zeta   [Nt x Ny x Nzeta]  w on uniform zeta grid
%     zeta_grid [Nzeta x 1]       the zeta coordinates used

[Nt, Ny, Nz] = size(v_meas);
z_piv  = z_piv(:);
eta_t  = eta_t(:);

% --- Default zeta grid: same spacing, starting at 0 ---
if nargin < 5 || isempty(zeta_grid)
    dz = abs(z_piv(2) - z_piv(1));
    % Max depth below surface that we can always resolve
    % (deepest z minus highest crest)
    zeta_max = z_piv(end) + min(eta_t);   % both positive-down
    zeta_grid = (0 : dz : zeta_max)';
end

Nzeta  = length(zeta_grid);
v_zeta = NaN(Nt, Ny, Nzeta);
w_zeta = NaN(Nt, Ny, Nzeta);

for it = 1:Nt
    % z_shifted = z_piv - (-eta_t(it))
    % In positive-down convention:
    %   eta > 0 (crest) → surface above mean → shift z down by eta
    %   zeta = z + eta(t)   (depth below instantaneous surface)
    %
    % Physical z of the surface = -eta_t(it) in positive-down coords,
    % so depth below surface = z_piv - (-eta_t(it)) = z_piv + eta_t(it).
    z_shifted = z_piv + eta_t(it);   % zeta values of each Cartesian z point

    for iy = 1:Ny
        v_col = squeeze(v_meas(it, iy, :));
        w_col = squeeze(w_meas(it, iy, :));

        valid = ~isnan(v_col);
        if sum(valid) < 3, continue; end

        % Interpolate from (z_shifted) onto uniform (zeta_grid)
        in_range = zeta_grid >= min(z_shifted(valid)) & ...
                   zeta_grid <= max(z_shifted(valid));
        if sum(in_range) < 1, continue; end

        v_zeta(it, iy, in_range) = interp1(z_shifted(valid), v_col(valid), ...
                                            zeta_grid(in_range), 'spline');
        w_zeta(it, iy, in_range) = interp1(z_shifted(valid), w_col(valid), ...
                                            zeta_grid(in_range), 'spline');
    end
end

end
