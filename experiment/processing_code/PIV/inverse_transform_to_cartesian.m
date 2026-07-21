function field_cart = inverse_transform_to_cartesian(field_wf, z_cart, transfo)
%INVERSE_TRANSFORM_TO_CARTESIAN  Map wave-following field back to Cartesian grid.
%
%   Adapted from inverseTransformVelField_decay.m for water side.
%
%   Inputs:
%     field_wf   [Nz_wf x Nx]  field on wave-following (x, zeta) grid
%     z_cart     [Nz_cart x 1]  physical z coordinates (m, negative downward)
%     transfo    struct from generate_transfo_water
%
%   Output:
%     field_cart  [Nz_cart x Nx]  field on Cartesian (x, z) grid

Z_grid = transfo.Z_grid;   % [Nz_wf x Nx]  z(zeta, x)
[Nz_wf, Nx] = size(Z_grid);
Nz_cart = length(z_cart);

z_cart = z_cart(:);
field_cart = NaN(Nz_cart, Nx);

for ix = 1:Nx
    col = field_wf(:, ix);
    valid = ~isnan(col);
    if sum(valid) < 3
        continue;
    end

    % "Initial sites": physical z at the valid zeta levels
    z_wf = Z_grid(valid, ix);
    col_valid = col(valid);

    % Target: regular Cartesian z points within the data range
    in_range = (z_cart >= min(z_wf)) & (z_cart <= max(z_wf));

    if sum(in_range) < 1
        continue;
    end

    field_cart(in_range, ix) = interp1(z_wf, col_valid, z_cart(in_range), 'spline');
end

end
