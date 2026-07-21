function field_wf = transform_field_to_wavefollowing(field_cart, z_cart, transfo)
%TRANSFORM_FIELD_TO_WAVEFOLLOWING  Interpolate Cartesian field onto wave-following grid.
%
%   Adapted from TransformVelField_decay.m for water side.
%   For each x-column, interpolates from physical z to constant-zeta levels
%   using the mapping z(x, zeta) from generate_transfo_water.
%
%   Inputs:
%     field_cart  [Nz_cart x Nx]  field on Cartesian (x, z) grid
%     z_cart      [Nz_cart x 1]   physical z coordinates (m, negative downward,
%                                  0 at mean surface)
%     transfo     struct from generate_transfo_water
%
%   Output:
%     field_wf    [Nz_wf x Nx]   field on wave-following (x, zeta) grid

Z_grid = transfo.Z_grid;   % [Nz_wf x Nx] physical z at each (zeta, x)
[Nz_wf, Nx] = size(Z_grid);

z_cart = z_cart(:);  % column, negative downward
field_wf = NaN(Nz_wf, Nx);

for ix = 1:Nx
    col = field_cart(:, ix);

    % Skip columns with no valid data
    valid = ~isnan(col);
    if sum(valid) < 3
        continue;
    end

    z_valid = z_cart(valid);
    col_valid = col(valid);

    % Target z values for this column (from the wave-following grid)
    z_target = Z_grid(:, ix);

    % Interpolate: spline within data range, NaN outside
    in_range = (z_target >= min(z_valid)) & (z_target <= max(z_valid));

    if sum(in_range) < 1
        continue;
    end

    % Use spline for smooth interpolation within range
    field_wf(in_range, ix) = interp1(z_valid, col_valid, z_target(in_range), 'spline');

    % Extrapolate one grid point near surface (zeta ~ 0) if needed
    first_valid = find(in_range, 1, 'last');  % closest to surface (z=0)
    if first_valid < Nz_wf
        next = first_valid + 1;
        if next <= Nz_wf
            field_wf(next, ix) = interp1(z_valid, col_valid, z_target(next), 'spline');
        end
    end
end

end
