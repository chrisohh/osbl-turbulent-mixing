function [v_turb, w_turb] = remove_wave_from_piv(v_meas, w_meas, method, varargin)

% v_meas, w_meas: [Nt × Ny × Nz]

[Nt, Ny, Nz] = size(w_meas);

v_turb = zeros(size(v_meas));
w_turb = zeros(size(w_meas));

switch method

    case 'hilbert'

        % Input: w_wave_profile [Nz × Nt]
        w_wave_profile = varargin{1};

        for it = 1:Nt
            for iz = 1:Nz

                w_wave_val = w_wave_profile(iz, it);

                w_turb(it,:,iz) = w_meas(it,:,iz) - w_wave_val;
                v_turb(it,:,iz) = v_meas(it,:,iz); % no correction

            end
        end

    case 'hmean'

        % Waves propagate in x (along-wind) and are uniform across y (cross-wind).
        % The y-mean of w therefore captures the wave orbital velocity.
        % Cross-wind velocity v has no wave orbital component; leave it uncorrected
        % so that compute_horizontal_avg_turbstats performs the Reynolds decomposition.
        for it = 1:Nt
            for iz = 1:Nz

                w_mean = mean(w_meas(it,:,iz), 2);

                w_turb(it,:,iz) = w_meas(it,:,iz) - w_mean;
                v_turb(it,:,iz) = v_meas(it,:,iz);  % no wave correction for v

            end
        end

    otherwise
        error('Unknown method');

end

end