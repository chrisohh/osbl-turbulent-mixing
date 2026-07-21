function [w_wave, u_wave, wave_info] = compute_wave_orbital_velocities(eta, x, t, z, params)

% Dimensions
Nx = length(x);
Nt = length(t);
Nz = length(z);

g     = params.g;
gamma = params.gamma;

dx = mean(diff(x));

% Preallocate
a      = zeros(Nx, Nt);
phi    = zeros(Nx, Nt);
k      = zeros(Nx, Nt);
omega  = zeros(Nx, Nt);

for it = 1:Nt

    eta_xt = eta(:, it);

    % --- Hilbert transform ---
    eta_hat = hilbert(eta_xt);

    a(:, it)   = abs(eta_hat);
    phi_raw    = unwrap(angle(eta_hat));

    % --- Smooth phase (VERY IMPORTANT) ---
    phi(:, it) = smoothn(phi_raw, 1e-4);   % reuse your smoothn

    % --- Wavenumber ---
    dphidx = gradient(phi(:, it), dx);

    % Smooth k to remove noise spikes
    k(:, it) = smoothn(dphidx, 1e-3);

    % Enforce physically reasonable minimum wavenumber.
    % Capillary waves here: k ~ 100-300 m^-1 (lambda ~ 2-6 cm).
    % Floor at 20 rad/m (lambda ~ 30 cm) rejects spurious near-zero values
    % that would produce exp(k*z) ~ 1 at all depths and inflate orbital velocities.
    k(:, it) = max(k(:, it), 20);

    % --- Dispersion relation (deep water) ---
    omega(:, it) = sqrt(k(:, it) .* (g + gamma .* k(:, it).^2));

end

% --- Allocate orbital velocities ---
w_wave = zeros(Nx, Nz, Nt);
u_wave = zeros(Nx, Nz, Nt);

for it = 1:Nt
    for ix = 1:Nx

        % Vectorize over depth: exp(k*z) for all z at once → [1 × Nz]
        decay = exp(k(ix,it) .* z(:)');

        u_wave(ix,:,it) = a(ix,it) * omega(ix,it) * cos(phi(ix,it)) .* decay;
        w_wave(ix,:,it) = a(ix,it) * omega(ix,it) * sin(phi(ix,it)) .* decay;

    end
end

% Pack info
wave_info.k       = k;
wave_info.omega   = omega;
wave_info.a       = a;
wave_info.phi     = phi;

end