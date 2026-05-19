function plot_E1_stokes_profile(data)
% E.1  Reconstructed u^S(z) profile from PIV-derived (a, k, omega).
% Time-averaged over full record (no stage windows yet).

eta   = data.eta;
x_eta = data.x_eta;
t_eta = data.t_eta;

info = piv_eta_hilbert(eta, x_eta, t_eta(:));
a     = double(info.a);          % [Nx x Nt]
k     = double(info.k);
omega = double(info.omega);

z = (0:-1:-100)' * 1e-3;          % m, surface to -100 mm
Nz = length(z);

% u_S(x,z,t) = a^2 * omega * k * exp(2 k z); average over x and t
uS_z = zeros(Nz, 1);
for iz = 1:Nz
    field = a.^2 .* omega .* k .* exp(2 .* k .* z(iz));
    uS_z(iz) = mean(field(:), 'omitnan');
end

% Monochromatic peak overlay
k_p   = median(k(:), 'omitnan');
om_p  = sqrt(9.81 * k_p + 7.4e-5 * k_p^3);
a_p   = sqrt(2 * mean(a(:).^2, 'omitnan'));   % rms amplitude
uS_mono = a_p^2 * om_p * k_p * exp(2 * k_p * z);

fig = figure('Name','E1 Stokes profile (PIV)','Position',[100 100 800 700],'Color','w');
semilogx(uS_z,    z*1e3, 'b-',  'LineWidth', 1.6); hold on;
semilogx(uS_mono, z*1e3, 'r--', 'LineWidth', 1.6);
xlabel('u^S(z)  (m/s)'); ylabel('z (mm)');
legend('Hilbert reconstruction','monochromatic at k_p','Location','best');
title(sprintf('E.1  Reconstructed Stokes drift profile  (k_p \\approx %.0f rad/m)', k_p));
grid on;
save_figure(fig, 'E1_stokes_profile');
end
