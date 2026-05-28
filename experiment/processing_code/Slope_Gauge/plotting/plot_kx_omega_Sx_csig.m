function plot_kx_omega_Sx_csig(data)
% (k_x, omega) spectrum of along-wind slope Sx with linear gravity-capillary
% dispersion overlay.  y-mean Sx is used before the 2D FFT.

Sx  = data.Sx;
dx  = data.dx;
Fs  = data.setup.frame_rate;

Sx_xt  = squeeze(mean(Sx, 1));      % [Nx x Nt], y-mean
[Nx, Nt] = size(Sx_xt);

% Detrend and apply 2D Hann taper
wx     = hann(Nx); wt = hann(Nt).';
Sx_xt  = (Sx_xt - mean(Sx_xt(:))) .* (wx * wt);

% 2D FFT
P       = fftshift(abs(fft2(Sx_xt)).^2);
kx_full = 2*pi * (-Nx/2 : Nx/2-1) / (Nx * dx);
om_full = 2*pi * (-Nt/2 : Nt/2-1) / (Nt / Fs);

% Positive-omega half
ip     = om_full > 0;
P_pos  = P(:, ip);
om_pos = om_full(ip);

fig = figure('Name', 'kx-omega Sx (CSIG)', 'Position', [100 100 900 700], 'Color', 'w');
imagesc(kx_full, om_pos, log10(P_pos.' + eps));
set(gca, 'YDir', 'normal'); axis xy; axis tight;
colormap(bone);
cb = colorbar;
cb.Label.String      = '$\log_{10}|S_x(k_x,\omega)|^2$';
cb.Label.Interpreter = 'latex';
cb.Label.FontSize    = 14;
xlabel('$k_x$ (rad/m)',    'Interpreter', 'latex');
ylabel('$\omega$ (rad/s)', 'Interpreter', 'latex');
title('$(k_x,\,\omega)$ spectrum of $S_x$', 'Interpreter', 'latex');
set(gca, 'fontsize', 14, 'fontname', 'times');
hold on;

% Linear gravity-capillary dispersion: omega = sqrt(g*|k| + (sigma/rho)*|k|^3)
g     = 9.81;
gamma = 7.4e-5;        % sigma/rho for clean water at ~20 C  (m^3 s^-2)
k_vec  = linspace(1, max(abs(kx_full)), 600);
om_gc  = sqrt(g .* k_vec + gamma .* k_vec.^3);

plot( k_vec, om_gc, 'r-', 'LineWidth', 1.8);
plot(-k_vec, om_gc, 'r-', 'LineWidth', 1.8);
legend('', 'gravity-capillary', 'Interpreter', 'none', 'Location', 'best');

% save_figure(fig, 'kx_omega_Sx_csig');
end
