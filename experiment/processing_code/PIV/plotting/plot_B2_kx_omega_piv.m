function plot_B2_kx_omega_piv(data)
% B.2  (k_x, omega) spectrum from PIV eta(x,t) with still-water dispersion overlay.
% Doppler curves skipped this pass.

eta   = data.eta;
x_eta = data.x_eta;
Fs    = data.Fs;
dx    = mean(diff(x_eta));

[Nx, Nt] = size(eta);
wx = hann(Nx); wt = hann(Nt).';
e = (eta - mean(eta(:))) .* (wx * wt);

P = fftshift(abs(fft2(e)).^2);
kx_full = 2*pi * (-Nx/2:Nx/2-1) / (Nx*dx);
om_full = 2*pi * (-Nt/2:Nt/2-1) / (Nt/Fs);

ip_om = om_full > 0;
P_pos = P(:, ip_om);
om_pos = om_full(ip_om);

fig = figure('Name','B2 kx-omega (PIV)','Position',[100 100 1000 700],'Color','w');
imagesc(kx_full, om_pos, log10(P_pos.' + eps));
set(gca,'YDir','normal'); axis xy;
colormap(parula); cb = colorbar; cb.Label.String = 'log_{10}|\eta|^2';
xlabel('k_x (rad/m)'); ylabel('\omega (rad/s)');
title('B.2  PIV (k_x, \omega) spectrum');
hold on;

k_pos = linspace(20, max(abs(kx_full)), 400);
om_disp = dispersion_curves(k_pos, 'still');
plot( k_pos,  om_disp, 'r-', 'LineWidth', 1.6);
plot(-k_pos,  om_disp, 'r-', 'LineWidth', 1.6);
legend('','still-water dispersion','Location','best');
save_figure(fig, 'B2_kx_omega_piv');
end
