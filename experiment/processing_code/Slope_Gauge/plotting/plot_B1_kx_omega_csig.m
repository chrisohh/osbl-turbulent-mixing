function plot_B1_kx_omega_csig(data)
% B.1  (k_x, omega) spectrum from y-mean eta with still-water dispersion overlay.
% TODO: pass U_E and U_S from PIV cross-validation file once available — Doppler
% curves are skipped this pass.

eta = data.eta;
roi = data.setup.roi;
dx  = data.dx;
Fs  = data.setup.frame_rate;

eta_roi = double(eta(roi.row_min:roi.row_max, roi.col_min:roi.col_max, :));
% y-mean
eta_xt = squeeze(mean(eta_roi, 1));   % [Nx x Nt]
[Nx, Nt] = size(eta_xt);

% Hann tapers
wx = hann(Nx); wt = hann(Nt).';
eta_xt = (eta_xt - mean(eta_xt(:))) .* (wx * wt);

% 2D FFT
P = fftshift(abs(fft2(eta_xt)).^2);

kx_full = 2*pi * (-Nx/2:Nx/2-1) / (Nx*dx);
om_full = 2*pi * (-Nt/2:Nt/2-1) / (Nt/Fs);
% positive omega only
ip_om = om_full > 0;
P_pos = P(:, ip_om);
om_pos = om_full(ip_om);

fig = figure('Name','B1 kx-omega (CSIG)','Position',[100 100 1000 700],'Color','w');
imagesc(kx_full, om_pos, log10(P_pos.' + eps));
set(gca,'YDir','normal'); axis xy;
colormap(parula); cb = colorbar; cb.Label.String = 'log_{10}|\eta|^2';
xlabel('k_x (rad/m)'); ylabel('\omega (rad/s)');
title('B.1  (k_x, \omega) spectrum, y-mean \eta');
hold on;

% Dispersion overlay (still water, both branches +/-)
k_pos = linspace(20, max(abs(kx_full)), 400);
om_disp = dispersion_curves(k_pos, 'still');
plot( k_pos,  om_disp, 'r-', 'LineWidth', 1.6);
plot(-k_pos,  om_disp, 'r-', 'LineWidth', 1.6);
legend('','still-water dispersion','Location','best');

save_figure(fig, 'B1_kx_omega_csig');
end
