function plot_B1_kx_omega_csig(data)
% B.1  (k_x, omega) spectrum via 3-D FFT in (x, y, t), integrated over k_y.
%
% Pipeline:
%   eta(x, y, t)  -->  3-D Hann window
%                 -->  fftn  ->  E(kx, ky, omega)
%                 -->  sum |E|^2 over ky  ->  S(kx, omega)
%
% Compared with the y-mean / 2-D FFT approach, this captures energy from
% ALL cross-wind wavenumbers (oblique / short-crested waves), not only ky=0.
%
% Memory note: requires ~2x the RAM of the eta cube in complex double.
% For a 500x1000x200 cube that is ~3 GB — use a frame subset if needed.

eta = data.eta;
roi = data.setup.roi;
dx  = data.dx;
dy  = data.dy;
Fs  = data.setup.frame_rate;

eta_roi = double(eta(roi.row_min:roi.row_max, roi.col_min:roi.col_max, :));
[Ny, Nx, Nt] = size(eta_roi);

fprintf('B1 3-D FFT: eta_roi is [%d x %d x %d]  (%.1f MB double)\n', ...
        Ny, Nx, Nt, Ny*Nx*Nt*8/1e6);

%% 3-D Hann taper in all dimensions
wy   = hann(Ny);                        % [Ny x 1]
wx   = hann(Nx).';                      % [1  x Nx]
wt   = reshape(hann(Nt), 1, 1, Nt);    % [1  x 1  x Nt]
win3 = wy .* wx .* wt;                  % broadcasts to [Ny x Nx x Nt]

eta_win = (eta_roi - mean(eta_roi(:))) .* win3;

%% 3-D FFT  ->  shift DC to centre  ->  power
F3 = fftshift(fftn(eta_win));           % [Ny x Nx x Nt], complex
P3 = abs(F3).^2;                        % [Ny x Nx x Nt], real

%% Wavenumber / frequency axes (two-sided, centred)
kx_full = 2*pi * (-Nx/2 : Nx/2-1) / (Nx * dx);   % rad/m
ky_full = 2*pi * (-Ny/2 : Ny/2-1) / (Ny * dy);   % rad/m  (for reference)
om_full = 2*pi * (-Nt/2 : Nt/2-1) / (Nt / Fs);   % rad/s

dky = ky_full(2) - ky_full(1);                    % rad/m

%% Integrate (sum x dky) over k_y  ->  S(kx, omega)  [Nx x Nt]
P_kx_om = squeeze(sum(P3, 1)) * dky;              % dim 1 = ky after fftshift

%% Keep positive omega only
ip_om = om_full > 0;
P_pos  = P_kx_om(:, ip_om);    % [Nx x n_om]
om_pos = om_full(ip_om);        % [1  x n_om]

%% Plot
fig = figure('Name','B1 kx-omega 3D FFT (CSIG)', ...
             'Position',[100 100 1000 700],'Color','w');

imagesc(kx_full, om_pos, log10(P_pos.' + eps));
set(gca, 'YDir', 'normal');
colormap(parula);
cb = colorbar;
cb.Label.String     = '$\log_{10} \int |\hat{\eta}(k_x,k_y,\omega)|^2 \, dk_y$';
cb.Label.Interpreter = 'latex';
cb.Label.FontSize    = 16;

xlabel('$k_x$ (rad/m)',  'Interpreter', 'latex');
ylabel('$\omega$ (rad/s)', 'Interpreter', 'latex');
title('B.1  $(k_x, \omega)$ spectrum, $k_y$-integrated 3-D FFT', ...
      'Interpreter', 'latex');
set(gca, 'fontsize', 14, 'fontname', 'times');
hold on;

%% Dispersion overlay (still water)
k_pos    = linspace(20, max(abs(kx_full)), 400);
om_disp  = dispersion_curves(k_pos, 'still');
plot( k_pos, om_disp, 'r-', 'LineWidth', 1.6);
plot(-k_pos, om_disp, 'r-', 'LineWidth', 1.6);
legend('', 'still-water dispersion', 'Location', 'best');

% save_figure(fig, 'B1_kx_omega_csig');
end
