function plot_kx_omega_Sx_csig(data, method, kx_lp)
% (k_x, f) spectrum of along-wind slope Sx using the entire time series.
%
% Usage:
%   plot_kx_omega_Sx_csig(data)
%   plot_kx_omega_Sx_csig(data, '2d')        % ky=0 slice (Leckler et al. 2015, default)
%   plot_kx_omega_Sx_csig(data, '3d')        % 3-D FFT, ky-integrated
%   plot_kx_omega_Sx_csig(data, '2d', kx_lp) % with Gaussian LP at kx_lp (rad/m)
%
% Suggested kx_lp cutoffs (rad/m):
%   k_pixel_Nyq = pi/dx        (absolute pixel Nyquist)
%   k_pattern   = 2*pi/p_LED   (LED pitch p_LED in metres — you supply)
%   k_snr       ~ 300-500      (empirical: where noise plateau begins)
%   k_temporal  = (pi*Fs)^2/g  (deep-water k at temporal Nyquist)

if nargin < 2 || isempty(method),  method = '2d'; end
if nargin < 3,                     kx_lp  = [];   end
method = lower(method);

Sx   = data.Sx;
dx   = data.dx;
dy   = data.dy;
Fs   = data.setup.frame_rate;
time = data.time;

[Ny, Nx, Nt] = size(Sx);

%% Print physical cutoffs
k_pixel_Nyq = pi / dx;
k_temporal  = (pi * Fs)^2 / 9.81;
T_total     = Nt / Fs;
fprintf('--- kx cutoffs (rad/m) ---\n');
fprintf('  Pixel Nyquist       : %6.0f  rad/m  (lambda = %.2f mm)\n', ...
        k_pixel_Nyq, 2*pi/k_pixel_Nyq*1e3);
fprintf('  Temporal Nyquist k  : %6.0f  rad/m  (deep-water at f=%.1f Hz)\n', ...
        k_temporal, Fs/2);
fprintf('  Empirical SNR ~      300-500 rad/m  (read from spectrum)\n');
fprintf('  LED pattern pitch    : supply p_LED (m) -> k = 2*pi/p_LED\n');
if ~isempty(kx_lp)
    fprintf('  LP filter applied at : %6.0f  rad/m  (lambda = %.2f mm)\n', ...
            kx_lp, 2*pi/kx_lp*1e3);
end
fprintf('  Full record          : %.1f s  (%d frames, df = %.3f Hz)\n', ...
        T_total, Nt, 1/T_total);
fprintf('--------------------------\n');

%% Optional Gaussian LP filter in x
if ~isempty(kx_lp)
    sigma_m  = 1 / kx_lp;
    sigma_px = sigma_m / dx;
    h_half   = ceil(4 * sigma_px);
    xk       = -h_half : h_half;
    kernel   = exp(-0.5 * (xk / sigma_px).^2);
    kernel   = kernel / sum(kernel);
    Sx_f = zeros(size(Sx), 'single');
    for t = 1:Nt
        for r = 1:Ny
            Sx_f(r,:,t) = conv(double(Sx(r,:,t)), kernel, 'same');
        end
    end
    Sx = Sx_f;  clear Sx_f;
end

%% Axes
kx_full = 2*pi * (-Nx/2 : Nx/2-1) / (Nx * dx);   % rad/m  [Nx x 1]

%% Compute spectrum over entire time series
if strcmp(method, '2d')
    Sx_xt = squeeze(mean(double(Sx), 1));         % [Nx x Nt]
    wx    = hann(Nx); wt = hann(Nt).';
    Sx_xt = (Sx_xt - mean(Sx_xt(:))) .* (wx * wt);
    P_raw = fftshift(abs(fft2(Sx_xt)).^2);        % [Nx x Nt]

elseif strcmp(method, '3d')
    wy   = hann(Ny);
    wx   = hann(Nx).';
    wt3  = reshape(hann(Nt), 1, 1, Nt);
    win3 = wy .* wx .* wt3;
    Sx_d = double(Sx) - mean(Sx(:), 'omitnan');
    F3   = fftshift(fftn(Sx_d .* win3));
    P3   = abs(F3).^2;
    ky_v = 2*pi * (-Ny/2 : Ny/2-1) / (Ny * dy);
    P_raw = squeeze(sum(P3, 1)) * (ky_v(2)-ky_v(1));

else
    error('method must be ''2d'' or ''3d''.');
end

f_full = (-Nt/2 : Nt/2-1) / (Nt / Fs);   % Hz
ip     = f_full > 0;
P_pos  = P_raw(:, ip);                     % [Nx x n_f]
f_pos  = f_full(ip);

%% Dispersion curve
g     = 9.81;  gamma = 7.4e-5;
k_vec = linspace(1, max(abs(kx_full)), 600);
f_gc  = sqrt(g.*k_vec + gamma.*k_vec.^3) / (2*pi);

lp_str = '';
if ~isempty(kx_lp), lp_str = sprintf(', LP%d', round(kx_lp)); end

%% Plot
fig = figure('Name', sprintf('kx-f Sx [%s%s] (CSIG)', method, lp_str), ...
             'Position', [100 80 900 700], 'Color', 'w');

imagesc(kx_full, f_pos, log10(P_pos.' + eps));
set(gca, 'YDir', 'normal'); axis xy; axis tight;
colormap(inferno);
cb = colorbar;
cb.Label.String      = '$\log_{10}|S_x(k_x,f)|^2$';
cb.Label.Interpreter = 'latex';
cb.Label.FontSize     = 14;
xlabel('$k_x$ (rad m$^{-1}$)', 'Interpreter', 'latex');
ylabel('$f$ (Hz)',              'Interpreter', 'latex');
title(sprintf('$(k_x, f)$ spectrum, $S_x$ [%s%s]', method, lp_str), ...
      'Interpreter', 'latex');
set(gca, 'fontsize', 14, 'fontname', 'times');
hold on;
plot( k_vec, f_gc, 'r-', 'LineWidth', 1.5);
plot(-k_vec, f_gc, 'r-', 'LineWidth', 1.5);
ylim([0, max(f_pos)]);
xlim([0, 2000]);

% save_figure(fig, 'kx_f_Sx_csig');
end
