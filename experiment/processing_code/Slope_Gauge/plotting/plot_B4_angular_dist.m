function plot_B4_angular_dist(data)
% B.4  Angular distribution at omega_p; cos^(2s) fit; spreading parameter s.

eta = data.eta;
roi = data.setup.roi;
dx  = data.dx;  dy = data.dy;
Fs  = data.setup.frame_rate;

eta_roi = double(eta(roi.row_min:roi.row_max, roi.col_min:roi.col_max, :));
[Ny, Nx, Nt] = size(eta_roi);

% Find omega_p from spatially-averaged eta
eta_bar = squeeze(mean(mean(eta_roi, 1), 2));
[Pf, ff] = pwelch(eta_bar - mean(eta_bar), hann(min(256,Nt)), [], [], Fs);
[~, ip] = max(Pf);
fp = ff(ip);
omega_p = 2*pi * fp;
fprintf('B.4: omega_p = %.2f rad/s (f_p = %.2f Hz)\n', omega_p, fp);

% Bandpass eta_xy by omega_p (FFT in time, keep narrow band)
df = 0.10 * fp;                  % +/- 10% band
mask_t = abs(ff - fp) <= df;
% Fourier in time, zero outside band (full + mirror)
eta_t_fft = fft(eta_roi - mean(eta_roi(:)), [], 3);
freq = (0:Nt-1) * Fs / Nt;
freq(freq > Fs/2) = freq(freq > Fs/2) - Fs;
keep = (abs(freq) >= fp - df) & (abs(freq) <= fp + df);
eta_t_fft(:,:,~keep) = 0;
eta_band = real(ifft(eta_t_fft, [], 3));

% Time-averaged 2D power |eta_hat(kx,ky)|^2
P_acc = zeros(Ny, Nx);
win = hann(Ny) * hann(Nx)';
for it = 1:Nt
    e = eta_band(:,:,it); e = e - mean(e(:));
    eh = fft2(e .* win);
    P_acc = P_acc + abs(eh).^2;
end
P_acc = fftshift(P_acc / Nt);
kx = fftshift(2*pi * (-Nx/2:Nx/2-1) / (Nx*dx));
ky = fftshift(2*pi * (-Ny/2:Ny/2-1) / (Ny*dy));
[KX, KY] = meshgrid(kx, ky);
Kmag  = sqrt(KX.^2 + KY.^2);
Theta = atan2(KY, KX);

% Estimate k_p from still-water dispersion at omega_p
k_p_est = fzero(@(k) dispersion_curves(k,'still') - omega_p, [10 2e3]);
dk = 0.20 * k_p_est;
ring = (Kmag >= k_p_est - dk) & (Kmag <= k_p_est + dk);

% Bin theta
nth = 36;
th_edges = linspace(-pi, pi, nth+1);
th_ctr   = 0.5*(th_edges(1:end-1) + th_edges(2:end));
I_th = zeros(nth, 1);
for it = 1:nth
    msk = ring & (Theta >= th_edges(it)) & (Theta < th_edges(it+1));
    if any(msk(:)), I_th(it) = mean(P_acc(msk)); end
end

% cos^(2s) fit
% I(theta) = A * |cos((theta - theta0)/2)|^(2s)
model = @(p, t) p(1) * abs(cos((t - p(2))/2)).^(2*p(3));
p0 = [max(I_th), 0, 5];
opts = optimset('Display','off');
try
    p_fit = lsqcurvefit(model, p0, th_ctr(:), I_th, [0 -pi 0.1], [Inf pi 100], opts);
catch
    p_fit = p0;
end
s_param = p_fit(3);

fig = figure('Name','B4 angular distribution','Position',[100 100 800 750],'Color','w');
polarplot(th_ctr, I_th, 'bo-', 'LineWidth', 1.4); hold on;
th_fine = linspace(-pi, pi, 360);
polarplot(th_fine, model(p_fit, th_fine), 'r-', 'LineWidth', 1.6);
title(sprintf('B.4  I(\\theta) at \\omega_p = %.1f rad/s,  s = %.2f', omega_p, s_param));
legend('data','cos^{2s} fit','Location','best');
save_figure(fig, 'B4_angular_dist');
end
