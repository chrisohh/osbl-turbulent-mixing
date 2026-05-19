function plot_C4_cross_stream_coh(data)
% C.4  Cross-stream coherence rho_eta(Delta y).

eta = data.eta;
roi = data.setup.roi;
dy  = data.dy;

eta_roi = double(eta(roi.row_min:roi.row_max, roi.col_min:roi.col_max, :));
[Ny, Nx, Nt] = size(eta_roi);

% Subsample x and t for speed
x_idx = round(linspace(1, Nx, min(20, Nx)));
t_idx = round(linspace(1, Nt, min(60, Nt)));

% FFT-based autocorrelation along y (dim 1)
nfft = 2 * Ny;
rho_acc = zeros(Ny, 1);
nL = 0;
for it = t_idx
    for ix = x_idx
        s = eta_roi(:, ix, it); s = s - mean(s);
        sf = fft(s, nfft);
        ac = real(ifft(abs(sf).^2));
        ac = ac(1:Ny) ./ ac(1);
        rho_acc = rho_acc + ac;
        nL = nL + 1;
    end
end
rho = rho_acc / nL;
dy_axis = (0:Ny-1) * dy * 100;       % cm

fig = figure('Name','C4 cross-stream coherence','Position',[100 100 900 500],'Color','w');
plot(dy_axis, rho, 'b-', 'LineWidth', 1.6); hold on;
yline(0.9, 'k--', '\rho = 0.9');
yline(0.5, 'k:',  '\rho = 0.5');
xlabel('\Delta y (cm)'); ylabel('\rho_\eta(\Delta y)');
title('C.4  Cross-stream coherence of \eta');
ylim([min(0, min(rho)*1.1), 1.05]);
grid on;
save_figure(fig, 'C4_cross_stream_coh');
end
