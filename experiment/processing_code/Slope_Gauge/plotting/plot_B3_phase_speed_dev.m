function plot_B3_phase_speed_dev(data)
% B.3  Phase-speed deviation Delta c(k) from ridge tracking on B.1 image.

eta = data.eta;
roi = data.setup.roi;
dx  = data.dx;
Fs  = data.setup.frame_rate;

eta_roi = double(eta(roi.row_min:roi.row_max, roi.col_min:roi.col_max, :));
eta_xt = squeeze(mean(eta_roi, 1));
[Nx, Nt] = size(eta_xt);
wx = hann(Nx); wt = hann(Nt).';
eta_xt = (eta_xt - mean(eta_xt(:))) .* (wx * wt);

P = fftshift(abs(fft2(eta_xt)).^2);
kx_full = 2*pi * (-Nx/2:Nx/2-1) / (Nx*dx);
om_full = 2*pi * (-Nt/2:Nt/2-1) / (Nt/Fs);
ip_om = om_full > 0;
P_pos  = P(:, ip_om);
om_pos = om_full(ip_om);

% Positive-k only (right-going waves)
ip_k = kx_full > 0;
kx_pos = kx_full(ip_k);
P_q = P_pos(ip_k, :);

om_ridge = nan(length(kx_pos), 1);
om_lo    = nan(length(kx_pos), 1);
om_hi    = nan(length(kx_pos), 1);
for ik = 1:length(kx_pos)
    slice = P_q(ik, :);
    if max(slice) <= 0, continue; end
    [pmax, im] = max(slice);
    om_ridge(ik) = om_pos(im);
    half = pmax / 2;
    above = slice >= half;
    if any(above)
        idxs = find(above);
        om_lo(ik) = om_pos(idxs(1));
        om_hi(ik) = om_pos(idxs(end));
    end
end

c_ridge = om_ridge ./ kx_pos(:);
c0      = dispersion_curves(kx_pos(:), 'still') ./ kx_pos(:);
dc      = c_ridge - c0;
dc_lo   = om_lo ./ kx_pos(:) - c0;
dc_hi   = om_hi ./ kx_pos(:) - c0;

fig = figure('Name','B3 phase-speed deviation','Position',[100 100 900 500],'Color','w');
fill([kx_pos, fliplr(kx_pos)], [dc_lo.', fliplr(dc_hi.')], [0.8 0.8 1], ...
     'EdgeColor','none','FaceAlpha',0.6); hold on;
plot(kx_pos, dc, 'b-', 'LineWidth', 1.6);
yline(0, 'k--');
xlabel('k_x (rad/m)'); ylabel('\Delta c = c_{ridge} - c_0 (m/s)');
title('B.3  Phase-speed deviation (full record; FWHM band)');
grid on;
save_figure(fig, 'B3_phase_speed_dev');
end
