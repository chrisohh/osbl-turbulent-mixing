function plot_C5_heterogeneity_idx(data)
% C.5  Heterogeneity indices H_a(t) = sigma_a / <a>, H_k(t) = sigma_kx / <kx>.

eta = data.eta;
roi = data.setup.roi;
dx  = data.dx;  dy = data.dy;
time = data.time;
Fs   = data.setup.frame_rate;

eta_roi = single(eta(roi.row_min:roi.row_max, roi.col_min:roi.col_max, :));
[~, ~, Nt] = size(eta_roi);

% Process in chunks to keep memory modest
chunk = 50;
H_a = nan(Nt, 1);  H_k = nan(Nt, 1);
for c0 = 1:chunk:Nt
    c1 = min(c0 + chunk - 1, Nt);
    [a, kx_loc, ~, ~] = cisg_local_wavenumber(eta_roi(:,:,c0:c1), dx, dy);
    for j = 1:(c1 - c0 + 1)
        a_j  = double(a(:,:,j));   a_j  = a_j(isfinite(a_j) & a_j > 0);
        kx_j = double(kx_loc(:,:,j)); kx_j = kx_j(isfinite(kx_j) & kx_j > 0);
        if ~isempty(a_j),  H_a(c0+j-1) = std(a_j)  / mean(a_j);  end
        if ~isempty(kx_j), H_k(c0+j-1) = std(kx_j) / mean(kx_j); end
    end
end

% 0.5-s rolling mean
nroll = max(1, round(0.5 * Fs));
H_a_s = movmean(H_a, nroll, 'omitnan');
H_k_s = movmean(H_k, nroll, 'omitnan');

fig = figure('Name','C5 heterogeneity indices','Position',[100 100 1000 500],'Color','w');
plot(time, H_a_s, 'b-', 'LineWidth', 1.6); hold on;
plot(time, H_k_s, 'r-', 'LineWidth', 1.6);
yline(0, 'k:');
xlabel('t (s)'); ylabel('heterogeneity index');
legend('H_a = \sigma_a / <a>', 'H_k = \sigma_{kx} / <k_x>', 'Location','best');
title('C.5  Heterogeneity indices  (0.5-s rolling mean)');
grid on;
save_figure(fig, 'C5_heterogeneity_idx');
end
