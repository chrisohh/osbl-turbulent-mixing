function plot_A1_eta_snapshots_piv(data)
% A.1 (PIV companion)  eta(x, t_*) at 5 evenly spaced times.

eta   = data.eta;
x_eta = data.x_eta;
t_eta = data.t_eta;
[~, Nt] = size(eta);

idx = round(linspace(1, Nt, 5));
clim = max(abs(eta(:, idx)), [], 'all') * 1000;     % mm
eta_mean_x = mean(eta, 2) * 1000;

fig = figure('Name','A1 eta snapshots (PIV)','Position',[100 100 1600 380],'Color','w');
tiledlayout(1,5,'TileSpacing','compact','Padding','compact');
for k = 1:5
    nexttile;
    plot(x_eta * 1e3, eta(:, idx(k)) * 1000, 'b-', 'LineWidth', 1.2); hold on;
    plot(x_eta * 1e3, eta_mean_x, 'k:', 'LineWidth', 0.8);
    ylim([-clim, clim]);
    xlabel('x (mm)'); ylabel('\eta (mm)');
    title(sprintf('t = %.2f s', t_eta(idx(k))));
    grid on;
end
sgtitle('A.1  PIV-extracted surface elevation snapshots');
save_figure(fig, 'A1_eta_snapshots_piv');
end
