function plot_C1_local_k_map(data)
% C.1  Local wavenumber map k_x(x,y,t) at 4 snapshots + histogram.

eta = data.eta;
roi = data.setup.roi;
dx  = data.dx;  dy = data.dy;
time = data.time;

eta_roi = double(eta(roi.row_min:roi.row_max, roi.col_min:roi.col_max, :));
[Ny, Nx, Nt] = size(eta_roi);
x_cm = ((1:Nx) - Nx/2) * dx * 100;
y_cm = ((1:Ny) - Ny/2) * dy * 100;

idx = round(linspace(1, Nt, 4));

% Compute local wavenumber for snapshots only
sel = single(eta_roi(:,:,idx));
[~, kx_sel, ~, ~] = cisg_local_wavenumber(sel, dx, dy);

% Compute kx histogram from a sparser set of snapshots across record
hist_idx = round(linspace(1, Nt, min(40, Nt)));
hist_sel = single(eta_roi(:,:,hist_idx));
[~, kx_h, ~, ~] = cisg_local_wavenumber(hist_sel, dx, dy);
kx_pool = double(kx_h(:));
kx_pool = kx_pool(isfinite(kx_pool) & kx_pool > 0 & kx_pool < 2000);

fig = figure('Name','C1 local k map','Position',[100 100 1700 800],'Color','w');
tiledlayout(2, 4, 'TileSpacing','compact','Padding','compact');
clim = quantile(kx_pool, 0.99);
for k = 1:4
    nexttile;
    imagesc(x_cm, y_cm, double(kx_sel(:,:,k)));
    axis equal tight; set(gca,'YDir','normal');
    caxis([0, clim]); colormap(gca, hot); colorbar;
    title(sprintf('k_x  (t = %.2f s)', time(idx(k))));
    xlabel('x (cm)'); ylabel('y (cm)');
end

% Histogram across bottom row
nexttile([1 4]);
histogram(kx_pool, 80, 'Normalization','pdf','FaceColor',[0.2 0.4 0.8]);
xlabel('k_x (rad/m)'); ylabel('pdf');
title('Histogram of local k_x (sampled across record)');
grid on;

sgtitle('C.1  Local along-x wavenumber from 2D Riesz transform');
save_figure(fig, 'C1_local_k_map');
end
