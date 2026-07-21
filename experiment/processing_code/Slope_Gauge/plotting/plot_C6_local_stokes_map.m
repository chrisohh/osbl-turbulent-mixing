function plot_C6_local_stokes_map(data)
% C.6  Local Stokes drift map u_S(x, y, z=0, t) = a^2 * omega * k_x.

eta = data.eta;
roi = data.setup.roi;
dx  = data.dx;  dy = data.dy;
time = data.time;
Fs   = data.setup.frame_rate;

eta_roi = single(eta(roi.row_min:roi.row_max, roi.col_min:roi.col_max, :));
[Ny, Nx, Nt] = size(eta_roi);
x_cm = ((1:Nx) - Nx/2) * dx * 100;
y_cm = ((1:Ny) - Ny/2) * dy * 100;

idx_snap = round(linspace(1, Nt, 4));

% Snapshot maps
[a_s, kx_s, ky_s, ~] = cisg_local_wavenumber(eta_roi(:,:,idx_snap), dx, dy);
kmag_s = sqrt(double(kx_s).^2 + double(ky_s).^2);
om_s   = sqrt(9.81 .* kmag_s + 7.4e-5 .* kmag_s.^3);
uS_s   = double(a_s).^2 .* om_s .* double(kx_s);

% Time series of mean and std (chunked)
chunk = 50;
mean_uS = nan(Nt, 1);
std_uS  = nan(Nt, 1);
for c0 = 1:chunk:Nt
    c1 = min(c0 + chunk - 1, Nt);
    [a, kx, ky, ~] = cisg_local_wavenumber(eta_roi(:,:,c0:c1), dx, dy);
    kmag = sqrt(double(kx).^2 + double(ky).^2);
    om   = sqrt(9.81 .* kmag + 7.4e-5 .* kmag.^3);
    uS   = double(a).^2 .* om .* double(kx);
    for j = 1:(c1 - c0 + 1)
        v = uS(:,:,j); v = v(isfinite(v));
        if ~isempty(v)
            mean_uS(c0+j-1) = mean(v);
            std_uS(c0+j-1)  = std(v);
        end
    end
end
nroll = max(1, round(Fs));
mean_s = movmean(mean_uS, nroll, 'omitnan');
std_s  = movmean(std_uS,  nroll, 'omitnan');
ratio  = std_s ./ max(abs(mean_s), eps);

fig = figure('Name','C6 local Stokes map','Position',[100 100 1700 900],'Color','w');
tiledlayout(3, 4, 'TileSpacing','compact','Padding','compact');
clim_uS = quantile(abs(uS_s(isfinite(uS_s))), 0.99);
for k = 1:4
    nexttile;
    imagesc(x_cm, y_cm, uS_s(:,:,k));
    axis equal tight; set(gca,'YDir','normal');
    caxis([-clim_uS, clim_uS]);
    colormap(gca, bwr_cmap()); colorbar;
    title(sprintf('u^S(x,y,0)  (t = %.2f s)', time(idx_snap(k))));
    xlabel('x (cm)'); ylabel('y (cm)');
end

% Bottom: time series row spanning full width
nexttile([1 4]);
yyaxis left
plot(time, mean_s, 'b-', 'LineWidth', 1.4); hold on;
plot(time, std_s,  'r-', 'LineWidth', 1.4);
ylabel('u^S  (m/s)');
yyaxis right
plot(time, ratio, 'k-', 'LineWidth', 1.4);
ylabel('\sigma_{u^S} / |<u^S>|');
xlabel('t (s)');
legend('<u^S>','\sigma_{u^S}','ratio','Location','best');
title('Time series and heterogeneity ratio');
grid on;

% Empty next row to make tiling 3x4 work
nexttile([1 4]); axis off;

sgtitle('C.6  Local Stokes drift heterogeneity');
save_figure(fig, 'C6_local_stokes_map');
end


function cmap = bwr_cmap()
n = 256; h = floor(n/2);
r = [linspace(0,1,h), ones(1,n-h)];
g = [linspace(0,1,h), linspace(1,0,n-h)];
b = [ones(1,h),       linspace(1,0,n-h)];
cmap = [r' g' b'];
end
