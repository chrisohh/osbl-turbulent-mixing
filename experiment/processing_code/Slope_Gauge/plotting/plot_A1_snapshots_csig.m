function plot_A1_snapshots_csig(data)
% A.1  Five eta(x,y,t_*) snapshots at evenly spaced times.
%%
eta   = data.eta;
time  = data.time;
setup = data.setup;
dx    = data.dx;  dy = data.dy;

[Ny, Nx, Nt] = size(eta);
roi = setup.roi;
x_cm = ((1:Nx) - Nx/2) * dx * 100;
y_cm = ((1:Ny) - Ny/2) * dy * 100;

idx = round(linspace(1, Nt, 5));

% Common color scale: 99th percentile across all 5 frames (in ROI).
abs_pool = [];
for k = 1:5
    f = abs(double(eta(roi.row_min:roi.row_max, roi.col_min:roi.col_max, idx(k))));
    abs_pool = [abs_pool; f(:)]; %#ok<AGROW>
end
clim = quantile(abs_pool, 0.99);
if clim == 0, clim = max(abs_pool) + eps; end

fig = figure('Name','A1 eta snapshots (CSIG)','Position',[100 100 1800 380],'Color','w');
tiledlayout(1, 5, 'TileSpacing','compact','Padding','compact');
for k = 1:5
    nexttile;
    f = double(eta(:,:,idx(k))) * 1000;            % mm
    imagesc(x_cm, y_cm, f);
    axis equal tight;
    set(gca,'YDir','normal');
    caxis([-clim*1000, clim*1000]);
    colormap(gca, bwr_cmap());
    title(sprintf('t = %.2f s', time(idx(k))));
    xlabel('x (cm)'); ylabel('y (cm)');
end
cb = colorbar; cb.Label.String = '\eta (mm)';
sgtitle('A.1  CSIG surface elevation snapshots');
% save_figure(fig, 'A1_snapshots_csig');
end


function cmap = bwr_cmap()
n = 256; h = floor(n/2);
r = [linspace(0,1,h), ones(1,n-h)];
g = [linspace(0,1,h), linspace(1,0,n-h)];
b = [ones(1,h),       linspace(1,0,n-h)];
cmap = [r' g' b'];
end
