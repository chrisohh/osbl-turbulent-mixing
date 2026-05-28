function plot_A1_snapshots_csig(data, frame_indices)
% A.1  eta(x,y,t_*) snapshots.
% Optional frame_indices: vector of TIME INDICES into data.eta (not raw frame
% numbers). Defaults to 5 evenly spaced indices.
%%
eta   = data.eta;
Sx = data.Sx;Sy = data.Sy;
time  = data.time;
setup = data.setup;
dx    = data.dx;  dy = data.dy;

[Ny, Nx, Nt] = size(eta);
x_cm = ((1:Nx) - Nx/2) * dx * 100;
y_cm = ((1:Ny) - Ny/2) * dy * 100;

if nargin < 2 || isempty(frame_indices)
    idx = round(linspace(1, Nt, 5));
else
    idx = frame_indices(:).';
end
n_plots = numel(idx);

% Common color scale: 99th percentile across selected frames (in ROI).
abs_pool = [];
for k = 1:n_plots
    f = abs(double(eta(:,:, idx(k))));
    abs_pool = [abs_pool; f(:)]; %#ok<AGROW>
end
clim = quantile(abs_pool, 0.99);
if clim == 0, clim = max(abs_pool) + eps; end

fig = figure('Name','A1 eta snapshots (CSIG)','Position',[100 100 1800 380],'Color','w');
tiledlayout(1, n_plots, 'TileSpacing','compact','Padding','compact');
for k = 1:n_plots
    nexttile;
    f = double(eta(:,:,idx(k))) * 100;            % cm
    imagesc(x_cm, y_cm, f);
    axis equal tight;
    set(gca,'YDir','normal');
    caxis([-clim*100, clim*100]);
    colormap(gca, brewermap([], 'Spectral'));
    title(sprintf('$t = %.2f$ s', time(idx(k))), 'Interpreter', 'latex');
  xlabel('$x$ (cm)', 'interpreter', 'latex');
ylabel('$y$ (cm)', 'interpreter', 'latex');
    set(gca, 'fontsize', 14, 'fontname', 'times');
end
cb1 = colorbar; cb1.Label.String = '$\eta$ (cm)';cb1.Label.Interpreter = 'latex'; cb1.Label.FontSize = 16;

fig = figure('Name','Sx Sy snapshots (CSIG)','Position',[100 100 1800 800],'Color','w');
tiledlayout(2, n_plots, 'TileSpacing','compact','Padding','compact');
for k = 1:n_plots
    nexttile;
    f = double(Sx(:,:,idx(k)));          
    imagesc(x_cm, y_cm, f);
    axis equal tight;
    set(gca,'YDir','normal');
    caxis([-0.4, 0.4]);
    colormap(gca, brewermap([], 'Spectral'));
    title(sprintf('$t = %.2f$ s', time(idx(k))), 'Interpreter', 'latex');
  xlabel('$x$ (cm)', 'interpreter', 'latex');
ylabel('$y$ (cm)', 'interpreter', 'latex');
    set(gca, 'fontsize', 14, 'fontname', 'times');
end
cb2 = colorbar; cb2.Label.String = '$S_x$'; cb2.Label.Interpreter = 'latex'; cb2.Label.FontSize = 16;
for k = 1:n_plots
    nexttile;
    f = double(Sy(:,:,idx(k)));     
    imagesc(x_cm, y_cm, f);
    axis equal tight;
    set(gca,'YDir','normal');
    caxis([-0.2, 0.2]);
    colormap(gca, brewermap([], 'Spectral'));
xlabel('$x$ (cm)', 'interpreter', 'latex');
ylabel('$y$ (cm)', 'interpreter', 'latex');
    set(gca, 'fontsize', 14, 'fontname', 'times');
end
cb3 = colorbar; cb3.Label.String = '$S_y$'; cb3.Label.Interpreter = 'latex'; cb3.Label.FontSize = 16;
% save_figure(fig, 'A1_snapshots_csig');
end
