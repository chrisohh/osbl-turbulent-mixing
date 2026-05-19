function plot_C2_local_k_hovmoller(data)
% C.2  Local wavenumber Hovmöller k(x,t) from 1D Hilbert.

eta   = data.eta;
x_eta = data.x_eta;
t_eta = data.t_eta;

info = piv_eta_hilbert(eta, x_eta, t_eta(:));
k = double(info.k);

% Median peak k_p
k_p = median(k(:), 'omitnan');

fig = figure('Name','C2 local k Hovmoller','Position',[100 100 1100 600],'Color','w');
imagesc(t_eta, x_eta * 1e3, k);
set(gca,'YDir','normal');
cb = colorbar; cb.Label.String = 'k(x,t)  (rad/m)';
colormap(hot);
caxis([0, quantile(k(:), 0.99)]);
xlabel('t (s)'); ylabel('x (mm)');
title(sprintf('C.2  Local wavenumber k(x,t)  — median k_p \\approx %.0f rad/m', k_p));
save_figure(fig, 'C2_local_k_hovmoller');
end
