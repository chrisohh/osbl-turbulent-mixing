function plot_A6_directional_B(data)
% A.6  Directional saturation B(k, theta) on a polar plot, late half of record.

eta = data.eta;
roi = data.setup.roi;
dx  = data.dx;  dy = data.dy;

eta_roi = double(eta(roi.row_min:roi.row_max, roi.col_min:roi.col_max, :));
[Ny, Nx, Nt] = size(eta_roi);

% Late half of record
it_start = round(Nt/2);
nF = Nt - it_start + 1;

% 2D FFT averaged over late frames
kx = 2*pi * ifftshift(((-Nx/2):(Nx/2-1)) / (Nx*dx));
ky = 2*pi * ifftshift(((-Ny/2):(Ny/2-1)) / (Ny*dy));
[KX, KY] = meshgrid(kx, ky);

P_acc = zeros(Ny, Nx);
win = hann(Ny) * hann(Nx)';
for it = it_start:Nt
    e = eta_roi(:,:,it); e = e - mean(e(:));
    eh = fft2(e .* win);
    P_acc = P_acc + abs(eh).^2;
end
P_acc = P_acc / nF;
P_acc = fftshift(P_acc);
KX_s = fftshift(KX); KY_s = fftshift(KY);

Kmag  = sqrt(KX_s.^2 + KY_s.^2);
Theta = atan2(KY_s, KX_s);                  % -pi..pi

% Polar binning
k_min = max(20, min(Kmag(Kmag > 0)));
k_max = max(Kmag(:));
nk = 30; nth = 36;
k_edges  = logspace(log10(k_min), log10(k_max), nk+1);
th_edges = linspace(-pi, pi, nth+1);
k_ctr  = sqrt(k_edges(1:end-1) .* k_edges(2:end));
th_ctr = 0.5*(th_edges(1:end-1) + th_edges(2:end));

E_kth = zeros(nk, nth);
for ik = 1:nk
    for it = 1:nth
        msk = (Kmag >= k_edges(ik)) & (Kmag < k_edges(ik+1)) & ...
              (Theta >= th_edges(it)) & (Theta < th_edges(it+1));
        if any(msk(:))
            dk = k_edges(ik+1) - k_edges(ik);
            dth = th_edges(it+1) - th_edges(it);
            E_kth(ik, it) = mean(P_acc(msk)) / max(dk*dth, eps);
        end
    end
end
B_kth = (k_ctr.').^3 .* E_kth;

% Plot
[TH, K] = meshgrid(th_ctr, k_ctr);
[Xp, Yp] = pol2cart(TH, K);

fig = figure('Name','A6 directional B(k,theta)','Position',[100 100 800 750],'Color','w');
pcolor(Xp, Yp, log10(B_kth + eps)); shading interp;
axis equal;
cb = colorbar; cb.Label.String = 'log_{10} B(k,\theta)';
% Wind direction arrow (assumed +x)
hold on;
ax_lim = max(k_ctr) * 0.6;
quiver(0, 0, ax_lim, 0, 'k', 'LineWidth', 2, 'MaxHeadSize', 0.3);
text(ax_lim*1.05, 0, 'wind', 'FontWeight','bold');
xlabel('k_x (rad/m)'); ylabel('k_y (rad/m)');
title('A.6  Directional saturation B(k,\theta), late half of record');
colormap(hot);
save_figure(fig, 'A6_directional_B');
end
