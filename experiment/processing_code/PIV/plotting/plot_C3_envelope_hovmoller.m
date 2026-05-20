function plot_C3_envelope_hovmoller(data)
% C.3  Wave envelope Hovmöller a(x,t) = |eta_analytic(x,t)|.

eta   = data.eta;
x_eta = data.x_eta;
t_eta = data.t_eta;

info = piv_eta_hilbert(eta, x_eta, t_eta(:));
a = double(info.a);

% Modulation length scale L_mod from envelope autocorrelation 1/e width
% along x for the mid time
[~, it_mid] = min(abs(t_eta - mean(t_eta)));
a_mid = a(:, it_mid) - mean(a(:, it_mid));
nfft = 2 * length(a_mid);
af = fft(a_mid, nfft);
ac = real(ifft(abs(af).^2));
ac = ac(1:length(a_mid));
ac = ac / max(ac);
dx = mean(diff(x_eta));
ix_e = find(ac < exp(-1), 1, 'first');
if isempty(ix_e), L_mod = NaN; else, L_mod = ix_e * dx * 1e3; end   % mm

fig = figure('Name','C3 envelope Hovmoller','Position',[100 100 1100 600],'Color','w');
imagesc(t_eta, x_eta * 1e3, a * 1000);
set(gca,'YDir','normal');
cb = colorbar; cb.Label.String = 'a(x,t)  (mm)';
colormap(hot);
xlabel('t (s)'); ylabel('x (mm)');
if isnan(L_mod)
    title('C.3  Envelope a(x,t)');
else
    title(sprintf('C.3  Envelope a(x,t)   L_{mod} \\approx %.1f mm  (1/e width at t_{mid})', L_mod));
end
save_figure(fig, 'C3_envelope_hovmoller');
end
