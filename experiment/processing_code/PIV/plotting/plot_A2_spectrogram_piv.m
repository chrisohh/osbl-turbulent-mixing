function plot_A2_spectrogram_piv(data)
% A.2 (PIV part)  Frequency spectrogram of spatially-averaged eta.

eta   = data.eta;
t_eta = data.t_eta;
Fs    = data.Fs;

eta_bar = mean(eta, 1);
eta_bar = eta_bar - mean(eta_bar);
Nt = length(eta_bar);
win = min(256, floor(Nt/4));
nover = floor(win/2);

[Seta, f, ts] = spectrogram(eta_bar, hann(win), nover, [], Fs);
P = abs(Seta).^2;

fig = figure('Name','A2 spectrogram (PIV)','Position',[100 100 1100 500],'Color','w');
imagesc(ts + t_eta(1), f, log10(P + eps));
set(gca,'YDir','normal'); axis xy;
cb = colorbar; cb.Label.String = 'log_{10}|\eta|^2';
xlabel('t (s)'); ylabel('\omega/2\pi (Hz)');
% Peak ridge
[~, ip] = max(P, [], 1);
hold on; plot(ts + t_eta(1), f(ip), 'k:', 'LineWidth', 1.2);
title('A.2  PIV |\eta(\omega,t)|^2  (spatial mean) — dotted: \omega_p(t)');
save_figure(fig, 'A2_spectrogram_piv');
end
