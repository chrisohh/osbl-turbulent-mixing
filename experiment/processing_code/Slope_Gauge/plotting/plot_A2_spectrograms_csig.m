function plot_A2_spectrograms_csig(data)
% A.2  Frequency spectrograms of S_x, S_y, and eta(t) (spatially averaged).

Sx    = data.Sx;
Sy    = data.Sy;
eta   = data.eta;
time  = data.time(:);
setup = data.setup;
roi   = setup.roi;
Fs    = setup.frame_rate;

[~, ~, Nt] = size(Sx);
win = min(256, floor(Nt/4));
nover = floor(win/2);

% --- spatially-averaged eta time series ---
eta_bar = squeeze(mean(mean(eta(roi.row_min:roi.row_max, ...
                              roi.col_min:roi.col_max, :), 1), 2));
eta_bar = double(eta_bar(:));

% --- per-row Welch averaged over y for Sx, Sy ---
[Sx_pxx, f_psd] = avg_welch_y(Sx, roi, win, nover, Fs);
Sy_pxx          = avg_welch_y(Sy, roi, win, nover, Fs);

% --- spectrogram of eta_bar ---
[Seta, f_eta, t_eta] = spectrogram(eta_bar - mean(eta_bar), hann(win), nover, [], Fs);
P_eta = abs(Seta).^2;

% Time axis for the per-frame Welch (use full record window only -> single column)
% Build a moving-window spectrogram for Sx, Sy too.
[Sx_spec, ~, t_spec] = spectrogram(squeeze(mean(mean(Sx(roi.row_min:roi.row_max,roi.col_min:roi.col_max,:),1),2)) - 0, ...
                                    hann(win), nover, [], Fs);
[Sy_spec, ~, ~]      = spectrogram(squeeze(mean(mean(Sy(roi.row_min:roi.row_max,roi.col_min:roi.col_max,:),1),2)) - 0, ...
                                    hann(win), nover, [], Fs);
P_Sx = abs(Sx_spec).^2;
P_Sy = abs(Sy_spec).^2;

t_off = time(1);
fig = figure('Name','A2 spectrograms (CSIG)','Position',[100 100 1100 1000],'Color','w');
tiledlayout(3,1,'TileSpacing','compact','Padding','compact');

ax1 = nexttile;
imagesc(t_spec + t_off, f_psd, log10(P_Sx + eps));
set(gca,'YDir','normal'); colormap(ax1, parula); cb = colorbar; cb.Label.String='log_{10}|S_x|^2';
ylabel('\omega/2\pi (Hz)'); title('|S_x(\omega,t)|^2');

ax2 = nexttile;
imagesc(t_spec + t_off, f_psd, log10(P_Sy + eps));
set(gca,'YDir','normal'); colormap(ax2, parula); cb = colorbar; cb.Label.String='log_{10}|S_y|^2';
ylabel('\omega/2\pi (Hz)'); title('|S_y(\omega,t)|^2');

ax3 = nexttile;
imagesc(t_eta + t_off, f_eta, log10(P_eta + eps));
set(gca,'YDir','normal'); colormap(ax3, parula); cb = colorbar; cb.Label.String='log_{10}|\eta|^2';
% peak frequency ridge
[~, ip] = max(P_eta, [], 1);
hold on; plot(t_eta + t_off, f_eta(ip), 'k:', 'LineWidth', 1.2);
ylabel('\omega/2\pi (Hz)'); xlabel('t (s)');
title('|\eta(\omega,t)|^2 (spatial mean) — dotted: \omega_p(t)');

sgtitle('A.2  Frequency spectrograms');
% save_figure(fig, 'A2_spectrograms_csig');
end


function [pxx, f] = avg_welch_y(S, roi, win, nover, Fs)
% Averaged Welch PSD across y rows, returned as a single column [Nf x 1].
[~, ~, Nt] = size(S);
% Sample at a few y rows (faster) and average.
y_idx = round(linspace(roi.row_min, roi.row_max, 8));
x_mid = round((roi.col_min + roi.col_max)/2);
acc = [];
for ii = 1:length(y_idx)
    s = squeeze(double(S(y_idx(ii), x_mid, :)));
    [p, f] = pwelch(s - mean(s), hann(win), nover, [], Fs);
    if isempty(acc), acc = zeros(length(f), 1); end
    acc = acc + p;
end
pxx = acc / length(y_idx);
end
