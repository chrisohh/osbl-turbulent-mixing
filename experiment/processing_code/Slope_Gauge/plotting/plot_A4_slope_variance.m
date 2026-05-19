function plot_A4_slope_variance(data)
% A.4  Time series of <S_x^2>(t) and <S_y^2>(t), area-averaged in ROI.

Sx   = data.Sx;
Sy   = data.Sy;
time = data.time(:);
roi  = data.setup.roi;
Fs   = data.setup.frame_rate;

Nt = size(Sx, 3);
sx2 = zeros(Nt, 1);
sy2 = zeros(Nt, 1);
for it = 1:Nt
    sx_f = double(Sx(roi.row_min:roi.row_max, roi.col_min:roi.col_max, it));
    sy_f = double(Sy(roi.row_min:roi.row_max, roi.col_min:roi.col_max, it));
    sx2(it) = mean(sx_f(:).^2);
    sy2(it) = mean(sy_f(:).^2);
end

% 1-second rolling mean
nroll = max(1, round(Fs));
sx2_s = movmean(sx2, nroll);
sy2_s = movmean(sy2, nroll);

% First crossover Sy2/Sx2 > 0.5
ratio = sy2_s ./ max(sx2_s, eps);
ix_cross = find(ratio > 0.5, 1, 'first');

fig = figure('Name','A4 slope variance','Position',[100 100 1000 500],'Color','w');
semilogy(time, sx2_s, 'b-', 'LineWidth', 1.5); hold on;
semilogy(time, sy2_s, 'r-', 'LineWidth', 1.5);
if ~isempty(ix_cross)
    xline(time(ix_cross), 'k--', sprintf('S_y^2/S_x^2 > 0.5  @ t=%.1fs', time(ix_cross)), ...
          'LabelOrientation','horizontal');
end
xlabel('t (s)'); ylabel('<S^2>');
legend('<S_x^2>','<S_y^2>','Location','best');
title('A.4  Slope variance time series  (1-s rolling mean)');
grid on;
save_figure(fig, 'A4_slope_variance');
end
