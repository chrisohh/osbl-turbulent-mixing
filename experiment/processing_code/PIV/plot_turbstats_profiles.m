function plot_turbstats_profiles(stats)

z = stats.z;
t = stats.t;

figure;
imagesc(t, z, stats.ww);
set(gca,'YDir','normal');
colorbar;
xlabel('t (s)');
ylabel('z (m)');
title('<w''w''> (z,t)');
caxis([0 max(stats.ww(:))]);

figure;
hold on;

% Times within Ramp 2 PIV window (t_JFM = 38–98 s, fully within LC turbulence phase).
% JFM stages 1–3 (viscous, LC onset, self-sharpening: t < 20 s) precede PIV recording.
times = [40 55 70 90];

for i = 1:length(times)
    [~,it] = min(abs(t - times(i)));
    plot(stats.ww(:,it), z, 'DisplayName', sprintf('t=%.1f',t(it)));
end

set(gca,'YDir','reverse');
xlabel('<w''w''>');
ylabel('z');
legend;
title('Vertical profiles');

end