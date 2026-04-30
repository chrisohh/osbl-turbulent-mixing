function plot_compare_methods(stats_A, stats_B, stats_C)
%PLOT_COMPARE_METHODS  Compare wave removal methods A (Hilbert), B (y-mean), C (phase avg).

z = stats_A.z;
t = stats_A.t;

times = [40 55 70 90];
colors = lines(3);
labels = {'A: Hilbert', 'B: y-mean', 'C: phase avg'};

%% --- <w'w'> vertical profiles at selected times ---
figure('Name','Compare methods: <w''w''> profiles');
for ip = 1:length(times)
    [~, it] = min(abs(t - times(ip)));

    subplot(1, length(times), ip); hold on;
    plot(stats_A.ww(:,it), z, '-',  'Color', colors(1,:), 'LineWidth', 1.5);
    plot(stats_B.ww(:,it), z, '--', 'Color', colors(2,:), 'LineWidth', 1.5);
    plot(stats_C.ww(:,it), z, ':',  'Color', colors(3,:), 'LineWidth', 1.5);
    set(gca, 'YDir', 'normal');
    xlabel('<w''w''>');
    if ip == 1, ylabel('z (m)'); end
    title(sprintf('t = %.0f s', t(it)));
    if ip == length(times), legend(labels, 'Location', 'best'); end
end

%% --- <v'w'> vertical profiles ---
figure('Name','Compare methods: <v''w''> profiles');
for ip = 1:length(times)
    [~, it] = min(abs(t - times(ip)));

    subplot(1, length(times), ip); hold on;
    plot(stats_A.vw(:,it), z, '-',  'Color', colors(1,:), 'LineWidth', 1.5);
    plot(stats_B.vw(:,it), z, '--', 'Color', colors(2,:), 'LineWidth', 1.5);
    plot(stats_C.vw(:,it), z, ':',  'Color', colors(3,:), 'LineWidth', 1.5);
    set(gca, 'YDir', 'normal');
    xlabel('<v''w''>');
    if ip == 1, ylabel('z (m)'); end
    title(sprintf('t = %.0f s', t(it)));
    if ip == length(times), legend(labels, 'Location', 'best'); end
end

%% --- Hovmoller comparison: <w'w'>(z,t) for all three methods ---
figure('Name','Compare methods: Hovmoller <w''w''>');
clim_max = max([stats_A.ww(:); stats_B.ww(:); stats_C.ww(:)]);

subplot(1,3,1);
imagesc(t, z, stats_A.ww); set(gca,'YDir','normal');
colorbar; caxis([0 clim_max]);
xlabel('t (s)'); ylabel('z (m)'); title('A: Hilbert');

subplot(1,3,2);
imagesc(t, z, stats_B.ww); set(gca,'YDir','normal');
colorbar; caxis([0 clim_max]);
xlabel('t (s)'); title('B: y-mean');

subplot(1,3,3);
imagesc(t, z, stats_C.ww); set(gca,'YDir','normal');
colorbar; caxis([0 clim_max]);
xlabel('t (s)'); title('C: phase avg');

end
