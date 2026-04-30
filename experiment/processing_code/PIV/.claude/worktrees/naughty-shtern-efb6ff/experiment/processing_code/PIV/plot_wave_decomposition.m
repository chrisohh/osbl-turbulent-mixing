function plot_wave_decomposition(u_cart, w_cart, u_wf, w_wf, ...
    u_phaseAvg, w_phaseAvg, u_ensembleAvg, w_ensembleAvg, ...
    u_turb, w_turb, u_wave, w_wave, ...
    transfo, zeta_vec, z_piv, x_piv, t_snap, ...
    Nbins, phase_edges)

phase_centers = 0.5 * (phase_edges(1:end-1) + phase_edges(2:end));

%% ===== FIGURE 1: Cartesian snapshot (raw velocity field) =====
figure('Name', sprintf('1. Cartesian velocity field  t=%.1fs', t_snap), ...
       'Position', [50 400 1200 400]);

subplot(1,2,1);
pcolor(x_piv*1e2, z_piv*1e2, u_cart); shading flat;
hold on;
plot(x_piv*1e2, transfo.eta*1e2, 'k', 'LineWidth', 1.5);
colorbar; colormap(gca, bluewhitered_cmap(256));
xlabel('x (cm)'); ylabel('z (cm)');
title('u (along-wind)  [m/s]');
set(gca, 'YDir', 'normal');
caxis_symmetric(gca);

subplot(1,2,2);
pcolor(x_piv*1e2, z_piv*1e2, w_cart); shading flat;
hold on;
plot(x_piv*1e2, transfo.eta*1e2, 'k', 'LineWidth', 1.5);
colorbar; colormap(gca, bluewhitered_cmap(256));
xlabel('x (cm)'); ylabel('z (cm)');
title('w (vertical)  [m/s]');
set(gca, 'YDir', 'normal');
caxis_symmetric(gca);

sgtitle(sprintf('Raw Cartesian velocity — t = %.1f s', t_snap));

%% ===== FIGURE 2: Wave-following coordinate system =====
figure('Name', sprintf('2. Wave-following coordinates  t=%.1fs', t_snap), ...
       'Position', [50 400 1200 400]);

subplot(1,2,1);
pcolor(x_piv*1e2, -zeta_vec*1e2, u_wf); shading flat;
colorbar; colormap(gca, bluewhitered_cmap(256));
xlabel('x (cm)'); ylabel('-\zeta (cm)   [0 = surface]');
title('u(x, \zeta)  [m/s]');
set(gca, 'YDir', 'normal');
caxis_symmetric(gca);

subplot(1,2,2);
pcolor(x_piv*1e2, -zeta_vec*1e2, w_wf); shading flat;
colorbar; colormap(gca, bluewhitered_cmap(256));
xlabel('x (cm)'); ylabel('-\zeta (cm)');
title('w(x, \zeta)  [m/s]');
set(gca, 'YDir', 'normal');
caxis_symmetric(gca);

sgtitle(sprintf('Wave-following coordinates — t = %.1f s  (top = surface, flattened)', t_snap));

% Also show the coordinate grid overlay on the Cartesian frame
figure('Name', '2b. Coordinate grid overlay', 'Position', [100 350 800 500]);
pcolor(x_piv*1e2, z_piv*1e2, w_cart); shading flat;
colorbar; colormap(bluewhitered_cmap(256));
hold on;
% Plot constant-zeta lines
zeta_show = [0, 0.002, 0.005, 0.01, 0.02, 0.05, 0.10];
colors_zeta = winter(length(zeta_show));
for k = 1:length(zeta_show)
    [~, iz] = min(abs(zeta_vec - zeta_show(k)));
    plot(x_piv*1e2, transfo.Z_grid(iz,:)*1e2, '-', ...
        'Color', colors_zeta(k,:), 'LineWidth', 1.2);
end
plot(x_piv*1e2, transfo.eta*1e2, 'k', 'LineWidth', 2);
xlabel('x (cm)'); ylabel('z (cm)');
title('Constant-\zeta lines on Cartesian field (black = surface)');
set(gca, 'YDir', 'normal');
legend([{'w field', '\eta(x)'}; ...
    arrayfun(@(z) sprintf('\\zeta = %.0f mm', z*1e3), zeta_show', 'Uni', 0)], ...
    'Location', 'southwest', 'FontSize', 7);
caxis_symmetric(gca);

%% ===== FIGURE 3: Phase vs x =====
figure('Name', '3. Wave phase vs x', 'Position', [100 300 900 300]);

plot(x_piv*1e2, transfo.phase, 'b-', 'LineWidth', 1);
hold on;
% Overlay surface elevation (scaled) for reference
eta_scaled = transfo.eta / max(abs(transfo.eta)) * pi;
plot(x_piv*1e2, eta_scaled, 'r--', 'LineWidth', 1);
xlabel('x (cm)');
ylabel('Phase (rad)');
title('Wave phase from Hilbert transform');
legend('\phi(x)', '\eta(x)  scaled', 'Location', 'best');
ylim([-pi pi]);
set(gca, 'YTick', [-pi -pi/2 0 pi/2 pi], ...
    'YTickLabel', {'-\pi', '-\pi/2', '0', '\pi/2', '\pi'});
grid on;

%% ===== FIGURE 4: Phase-averaged fields =====
figure('Name', '4. Phase-averaged velocity', 'Position', [50 250 1200 400]);

subplot(1,2,1);
pcolor(phase_centers, -zeta_vec*1e2, u_phaseAvg); shading flat;
colorbar; colormap(gca, bluewhitered_cmap(256));
xlabel('Phase (rad)'); ylabel('-\zeta (cm)');
title('<u>(\zeta, \phi)  [m/s]');
set(gca, 'YDir', 'normal');
set(gca, 'XTick', [-pi -pi/2 0 pi/2 pi], ...
    'XTickLabel', {'-\pi', '-\pi/2', '0', '\pi/2', '\pi'});
caxis_symmetric(gca);

subplot(1,2,2);
pcolor(phase_centers, -zeta_vec*1e2, w_phaseAvg); shading flat;
colorbar; colormap(gca, bluewhitered_cmap(256));
xlabel('Phase (rad)'); ylabel('-\zeta (cm)');
title('<w>(\zeta, \phi)  [m/s]');
set(gca, 'YDir', 'normal');
set(gca, 'XTick', [-pi -pi/2 0 pi/2 pi], ...
    'XTickLabel', {'-\pi', '-\pi/2', '0', '\pi/2', '\pi'});
caxis_symmetric(gca);

sgtitle('Phase-averaged velocity (wave-following coords)');

%% ===== FIGURE 5: Ensemble average =====
figure('Name', '5. Ensemble-averaged profiles', 'Position', [150 200 600 450]);

subplot(1,2,1);
plot(u_ensembleAvg, -zeta_vec*1e2, 'b-', 'LineWidth', 1.5);
xlabel('<u>_{ens} (m/s)'); ylabel('-\zeta (cm)');
title('Ensemble-averaged u');
set(gca, 'YDir', 'normal');
grid on;

subplot(1,2,2);
plot(w_ensembleAvg, -zeta_vec*1e2, 'r-', 'LineWidth', 1.5);
xlabel('<w>_{ens} (m/s)'); ylabel('-\zeta (cm)');
title('Ensemble-averaged w');
set(gca, 'YDir', 'normal');
grid on;

sgtitle('Ensemble average (all phases, all frames)');

%% ===== FIGURE 6: Turbulent & wave-coherent decomposition =====
figure('Name', sprintf('6. Decomposition  t=%.1fs', t_snap), ...
       'Position', [50 50 1400 700]);

% --- Top row: turbulent (u' = u - phase avg) ---
subplot(2,2,1);
pcolor(x_piv*1e2, -zeta_vec*1e2, u_turb); shading flat;
colorbar; colormap(gca, bluewhitered_cmap(256));
xlabel('x (cm)'); ylabel('-\zeta (cm)');
title("u' = u - <u>_\phi   (turbulent)");
set(gca, 'YDir', 'normal');
caxis_symmetric(gca);

subplot(2,2,2);
pcolor(x_piv*1e2, -zeta_vec*1e2, w_turb); shading flat;
colorbar; colormap(gca, bluewhitered_cmap(256));
xlabel('x (cm)'); ylabel('-\zeta (cm)');
title("w' = w - <w>_\phi   (turbulent)");
set(gca, 'YDir', 'normal');
caxis_symmetric(gca);

% --- Bottom row: wave-coherent (u~ = phase avg - ensemble avg) ---
subplot(2,2,3);
pcolor(x_piv*1e2, -zeta_vec*1e2, u_wave); shading flat;
colorbar; colormap(gca, bluewhitered_cmap(256));
xlabel('x (cm)'); ylabel('-\zeta (cm)');
title('u~ = <u>_\phi - <u>_{ens}   (wave-coherent)');
set(gca, 'YDir', 'normal');
caxis_symmetric(gca);

subplot(2,2,4);
pcolor(x_piv*1e2, -zeta_vec*1e2, w_wave); shading flat;
colorbar; colormap(gca, bluewhitered_cmap(256));
xlabel('x (cm)'); ylabel('-\zeta (cm)');
title('w~ = <w>_\phi - <w>_{ens}   (wave-coherent)');
set(gca, 'YDir', 'normal');
caxis_symmetric(gca);

sgtitle(sprintf('Triple decomposition — t = %.1f s', t_snap));

end


%% ========== HELPER FUNCTIONS ==========

function cmap = bluewhitered_cmap(n)
%BLUEWHITERED_CMAP  Diverging blue-white-red colormap.
if nargin < 1, n = 256; end
half = floor(n/2);
r = [linspace(0, 1, half), ones(1, n-half)];
g = [linspace(0, 1, half), linspace(1, 0, n-half)];
b = [ones(1, half), linspace(1, 0, n-half)];
cmap = [r' g' b'];
end

function caxis_symmetric(ax)
%CAXIS_SYMMETRIC  Set symmetric color limits around zero.
cl = caxis(ax);
mx = max(abs(cl));
if mx > 0
    caxis(ax, [-mx mx]);
end
end
