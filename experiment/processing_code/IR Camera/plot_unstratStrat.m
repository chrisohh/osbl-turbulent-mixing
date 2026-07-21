%% Plot unstratified vs stratified SST anomaly (ΔT)
% Loads pre-processed .mat files saved by distortion_correction.m
% Rec-000014: unstratified
% Rec-000026: stratified, 10 cm freshwater layer

configs = {
    'D:\HLAB_2026\IR_camera\Processed\IR_Rec-000014_params.mat', 'Unstratified';
    'D:\HLAB_2026\IR_camera\Processed\IR_Rec-000026_params.mat', 'Stratified (10 cm FW)';
};

clim_sym = 0.3;  % ± °C — adjust after first run

figure('Position', [50 50 1100 520], 'Color', 'w');

for k = 1:2
    d = load(configs{k,1});
    % dT is already computed and saved by distortion_correction.m
    dT = d.dT;

    subplot(1, 2, k);
    imagesc(d.x*100, d.y*100, dT);

    % Diverging colormap: cmocean('balance') if available, else inline BWR
    if exist('cmocean', 'file')
        cmocean('balance');
    else
        colormap(bwr_inline());
    end

    c = colorbar;
    clim([-clim_sym, clim_sym]);
    axis equal tight;

    title(sprintf('%s\n$t=%.2f$ s', configs{k,2}, d.frameToPlot/d.frameRate), ...
        'Interpreter', 'latex', 'FontSize', 14);
    xlabel('Along-wind (cm)', 'FontSize', 13);
    if k == 1
        ylabel('Cross-wind (cm)', 'FontSize', 13);
    end
    c.Label.String = '\DeltaT (°C)';
    set(gca, 'FontSize', 12, 'FontName', 'Times');
end

sgtitle('SST Anomaly \DeltaT — Unstratified vs Stratified', ...
    'FontSize', 15, 'FontWeight', 'bold');

%% -----------------------------------------------------------------------
function cmap = bwr_inline()
% Blue-white-red diverging colormap (fallback if cmocean not installed)
    n = 256; h = n/2;
    cmap = [linspace(0,1,h)', linspace(0,1,h)', ones(h,1); ...
            ones(h,1), linspace(1,0,h)', linspace(1,0,h)'];
end
