function cisg_plot_results(Sx_true, Sy_true, Sx_meas, Sy_meas, eta, x, y)
% CISG_PLOT_RESULTS Visualize slope measurement results
%
% Inputs:
%   Sx_true, Sy_true - true slopes (if available, else pass [])
%   Sx_meas, Sy_meas - measured slopes
%   eta - wave elevation (if available, else pass [])
%   x, y - coordinate vectors

    have_truth = ~isempty(Sx_true);
    have_elevation = ~isempty(eta);
    
    if have_truth
        % Full comparison plot
        figure('Name', 'CISG Results', 'Position', [100 100 1400 900]);
        
        % True Sx
        subplot(3,3,1);
        imagesc(x, y, Sx_true);
        axis equal tight; colorbar;
        title('True S_x');
        xlabel('X (cm)'); ylabel('Y (cm)');
        colormap(gca, 'jet');
        
        % Measured Sx
        subplot(3,3,2);
        imagesc(x, y, Sx_meas);
        axis equal tight; colorbar;
        title('Measured S_x');
        xlabel('X (cm)'); ylabel('Y (cm)');
        colormap(gca, 'jet');
        
        % Error Sx
        subplot(3,3,3);
        Sx_err = Sx_meas - Sx_true;
        imagesc(x, y, Sx_err);
        axis equal tight; colorbar;
        title(sprintf('Error S_x (RMS=%.4f)', std(Sx_err(:))));
        xlabel('X (cm)'); ylabel('Y (cm)');
        colormap(gca, 'jet');
        
        % True Sy
        subplot(3,3,4);
        imagesc(x, y, Sy_true);
        axis equal tight; colorbar;
        title('True S_y');
        xlabel('X (cm)'); ylabel('Y (cm)');
        colormap(gca, 'jet');
        
        % Measured Sy
        subplot(3,3,5);
        imagesc(x, y, Sy_meas);
        axis equal tight; colorbar;
        title('Measured S_y');
        xlabel('X (cm)'); ylabel('Y (cm)');
        colormap(gca, 'jet');
        
        % Error Sy
        subplot(3,3,6);
        Sy_err = Sy_meas - Sy_true;
        imagesc(x, y, Sy_err);
        axis equal tight; colorbar;
        title(sprintf('Error S_y (RMS=%.4f)', std(Sy_err(:))));
        xlabel('X (cm)'); ylabel('Y (cm)');
        colormap(gca, 'jet');
        
        % Slope magnitude scatter
        subplot(3,3,7);
        slope_mag_true = sqrt(Sx_true.^2 + Sy_true.^2);
        slope_mag_meas = sqrt(Sx_meas.^2 + Sy_meas.^2);
        plot(slope_mag_true(:), slope_mag_meas(:), '.', 'MarkerSize', 2);
        hold on;
        lims = [0 max(slope_mag_true(:))];
        plot(lims, lims, 'r--', 'LineWidth', 2);
        xlabel('True |S|');
        ylabel('Measured |S|');
        title('Slope Magnitude');
        axis equal; grid on;
        
        % Error histogram
        subplot(3,3,8);
        histogram(Sx_err(:), 50, 'Normalization', 'pdf');
        hold on;
        histogram(Sy_err(:), 50, 'Normalization', 'pdf');
        xlabel('Slope Error');
        ylabel('PDF');
        legend('S_x', 'S_y');
        title('Error Distribution');
        grid on;
        
        % Wave elevation
        if have_elevation
            subplot(3,3,9);
            imagesc(x, y, eta);
            axis equal tight; colorbar;
            title(sprintf('Elevation (RMS=%.3f cm)', std(eta(:))));
            xlabel('X (cm)'); ylabel('Y (cm)');
            colormap(gca, 'jet');
        end
        
    else
        % Simple visualization without ground truth
        figure('Name', 'CISG Measured Slopes', 'Position', [100 100 1000 400]);
        
        subplot(1,3,1);
        imagesc(x, y, Sx_meas);
        axis equal tight; colorbar;
        title('Measured S_x');
        xlabel('X (cm)'); ylabel('Y (cm)');
        colormap(gca, 'jet');
        
        subplot(1,3,2);
        imagesc(x, y, Sy_meas);
        axis equal tight; colorbar;
        title('Measured S_y');
        xlabel('X (cm)'); ylabel('Y (cm)');
        colormap(gca, 'jet');
        
        subplot(1,3,3);
        slope_mag = sqrt(Sx_meas.^2 + Sy_meas.^2);
        imagesc(x, y, slope_mag);
        axis equal tight; colorbar;
        title('Slope Magnitude |S|');
        xlabel('X (cm)'); ylabel('Y (cm)');
        colormap(gca, 'jet');
    end
end
