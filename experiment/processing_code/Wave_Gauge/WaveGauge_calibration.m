% Calibration data collected
% calibration.date = datetime(2025,12,22,12,45,00);
% height_zero=147;
% height_calib=[160,170,180,190,150,140,130,120,110,104];
% voltage_calib=[0,1.6,4.2,7.9,-0.92,-1.7,-2.3,-2.8,-3.22,-3.5];

calibration.date = datetime(2025,12,22,14,38,00);
height_zero=147;
height_calib=[161:-2:133];
voltage_calib=[7.3,6,4.6,3,1.4,0,-1.4,-2.6,-3.8,-4.9,-6,-7.1,-8.2,-9.1,-10];


[height_calib,idx]=sort(height_calib);
voltage_calib=voltage_calib(idx);

rel_height_calib=-(height_calib-height_zero);

% Use cubic polynomial
p = polyfit(voltage_calib, rel_height_calib, 3);
height_fit = polyval(p, voltage_calib);

residuals = rel_height_calib - height_fit;

fprintf('Max residual: %.2f cm\n', max(abs(residuals)));
fprintf('RMS residual: %.2f cm\n', rms(residuals));

% plot the fit
figure;
scatter(voltage_calib, rel_height_calib, 'filled')
hold on
plot(voltage_calib, height_fit, 'r-', 'LineWidth', 2)
xlabel('Voltage (V)')
ylabel('Translated Height (cm)')
legend('Calibration data', 'Cubic fit')
set(gca,'fontsize',14,'fontname','times')

%% Save calibration data
calibration.height_zero = height_zero;
calibration.height_calib = height_calib;
calibration.voltage_calib = voltage_calib;
calibration.rel_height_calib = rel_height_calib;
calibration.polynomial_coefficients = p;
calibration.polynomial_order = 3;
calibration.max_residual_cm = max(abs(residuals));
calibration.rms_residual_cm = rms(residuals);

calibration_date=datestr(calibration.date, 'yyyymmdd_HHMM');
save(strcat('G:\Shared drives\AirSeaLab\Projects\SOARS\Data\20251222\WaveGauge\',...
    sprintf('wave_gauge_calibration_%s.mat',calibration_date)), 'calibration');
fprintf('Calibration saved to: wave_gauge_calibration_%s.mat\n', calibration_date);

