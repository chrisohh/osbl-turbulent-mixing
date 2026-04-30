%% 
clear all
folder = 'G:\Shared drives\AirSeaLab\Projects\SOARS\Data\20251222\WaveGauge\';
calibrationFile = 'wave_gauge_calibration_20251222_1438.mat';

% amplitude=[0.2,]
dataFile = 'wave_gauge_20251222_144838.mat';
load(strcat(folder,calibrationFile))
load(strcat(folder,dataFile))
voltage=save_data.voltage;
time=save_data.time;
wave_height_cm=polyval(calibration.polynomial_coefficients, voltage);
figure;plot(time,wave_height_cm)
xlabel('Time (s)')
ylabel('Wave Height (cm)')
set(gca,'fontsize',14,'fontname','times')
%% Chop one time series into increments
fs = 1000;
time_start_index=[0:2:10]*60*fs+1;
figure;hold on;
for i=1:length(time_start_index)-1
    range=(time_start_index(i):time_start_index(i+1)-1);
    plot(time(range),wave_height_cm(range))
end
    xlabel('Time (s)')
    ylabel('Wave Height (cm)')
    set(gca,'fontsize',14,'fontname','times')
legend('50rpm','200rpm','400rpm','600rpm','800rpm')
%% plot all
files = dir(strcat(folder,'wave_gauge_*.mat'));
file_index=[8:13];
legend_str=[50:10:100];

figure;
for i = 1:length(file_index)
    load(strcat(folder,files(file_index(i)).name));

    % Extract your data (adjust variable names as needed)
    voltage = save_data.voltage;
    time = save_data.time;
    wave_height_cm = polyval(calibration.polynomial_coefficients, voltage);

    subplot(length(file_index), 1, i)
    hold on
    plot(time, wave_height_cm)
    % title(files(i).name, 'Interpreter', 'none')
    if i==length(file_index)
    xlabel('Time (s)')
    else
        set(gca,'xticklabel',[])
    end
    if i==round(length(file_index)/2)
        ylabel('Height (cm)')
    end
    legend(sprintf('%d%%',legend_str(i)))
    ylim([-30,30])
    set(gca,'fontsize',14,'fontname','times')

    %save
    % save(strcat(folder,sprintf('sin-sin_A%.1f_P1.mat',legend_str(i)/100)),'wave_height_cm','time')
end

%%
analytic_signal = hilbert(wave_height_cm);
envelope = abs(analytic_signal);
plot(time, wave_height_cm, time, envelope, 'r--', time, -envelope, 'r--')

%%
figure; hold on
for i=1:length(time_start_index)-1
    range=(time_start_index(i):time_start_index(i+1)-1);
    [Pxx, f] = pwelch(wave_height_cm(range), [], [], [], fs);
    plot(f, Pxx, 'LineWidth', 1)
    set(gca,'xscale','log','yscale','log')
    xlabel('Frequency (Hz)')
    ylabel('PSD (cm^2/Hz)')

    time_centers(i) = time(round((range(1)+range(end))/2));
    % [pks,idx]=findpeaks(wave_height_cm(range));
    % H_rms(i) = sqrt(mean(pks.^2)); 

    [~, idx_peak] = max(Pxx(f>1));
    f_above1=f(f>1);
    peak_freq(i) = f_above1(idx_peak);
end
legend('50rpm=0.25m/s','200rpm=2.5m/s','400rpm=5.5m/s','600rpm=8.5m/s','800rpm=11.5m/s')

% figure;
% subplot(2,1,1)
% plot(time_centers, H_rms)
% ylabel('H_{rms} (cm)')
% title('Wave Statistics Evolution')
% 
% subplot(2,1,2)
% plot(time_centers, peak_freq)
% ylabel('f_{peak} (Hz)')
% xlabel('Time (s)')

%% Calculate wave slope (important for breaking)
figure; hold on
for i=1:length(time_start_index)-1
    range=(time_start_index(i):time_start_index(i+1)-1);
wave_slope(i,:) = gradient(wave_height_cm(range), mean(diff(time(range))));
time_range(i,:)=time(range);
end
figure; hold on
for i=1:length(time_start_index)-1
histogram(wave_slope(i,:), 50)
end
xlabel('Wave Slope (cm/s)')
ylabel('Count')
title('Wave Slope Distribution')

% Or plot vs time
figure;hold on
for i=1:length(time_start_index)-1
plot(time_range(i,:), abs(wave_slope(i,:)))
end
ylabel('|d\eta/dt| (cm/s)')
xlabel('Time (s)')
