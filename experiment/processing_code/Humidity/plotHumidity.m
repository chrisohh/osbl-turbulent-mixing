load('20251223_Gill_SOARS_Data.mat')
%%
[AH] = humidity_conversion(GILL.rh,GILL.temp);
%% Plot raw
figure;
subplot(3,1,1);hold on;plot(GILL.time_str,GILL.temp);
subplot(3,1,2);hold on;plot(GILL.time_str,GILL.rh)
subplot(3,1,3);hold on;plot(GILL.time_str,AH)

subplot(3,1,1);set(gca,'fontsize',14,'fontname','times'),xlabel('Time'),ylabel('Temperature (C^\circ)')
subplot(3,1,2);set(gca,'fontsize',14,'fontname','times'),xlabel('Time'),ylabel('Relative humidity %')
subplot(3,1,3);set(gca,'fontsize',14,'fontname','times'),xlabel('Time'),ylabel('Absolute humidity (g/cm^3)')

subplot(3,1,1);xlim([datetime(2025,12,23,11,53,00),datetime(2025,12,23,12,20,00)])
subplot(3,1,2);xlim([datetime(2025,12,23,11,53,00),datetime(2025,12,23,12,20,00)])
subplot(3,1,3);xlim([datetime(2025,12,23,11,53,00),datetime(2025,12,23,12,20,00)])


%%
