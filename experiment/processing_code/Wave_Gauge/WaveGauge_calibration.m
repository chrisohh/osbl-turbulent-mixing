height_zero=147;
height_calib=[160,170,180,190,150,140,130,120,110,104];
voltage_calib=[0,1.6,4.2,7.9,-0.92,-1.7,-2.3,-2.8,-3.22,-3.5];
[height_calib,idx]=sort(height_calib);
voltage_calib=voltage_calib(idx);
% height_calib=[161:-2:133];
% voltage_calib=[7.3,6,4.6,3,1.4,0,-1.4,-2.6,-3.8,-4.9,-6,-7.1,-8.2,-9.1,-10];

rel_height_calib=-(height_calib-height_zero);

figure;scatter(voltage_calib,rel_height_calib)
p=polyfit(voltage_calib,rel_height_calib,2);
height_fit=p(1)*voltage_calib.^2+p(2)*voltage_calib+p(3);
hold on
plot(voltage_calib,height_fit)

% wave_height_cm=(height-height_zero,voltage_calib,2);

%%
load('wave_gauge_20251222_134754.mat')
voltage=save_data.voltage;
time=save_data.time;
wave_height_cm=p(1)*voltage.^2+p(2)*voltage+p(3);
figure;plot(time,wave_height_cm)