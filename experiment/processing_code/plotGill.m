load('G:\My Drive\OSBL\GILL\20250923_Gill_SOARS_Data.mat')
duration_sec=(datenum(GILL.time_str)-datenum(GILL.time_str(1)))*24*3600;
time_rampup=1:16196;
figure;hold on
plot(duration_sec(time_rampup),GILL.ux(time_rampup))
