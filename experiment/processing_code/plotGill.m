load('G:\My Drive\OSBL\GILL\20250923_Gill_SOARS_Data.mat')
%% Plot raw
figure;
subplot(2,1,1);hold on;plot(GILL.time_str,GILL.ux);plot(GILL.time_str,GILL.uy);plot(GILL.time_str,GILL.uz)
subplot(2,1,2);hold on;plot(GILL.time_str,GILL.temp)

subplot(2,1,1);set(gca,'fontsize',12),xlabel('Time'),ylabel('Velocity (m/s)')
subplot(2,1,2);set(gca,'fontsize',12),xlabel('Time'),ylabel('Sonic Temperature (\circC)')

%%
duration_sec=(datenum(GILL.time_str)-datenum(GILL.time_str(1)))*24*3600;
time_rampup=1:16196;
figure;hold on
plot(duration_sec(time_rampup),GILL.ux(time_rampup))
%%
figure;subplot(2,1,1);plot(duration_sec,GILL.ux);
xlabel('time');ylabel('u_x(m/s)')
subplot(2,1,2);plot(duration_sec,GILL.temp);
xlabel('time');ylabel('T (^\circC)')

%%
% Take steady-state portion only
time_raw=duration_sec;
fs = 1/mean(diff(time_raw),'omitnan');
[~,start_idx]=min(abs(time_raw-0));
[~,end_idx]=min(abs(time_raw-130));
steady_state_idx=[start_idx:end_idx];
U_steady=GILL.ux(start_idx:end_idx);
V_steady=GILL.uy(start_idx:end_idx);
W_steady=GILL.uz(start_idx:end_idx);
T_steady=GILL.temp(start_idx:end_idx);
time_steady = time_raw(steady_state_idx);
n=6;
sonicMean.RPM(n)=50;
sonicMean.U(n)=mean(U_steady,'omitnan');sonicMean.V(n)=mean(V_steady,'omitnan');sonicMean.W(n)=mean(W_steady,'omitnan');

%%
RPM=[sonicMean.RPM(6),sonicMean.RPM(1:end-2)];
U_Mag=sqrt([sonicMean.U(6),sonicMean.U(1:end-2)].^2+[sonicMean.U(6),sonicMean.U(1:end-2)].^2+[sonicMean.U(6),sonicMean.U(1:end-2)].^2);
figure;plot(RPM,U_Mag,'k-o')
xlabel('RPM');ylabel('U')
%% Linear fit
coefficients = polyfit(RPM, U_Mag, 1);
a = coefficients(1);  % slope
b = coefficients(2);  % intercept
fprintf('Linear fit: U = %.4f * rpm + %.4f\n', a, b);