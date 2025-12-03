load('D:\Scripps\osbl-turbulent-mixing\experiment\data\251023\20251023_Gill_SOARS_Data.mat')
%% Plot raw
figure;
subplot(2,1,1);hold on;plot(GILL.time_str,GILL.ux);plot(GILL.time_str,GILL.uy);plot(GILL.time_str,GILL.uz)
subplot(2,1,2);hold on;plot(GILL.time_str,GILL.temp)

subplot(2,1,1);set(gca,'fontsize',12),xlabel('Time'),ylabel('Velocity (m/s)')
subplot(2,1,2);set(gca,'fontsize',12),xlabel('Time'),ylabel('Sonic Temperature (\circC)')

%%
duration_sec=(datenum(GILL.time_str)-datenum(GILL.time_str(1)))*24*3600;
time_rampup=1:length(duration_sec);
figure;hold on
plot(duration_sec(time_rampup),GILL.ux(time_rampup))
%%
figure;subplot(2,1,1);plot(duration_sec,GILL.ux);
xlabel('time');ylabel('u_x(m/s)')
subplot(2,1,2);plot(duration_sec,GILL.temp);
xlabel('time');ylabel('T (^\circC)')

%%
% Take steady-state portion only
rpm=[50,200,400,600,800];
start_time=[17,1359,1518,1659,1873];
end_time=[1303,1492,1643,1855,2548];
for n=1:length(rpm)
time_raw=duration_sec;
fs = 1/mean(diff(time_raw),'omitnan');
[~,start_idx]=min(abs(time_raw-start_time(n)));
[~,end_idx]=min(abs(time_raw-end_time(n)));
steady_state_idx=[start_idx:end_idx];
U_steady=GILL.ux(start_idx:end_idx);
V_steady=GILL.uy(start_idx:end_idx);
W_steady=GILL.uz(start_idx:end_idx);
T_steady=GILL.temp(start_idx:end_idx);
time_steady = time_raw(steady_state_idx);

sonicMean.RPM(n)=rpm(n);
sonicMean.U(n)=mean(U_steady,'omitnan');sonicMean.V(n)=mean(V_steady,'omitnan');sonicMean.W(n)=mean(W_steady,'omitnan');
end
%%
save('D:\Scripps\osbl-turbulent-mixing\experiment\data\250923\rpm2u.mat','start_time','end_time','sonicMean')
%%
figure;plot(sonicMean.RPM,sonicMean.U,'k-o')
hold on;plot(sonicMean.RPM,sonicMean.V,'b-o')
hold on;plot(sonicMean.RPM,sonicMean.W,'r-o')
xlabel('RPM');ylabel('U')
%% Linear fit
coefficients = polyfit(sonicMean.RPM, sonicMean.U, 1);
a = coefficients(1);  % slopef
b = coefficients(2);  % intercept
fprintf('Linear fit: U = %.4f * rpm + %.4f\n', a, b);