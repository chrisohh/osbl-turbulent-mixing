% load('C:\Users\asila\OneDrive\Documents\osbl-turbulent-mixing\experiment\processing_code\CT\StratificationVar.mat')

clear; clc;
%%
DEV_ID = "Dev4";   % <-- confirm this matches NI MAX for the USB-6451
Duration = 5; %sec

d = daq("ni");
ch = addinput(d, DEV_ID, "ai7", "Voltage");
ch.TerminalConfig = "SingleEnded";
ch.Range = [-10 10];
d.Rate = 100;   % conservative; raise only if you know the probe's conditioner bandwidth
height = [];
voltage = [];
clf(1);clf(2)
%% Manual depth traverse
% Move probe so it just breaks the water surface = height 0 cm.
% At each depth, enter the height and press Enter to record a reading.
% Height increases as the probe goes deeper. Type "q" to stop.

figure(1);hold on
plot(voltage,height,'o','LineWidth',1.5)
set(gca,'Ydir','reverse')
xlabel('Voltage')
ylabel('z (cm)')
set(gca,'fontsize',14,'fontname','times')
set(gcf,'color','white')

while true
    resp = input('Height (cm), or "q" to stop: ', 's');
    if strcmpi(resp, 'q') || isempty(resp)
        break;
    end
    h = str2double(resp);
    if isnan(h)
        disp('Invalid input, enter a number or "q".');
        continue;
    end

    data = read(d, seconds(5));   % short averaged reading at this depth
    v = mean(data.Variables);
    figure(2);plot(data.Time,data.Variables)
    figure(1);plot(v,h,'o','LineWidth',1.5)
    height(end+1) = h;
    voltage(end+1) = v;
    fprintf('  recorded voltage = %.4f V at height = %.2f cm\n', v, h);

end

%%
h0=29.5;%when black box is submerged below the surface
depth=height-h0;
save(sprintf('densityProfileData_%s.mat',datestr(now,'yyyymmdd_HHMM')), 'depth', 'voltage');
%%

figure;hold on
plot(voltage,depth,'o','LineWidth',1.5)
set(gca,'Ydir','reverse')
xlabel('Voltage')
ylabel('z (cm)')
plot([-0.5,5],[10,10],'--')
text(1,3,'Fresh Water 10 cm','fontsize',12,'fontname','times')
text(0.75,12.5,'Salt Water','fontsize',12,'fontname','times')
set(gca,'fontsize',14,'fontname','times')
set(gcf,'color','white')
xlim([0,3])
% figure;
% plot(voltThurs/10,adjusted_depth5cm,'linewidth',1.5)
% set(gca,'Ydir','reverse')
% hold on
% plot(voltFri,adjusted_depth10cm,'linewidth',1.5)
% plot(voltEODFriday,adjusted_depth10cm,'linewidth',1.5)
% xlabel('Voltage')
% ylabel('z (cm)')
% plot([0,2],[5,5],'--')
% plot([0,2],[10,10],'--')
% legend('Strat with 5cm FW, Before Wind', 'Strat with 10cm FW, Before Wind','Strat with 10cm FW, After Wind','','','location','southwest')
% ylim([0,20])
% text(1.5,7.5,'Fresh 5cm','fontsize',12,'fontname','times')
% text(1,3,'Additional Fresh Water 5cm','fontsize',12,'fontname','times')
% text(0.75,12.5,'Salt Water','fontsize',12,'fontname','times')
% set(gca,'fontsize',14,'fontname','times')
% set(gcf,'color','white')