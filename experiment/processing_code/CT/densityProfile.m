load('C:\Users\asila\OneDrive\Documents\osbl-turbulent-mixing\experiment\processing_code\CT\StratificationVar.mat')

%%
figure;
plot(voltThurs/10,adjusted_depth5cm,'linewidth',1.5)
set(gca,'Ydir','reverse')
hold on
plot(voltFri,adjusted_depth10cm,'linewidth',1.5)
plot(voltEODFriday,adjusted_depth10cm,'linewidth',1.5)
xlabel('Voltage')
ylabel('z (cm)')
plot([0,2],[5,5],'--')
plot([0,2],[10,10],'--')
legend('Strat with 5cm FW, Before Wind', 'Strat with 10cm FW, Before Wind','Strat with 10cm FW, After Wind','','','location','southwest')
ylim([0,20])
text(1.5,7.5,'Fresh 5cm','fontsize',12,'fontname','times')
text(1,3,'Additional Fresh Water 5cm','fontsize',12,'fontname','times')
text(0.75,12.5,'Salt Water','fontsize',12,'fontname','times')
set(gca,'fontsize',14,'fontname','times')
set(gcf,'color','white')