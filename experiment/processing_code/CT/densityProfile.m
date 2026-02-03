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
legend('Strat 5cm Fresh, Before Wind', 'Strat 10cm Fresh, Before Wind','Strat 10cm Fresh, After Wind','location','southwest')
plot([0,2],[5,5],'--')
plot([0,2],[10,10],'--')
ylim([0,20])
text(1.5,7.5,'Fresh 5cm')
text(1.25,3,'Additional Fresh Water 5cm')
text(0.75,15.5,'Salt Water')