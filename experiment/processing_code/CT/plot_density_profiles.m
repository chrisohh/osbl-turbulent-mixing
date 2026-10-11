% Overlay the post-run density (conductivity voltage) profiles,
% starting with densityProfileData_20261007_2340.
clear; clc;

here = fileparts(mfilename('fullpath'));
if isempty(here), here = pwd; end

startStamp = "20261007_2340";   % first run to include

f = dir(fullfile(here,'densityProfileData_*.mat'));
stamps = string(regexprep({f.name},'^densityProfileData_|\.mat$',''));
keep = strlength(stamps) == 13 & stamps >= startStamp;
f = f(keep);  stamps = stamps(keep);
[stamps,ix] = sort(stamps);  f = f(ix);

figure('Color','white'); hold on
co = lines(numel(f));
for k = 1:numel(f)
    S = load(fullfile(here,f(k).name),'depth','voltage');
    t = datetime(stamps(k),'InputFormat','yyyyMMdd_HHmm');
    lbl = string(datestr(t,'dd mmm HH:MM'));
    plot(S.voltage,S.depth,'o-','LineWidth',1.5,'Color',co(k,:), ...
         'MarkerFaceColor',co(k,:),'DisplayName',lbl);
    fprintf('%s : %2d points, z = %.1f to %.1f cm, V = %.3f to %.3f\n', ...
        lbl, numel(S.depth), min(S.depth), max(S.depth), ...
        min(S.voltage), max(S.voltage));
end

set(gca,'YDir','reverse','FontSize',14,'FontName','times')
xlabel('Voltage (V)')
ylabel('z (cm)')
title('Density profiles after each run')
legend('Location','southwest')
grid on; box on

print(gcf,fullfile(here,'density_profiles_overlay.png'),'-dpng','-r150')
