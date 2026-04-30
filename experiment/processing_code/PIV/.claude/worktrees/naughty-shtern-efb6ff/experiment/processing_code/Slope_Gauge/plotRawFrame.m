coreview_num=22;
obs_folderName = sprintf('Y:\\HLAB_2026\\SlopeGauge\\CoreView_%d\\Flare 12M125 CCL C1112A00\\', coreview_num);
n=1;obs_fileName = sprintf("CoreView_%d_Flare 12M125 CCL C1112A00_%04d.raw", coreview_num, n);[obs_image, ~] = cisg_load_coreview(strcat(obs_folderName, obs_fileName));imshow(obs_image);
[img_height, img_width] = size(obs_image, [1 2]);

load("Y:\HLAB_2026\SlopeGauge\Processed\CISG_Calibration_CoreView_6.mat")
setup=cal;
nx = img_width;
ny = img_height;
x_cm = ((1:nx) - nx/2) * setup.cm_per_pixel_x;
y_cm = ((1:ny) - ny/2) * setup.cm_per_pixel_y;

%%
n=2271;
obs_fileName = sprintf("CoreView_%d_Flare 12M125 CCL C1112A00_%04d.raw", coreview_num, n);
[obs_image, ~] = cisg_load_coreview(strcat(obs_folderName, obs_fileName));
imshow(obs_image);
[img_height, img_width] = size(obs_image, [1 2]);

imagesc(x_cm,y_cm,obs_image)
axis equal tight
fs = 50;  % Hz
dt = 1/fs;  % 0.02 seconds between frames
time=(n - 1) * dt;
title(sprintf('$t = %.2f$ s',time),'Interpreter','latex','fontsize',16)
set(gca,'fontsize',16,'fontname','times')
% xlim([-14,14])
set(gcf,'color','white')
xlabel('$x$ (cm)','Interpreter','latex');ylabel('$y$ (cm)','Interpreter','latex')