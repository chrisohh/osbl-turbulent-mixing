% VideoPlotter_SurfaceDetection.m
% Generates a video of surface detection results (scaled surface image
% with filtered and raw surface overlays) for all image pairs.
clear
clc

LONG = '/media/surflab/LC_Working24/LC/FabMarcNovDec2014/data/Longitudinal/PIVdt10ms_IRlas1_8hz/';
DIRS = dir(LONG);
DIRS = DIRS(3:end);

for ii = 1%:length(DIRS)

exp_name = DIRS(ii).name;
num_of_digits = 3;
load_path = [LONG exp_name];
files = dir([load_path '/PIVRaw/PIV/*.mat']);
number_of_pair = length(files)/2;

% Video output setup
videoPath = [load_path '/Results_Surflab/'];
if ~exist(videoPath, 'dir')
    mkdir(videoPath);
end
videoFile = [videoPath exp_name '_SurfaceDetection.mp4'];
v = VideoWriter(videoFile, 'MPEG-4');
v.FrameRate = 8;
v.Quality = 95;
open(v);

% Figure setup
fig = figure('Position', [100 100 1200 800]);
colormap bone
cl = [1, 0.4, 0.4];

for image_pair_number = 0:number_of_pair-1
    imSurfa = FindSurfaceCapillary( ...
        [load_path '/PIVRaw/PIVSURF/' exp_name '_Pivsurf_' ...
         sprintf(['%0' num2str(num_of_digits) 'd'], image_pair_number) '_a.mat'], ...
        findMask = false);

    hold off
    imagesc(imSurfa.ImgScaledToPIVSmallCrop, [0, 300])
    colormap bone
    hold on
    plot(imSurfa.surfaceSurfImgScaled, '-', 'Color', cl, 'LineWidth', 1.5)
    plot(imSurfa.surface_raw, '-r', 'LineWidth', 1.5)
    daspect([1, 1, 1])
    title(sprintf('%s  Pair %d / %d', exp_name, image_pair_number, number_of_pair-1), ...
        'Interpreter', 'none')

    drawnow
    Frame = getframe(gcf);
    writeVideo(v, Frame);

    disp(['Frame ' num2str(image_pair_number) ' of ' num2str(number_of_pair-1)])
end

close(v);
close(fig);
disp(['Video saved: ' videoFile])

end
