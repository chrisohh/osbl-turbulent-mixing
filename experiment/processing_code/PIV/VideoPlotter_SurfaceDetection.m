% VideoPlotter_SurfaceDetection.m
% Generates a video of surface detection results (scaled surface image
% with filtered and raw surface overlays) for all image pairs.
LONG = 'D:\DelawareDataBackup\Longitudinal\PIV\';
DIRS = dir(LONG);
DIRS = DIRS(3:end);

Fs_PIV       = 7.2;        % Hz — frame rate between consecutive pairs
delay        = 73 + 5;     % s  — imaging starts this many s after trigger
t_wind_onset = 40;         % s  — trigger to wind onset
DX           = 1/17697.69; % m per pixel

for ii = 4%:length(DIRS)

exp_name = DIRS(ii).name;
num_of_digits = 3;
load_path = [LONG exp_name];
files = dir([load_path '/PIVRaw/PIV/*.mat']);
number_of_pair = length(files)/2;

% Video output setup
videoPath = ['D:\DelawareDataResult\'];
if ~exist(videoPath, 'dir')
    mkdir(videoPath);
end
videoFile = [videoPath exp_name '_SurfaceDetection.mp4'];
v = VideoWriter(videoFile, 'MPEG-4');
v.FrameRate = Fs_PIV;
v.Quality = 100;
open(v);

% Figure setup
fig = figure('Position', [100 100 1000 800]);
colormap bone
cl = [1, 0.4, 0.4];

for image_pair_number = 0:number_of_pair-1
    imSurfa = FindSurfaceCapillary( ...
        [load_path '/PIVRaw/PIVSURF/' exp_name '_Pivsurf_' ...
         sprintf(['%0' num2str(num_of_digits) 'd'], image_pair_number) '_a.mat'], ...
        findMask = false);

    [h_s, w_s] = size(imSurfa.ImgScaledToPIVSmallCrop);
    t_since_wind = image_pair_number / Fs_PIV + delay - t_wind_onset;

    hold off
    imagesc((1:w_s)*DX*1e3, (1:h_s)*DX*1e3, imSurfa.ImgScaledToPIVSmallCrop, [0, 300])
    colormap bone
    hold on
    plot((1:length(imSurfa.surfaceSurfImgScaled))*DX*1e3, ...
         imSurfa.surfaceSurfImgScaled*DX*1e3, '-', 'Color', cl, 'LineWidth', 1)
    plot((1:length(imSurfa.surface_raw))*DX*1e3, ...
         imSurfa.surface_raw*DX*1e3, '-r', 'LineWidth', 1)
    daspect([1, 1, 1])
    xlabel('x (mm)'); ylabel('z (mm)')
    title(sprintf('%s  Frame %d_a,  t = %.1f s from wind start', ...
        exp_name, image_pair_number, t_since_wind), 'Interpreter', 'none')
    set(gca,'fontsize',15);axis equal;axis tight
    drawnow
    Frame = getframe(gcf);
    % Lock frame size to first frame (ImgScaledToPIVSmallCrop size varies
    % slightly between frames due to surface-dependent cropping)
    if image_pair_number == 0
        frame_size = size(Frame.cdata(:,:,1));
    else
        Frame.cdata = imresize(Frame.cdata, frame_size);
    end
    writeVideo(v, Frame);

    disp(['Frame ' num2str(image_pair_number) ' of ' num2str(number_of_pair-1)])
end

close(v);
disp(['Video saved: ' videoFile])

end
