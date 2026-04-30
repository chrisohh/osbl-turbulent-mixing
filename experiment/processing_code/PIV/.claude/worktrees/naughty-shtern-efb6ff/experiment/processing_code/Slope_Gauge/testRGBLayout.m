%% Quick test of channel reordering
% Test if GRB→RGB swap fixes your colors

filename = strcat(ref_folderName,ref_fileName);  % <-- CHANGE THIS

% Load with corrected channel order
[img, meta] = cisg_load_coreview(filename);

% Display
figure('Position', [100 100 1200 800]);

% Full image
subplot(1,4,1);
imshow(img);
title('Full Image (GBR→RGB corrected)');

% Individual channels
subplot(1,4,2);
imagesc(img(:,:,1)); axis equal tight; colorbar;
title('Red Channel \n(should decrease left→right)');
colormap(gca, 'gray');

subplot(1,4,3);
imagesc(img(:,:,2)); axis equal tight; colorbar;
title('Green Channel \n(should increase bottom→top)');
colormap(gca, 'gray');

subplot(1,4,4);
imagesc(img(:,:,3)); axis equal tight; colorbar;
title('Blue Channel \n(should be constant ~128)');
colormap(gca, 'gray');
