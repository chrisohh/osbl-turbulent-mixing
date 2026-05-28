function roi = pick_roi_interactive(img, prior_roi, scale_x, scale_y, title_text)
% PICK_ROI_INTERACTIVE  Two-click ROI picker. Returns roi in pixel indices.

[ny, nx, ~] = size(img);
x_axis = ((1:nx) - nx/2) * scale_x;
y_axis = ((1:ny) - ny/2) * scale_y;

figure('Name', title_text, 'Position', [100 100 1200 800]);
imagesc(x_axis, y_axis, img);
axis image; set(gca, 'YDir', 'normal');
hold on;

if ~isempty(prior_roi)
    x_box = ([prior_roi.col_min, prior_roi.col_max] - nx/2) * scale_x;
    y_box = ([prior_roi.row_min, prior_roi.row_max] - ny/2) * scale_y;
    rectangle('Position', [x_box(1), y_box(1), diff(x_box), diff(y_box)], ...
              'EdgeColor', 'r', 'LineWidth', 1.5);
    title({title_text, 'Red = current ROI. Click two opposite corners of the new ROI.'});
else
    title({title_text, 'Click two opposite corners of the ROI.'});
end

fprintf('Click two opposite corners of the ROI...\n');
[cx, cy] = ginput(2);

col_pix = cx / scale_x + nx/2;
row_pix = cy / scale_y + ny/2;

col_min = max(1,  ceil(min(col_pix)));
col_max = min(nx, floor(max(col_pix)));
row_min = max(1,  ceil(min(row_pix)));
row_max = min(ny, floor(max(row_pix)));

roi = struct('row_min', row_min, 'row_max', row_max, ...
             'col_min', col_min, 'col_max', col_max);

fprintf('ROI: rows [%d %d], cols [%d %d]\n', row_min, row_max, col_min, col_max);
end
