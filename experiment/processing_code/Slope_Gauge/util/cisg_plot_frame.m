function info = cisg_plot_frame(fig, obs_image, Sx, Sy, setup, cal, opts)
%CISG_PLOT_FRAME 2x2 slope figure of one CISG frame: raw | Sx | Sy | eta.
%
%   info = cisg_plot_frame(fig, obs_image, Sx, Sy, setup, cal, opts)
%
% obs_image - raw RGB frame (full resolution, or binned by opts.bin)
% Sx, Sy    - slopes inside the ROI (setup.roi, in the pixels of obs_image)
% setup     - needs roi, cm_per_pixel_x/y (of obs_image's pixels), colour_map, L_eff_cm
% opts fields (defaults):
%   head      ''     first part of the sgtitle (run, frame, time)
%   bin       1      obs_image is binned by this factor (pixel coordinates are mapped back for the cal)
%   lim_s     []     symmetric slope colour limit; [] = 3 * max(rms Sx, rms Sy) of this frame
%   lim_eta   []     symmetric eta limit (mm);      [] = 3 * rms eta of this frame
%   q         16     display / integration block in full-resolution px
%
% Maps are in tank cm (x right, y up), eta (mm, relative) is integrated by
% cisg_slope_to_eta.m.  info: rms_sx, rms_sy, rms_eta, masked, lim_s, lim_eta, nan_fraction.

    if ~isfield(opts, 'head'), opts.head = ''; end
    if ~isfield(opts, 'bin') || isempty(opts.bin), opts.bin = 1; end
    if ~isfield(opts, 'lim_s'), opts.lim_s = []; end
    if ~isfield(opts, 'lim_eta'), opts.lim_eta = []; end
    if ~isfield(opts, 'q') || isempty(opts.q), opts.q = 16; end
    b = opts.bin;  q = max(2, round(opts.q / b));            % block in binned px (even)
    q = q + mod(q, 2);
    r = setup.roi;

    [eta, ei] = cisg_slope_to_eta(Sx, Sy, 10*setup.cm_per_pixel_x, 10*setup.cm_per_pixel_y, q);
    [nrq, ncq] = size(eta);
    [Cq, Rq] = meshgrid(r.col_min + (0:ncq-1)*q + (q-1)/2, r.row_min + (0:nrq-1)*q + (q-1)/2);
    xy = surf_px2world(cal, [(Cq(:) - 0.5)*b + 0.5, (Rq(:) - 0.5)*b + 0.5]) / 10;   % cm
    Xq = reshape(xy(:,1), size(Cq));  Yq = reshape(xy(:,2), size(Cq));
    rsx = rms(Sx, 'all', 'omitnan');  rsy = rms(Sy, 'all', 'omitnan');  re = rms(eta(:));
    lim_s = opts.lim_s;    if isempty(lim_s),   lim_s = 3 * max(rsx, rsy); end
    lim_e = opts.lim_eta;  if isempty(lim_e),   lim_e = 3 * re; end

    figure(fig);  clf(fig);
    ds = max(1, round(4 / b));                                % raw display stride
    subplot(2,2,1);
    imagesc(1:ds:size(obs_image, 2), 1:ds:size(obs_image, 1), obs_image(1:ds:end, 1:ds:end, :));  axis image;  hold on;
    rectangle('Position', [r.col_min r.row_min r.col_max-r.col_min r.row_max-r.row_min], 'EdgeColor', 'r', 'LineWidth', 1.5);
    title('Raw frame (ROI in red)');
    maps = {Sx(q/2:q:nrq*q, q/2:q:ncq*q), Sy(q/2:q:nrq*q, q/2:q:ncq*q), eta};
    ttl = {'S_x', 'S_y', '\eta (mm, relative)'};
    lims = {lim_s, lim_s, lim_e};
    cmap = interp1([0 0.5 1], [0.2 0.3 0.8; 1 1 1; 0.8 0.2 0.2], linspace(0, 1, 256));
    for k = 1:3
        subplot(2,2,k+1);
        pcolor(Xq, Yq, double(maps{k}));  shading flat;  axis image;
        caxis(lims{k} * [-1 1]);  colormap(gca, cmap);  colorbar;
        set(gca, 'Color', [0.8 0.8 0.8]);  xlabel('x (cm)');  ylabel('y (cm)');  title(ttl{k});
    end
    masked = 100 * mean(isnan(Sx), 'all');
    sgtitle(sprintf('%s | rms S_x %.4f, S_y %.4f | masked %.1f%% | rms \\eta %.2f mm | %s map, L_{eff} %.1f cm', ...
        opts.head, rsx, rsy, masked, re, setup.colour_map, setup.L_eff_cm));
    info = struct('rms_sx', rsx, 'rms_sy', rsy, 'rms_eta', re, 'masked', masked, ...
                  'lim_s', lim_s, 'lim_eta', lim_e, 'nan_fraction', ei.nan_fraction);
end
