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
%   blank     false  missing frame: all-NaN data, title is opts.head only
%   sigma     []     struct(sx, sy, thr): ROI-sized slope uncertainties (info.sig_sx, sig_sy) and the mask threshold: draws the
%                    threshold contour on the slope panels (off by default)
%   ts        []     time series under the maps (cisg_wind_wave_series: wind, wave), as in the LIF videos
%   t_now     []     s since wind start of this frame: red cursor on the time series
%   ts_xlim   []     time-axis limits (s); [] = whole record
%   eta_hist  []     struct t, v: earlier frames' CISG eta at (0, 0) (mm, info.eta_centre_mm), drawn in orange on the
%                    wave-gauge panel with this frame's value.  CISG eta is relative (mean over the view removed), so
%                    it keeps waves shorter than the view only; the gauge is 1 m downstream: compare amplitudes, not phase
%
% Maps are in tank cm (x right, y up), eta (mm, relative) is integrated by
% cisg_slope_to_eta.m.  info: rms_sx, rms_sy, rms_eta, masked, lim_s, lim_eta, nan_fraction, ax_raw (raw-frame axes),
% eta_centre_mm (CISG eta at tank (0, 0), mean within 1 cm).

    if ~isfield(opts, 'head'), opts.head = ''; end
    if ~isfield(opts, 'bin') || isempty(opts.bin), opts.bin = 1; end
    if ~isfield(opts, 'lim_s'), opts.lim_s = []; end
    if ~isfield(opts, 'lim_eta'), opts.lim_eta = []; end
    if ~isfield(opts, 'q') || isempty(opts.q), opts.q = 16; end
    b = opts.bin;  q = max(2, round(opts.q / b));            % block in binned px (even)
    q = q + mod(q, 2);
    r = setup.roi;

    [eta, ei] = cisg_slope_to_eta(Sx, Sy, 10*setup.cm_per_pixel_x, 10*setup.cm_per_pixel_y, q);
    eta(ei.valid_frac < 0.5) = NaN;                                  % gaps are filled only for the solve: show them empty
    eta = eta - mean(eta, 'all', 'omitnan');
    [nrq, ncq] = size(eta);
    [Cq, Rq] = meshgrid(r.col_min + (0:ncq-1)*q + (q-1)/2, r.row_min + (0:nrq-1)*q + (q-1)/2);
    xy = surf_px2world(cal, [(Cq(:) - 0.5)*b + 0.5, (Rq(:) - 0.5)*b + 0.5]) / 10;   % cm
    y0 = 0;  if isfield(setup, 'y0_cm'), y0 = setup.y0_cm; end             % channel centre line = y 0 (cisg_channel_center)
    Xq = reshape(xy(:,1), size(Cq));  Yq = reshape(xy(:,2), size(Cq)) - y0;
    rsx = rms(Sx, 'all', 'omitnan');  rsy = rms(Sy, 'all', 'omitnan');  re = rms(eta(:), 'omitnan');
    lim_s = opts.lim_s;    if isempty(lim_s),   lim_s = 3 * max(rsx, rsy); end
    if ~isfinite(lim_s) || lim_s <= 0, lim_s = 0.3; end              % blank (all-NaN) frame before limits are known
    lim_e = opts.lim_eta;  if isempty(lim_e),   lim_e = 3 * re; end
    if ~isfinite(lim_e) || lim_e <= 0, lim_e = 5; end

    figure(fig);  clf(fig);
    series = {};                                              % wind / wave panels under the maps
    if isfield(opts, 'ts') && ~isempty(opts.ts)
        series = {opts.ts.wind, opts.ts.wave};  series = series(~cellfun(@isempty, series));
    end
    % tiled layout (compact spacing keeps titles and axis labels apart).  With series, a 4 x 3 grid:
    %   raw | Sx | Sy        (rows 1-2)
    %   eta | wind           (rows 3-4: eta map left; wind above wave gauge, same time axis, so the
    %       | wave gauge      red frame lines sit one above the other)
    if isempty(series)
        tl = tiledlayout(fig, 2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
        ax_m = arrayfun(@(k) nexttile(tl, k), 1:4);
    else
        tl = tiledlayout(fig, 4, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
        ax_m = [nexttile(tl, 1, [2 1]), nexttile(tl, 2, [2 1]), nexttile(tl, 3, [2 1]), nexttile(tl, 7, [2 1])];
    end
    sp = @(k) axes_(ax_m(k));
    ds = max(1, round(4 / b));                                % raw display stride
    ax_raw = sp(1);
    imagesc(1:ds:size(obs_image, 2), 1:ds:size(obs_image, 1), obs_image(1:ds:end, 1:ds:end, :));  axis image;  hold on;
    rectangle('Position', [r.col_min r.row_min r.col_max-r.col_min r.row_max-r.row_min], 'EdgeColor', 'r', 'LineWidth', 1.5);
    title('Raw frame (ROI in red)');
    maps = {Sx(q/2:q:nrq*q, q/2:q:ncq*q), Sy(q/2:q:nrq*q, q/2:q:ncq*q), eta};
    ttl = {'S_x', 'S_y', '\eta (mm, relative)'};
    lims = {lim_s, lim_s, lim_e};
    cmap = interp1([0 0.5 1], [0.2 0.3 0.8; 1 1 1; 0.8 0.2 0.2], linspace(0, 1, 256));
    for k = 1:3
        sp(k + 1);
        pcolor(Xq, Yq, double(maps{k}));  shading flat;  axis image;
        caxis(lims{k} * [-1 1]);  colormap(gca, cmap);  colorbar;
        set(gca, 'Color', [0.8 0.8 0.8]);  xlabel('x (cm)');  ylabel('y (cm)');  title(ttl{k});
    end
    if isfield(opts, 'sigma') && ~isempty(opts.sigma)                  % optional: contour of the slope uncertainty at the mask threshold (black)
        fn = {'sx', 'sy'};
        for k = 1:2
            sp(k + 1);  hold on;
            Z = opts.sigma.(fn{k})(q/2:q:nrq*q, q/2:q:ncq*q);
            contour(Xq, Yq, double(Z), opts.sigma.thr * [1 1], 'k-', 'LineWidth', 1);
        end
    end
    t_now = [];  if isfield(opts, 't_now'), t_now = opts.t_now; end
    xl_ts = [];                                               % one time axis for all series: wind start to the last sample
    if ~isempty(series), xl_ts = [0, max(cellfun(@(s) s.t(end), series))]; end
    if isfield(opts, 'ts_xlim') && ~isempty(opts.ts_xlim), xl_ts = opts.ts_xlim; end
    % CISG eta at the tank centre (x, y) = (0, 0): mean of the coarse cells within 1 cm, marked on the eta map
    eta_c = mean(eta(hypot(Xq, Yq) <= 1), 'omitnan');
    if isempty(eta_c), eta_c = NaN; end
    eh = [];  if isfield(opts, 'eta_hist'), eh = opts.eta_hist; end
    for k = 1:numel(series)                                   % wind above wave gauge, red line = this frame
        s = series{k};
        if numel(series) == 1, nexttile(tl, 8, [2 2]); else, nexttile(tl, 8 + 3 * (k - 1), [1 2]); end
        is_wave = isfield(s, 'symbol') && strcmp(s.symbol, '$\eta$');
        plot(s.t, s.Y, 'LineWidth', 0.7);  hold on;  grid on;  xlim(xl_ts);
        ttl_k = s.names{1};
        if is_wave                                            % CISG eta at (0, 0) so far (mm -> cm), relative to the frame mean
            tc = t_now;  vc = eta_c;
            if ~isempty(eh), tc = [eh.t(:); t_now];  vc = [eh.v(:); eta_c]; end
            if ~isempty(tc), plot(tc, vc / 10, '-', 'Color', [0.95 0.5 0.1 0.55], 'LineWidth', 1); end   % 4th value = alpha: gauge shows through
            plot(t_now, eta_c / 10, 'o', 'Color', [0.95 0.5 0.1], 'MarkerFaceColor', [0.95 0.5 0.1], 'MarkerSize', 6);   % this frame (visible with one frame too)
        end
        if ~isempty(t_now) && isfinite(t_now)
            xline(t_now, 'r-', 'LineWidth', 2);
            [~, i] = min(abs(s.t - t_now));
            if isfield(s, 'symbol')
                ttl_k = sprintf('%s = %.2f %s%s', s.symbol, s.Y(i, 1), s.unit, s.note);
            else                                                  % series cached before symbol / unit / note existed
                ttl_k = sprintf('%s = %.2f', s.names{1}, s.Y(i, 1));
            end
            if is_wave, ttl_k = sprintf('%s $|$ CISG $\\eta(0,0)$ = %.2f cm (orange)', ttl_k, eta_c / 10); end
        end
        ylabel(s.ylabel);  title(ttl_k, 'Interpreter', 'latex');
        if k == numel(series), xlabel('t since wind start (s)'); else, set(gca, 'XTickLabel', []); end
    end
    masked = 100 * mean(isnan(Sx), 'all');
    if isfield(opts, 'blank') && opts.blank
        title(tl, opts.head);                                      % missing frame: title only
    else
        title(tl, {opts.head, sprintf('rms S_x %.4f, S_y %.4f  |  masked %.1f%%  |  rms \\eta %.2f mm  |  %s map, L_{eff} %.1f cm', ...
                   rsx, rsy, masked, re, setup.colour_map, setup.L_eff_cm)});
    end
    info = struct('rms_sx', rsx, 'rms_sy', rsy, 'rms_eta', re, 'masked', masked, ...
                  'lim_s', lim_s, 'lim_eta', lim_e, 'nan_fraction', ei.nan_fraction, 'ax_raw', ax_raw, 'eta_centre_mm', eta_c);
end

function ax = axes_(ax)
    set(ancestor(ax, 'figure'), 'CurrentAxes', ax);          % make the tile current for the plot calls (keeps an invisible figure invisible)
end
