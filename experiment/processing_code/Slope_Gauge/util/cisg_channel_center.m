function cc = cisg_channel_center(ref_image, cal, setup, r)
%CISG_CHANNEL_CENTER Channel centre line (y = 0) from the dark wall edges in the reference image.
%
%   cc = cisg_channel_center(ref_image, cal, setup, r)
%
% The walls block the view of the bottom pattern, which leaves dark bands at the top and bottom of the
% reference image.  For each of ~15 columns across the target the row profile of the mean RGB is smoothed and
% the steepest rise (top band -> bright) in the first 30% of the rows and the steepest fall (bright -> bottom
% band) in the last 30% are taken as the edges; the median over columns handles a wavy edge, and an edge counts
% only if the step is >= 25% of the bright level and found in >= 60% of the columns within +-40 px.
%
% Geometry: the edges lie at the pattern depth, seen through water, so their surface-plane coordinates are
% shrunk by 1/r about the optical axis.  With y_a the tank y of the optical axis (cal.center) and y_t, y_b the
% apparent surface-plane y of the top / bottom edge,  y0 = y_a + r * ((y_t + y_b)/2 - y_a)  (cm).
% If either edge is missing, y0 = NaN (the caller falls back to the plate centre and warns).
%
% setup: needs target_cols, pattern_scale (= r).   r: ignored name kept for the signature (use setup.pattern_scale).
% cc fields: top_row, bottom_row (full-res rows, median), found_top, found_bottom, strength, y_top_cm, y_bottom_cm,
%   y_axis_cm, y0_cm, mid_apparent_cm, plate_centre_y_cm (tank y of the image centre row for comparison)

    if nargin < 4, r = setup.pattern_scale; end %#ok<NASGU>
    rr = setup.pattern_scale;
    S = mean(double(ref_image), 3);
    [h, w] = size(S);
    tc = round(setup.target_cols);
    cols = round(linspace(tc(1) + 0.2 * diff(tc), tc(2) - 0.2 * diff(tc), 15));
    top = nan(1, numel(cols));  bot = top;  st_top = top;  st_bot = top;
    for j = 1:numel(cols)
        p = mean(S(:, max(1, cols(j) - 15):min(w, cols(j) + 15)), 2);
        p = imgaussfilt(p, 8);
        dp = gradient(p);
        n1 = round(0.3 * h);
        [m, e] = max(dp(1:n1));                                % rising edge, top band -> bright
        lvl = median(p(e + 20:min(h, e + 120)));
        if p(max(1, e - 40)) < lvl * 0.75, top(j) = e;  st_top(j) = (lvl - p(max(1, e - 40))) / lvl; end
        [m, e2] = min(dp(round(0.7 * h):h));  e2 = e2 + round(0.7 * h) - 1;   %#ok<ASGLU>
        lvl2 = median(p(max(1, e2 - 120):e2 - 20));
        if p(min(h, e2 + 40)) < lvl2 * 0.75, bot(j) = e2;  st_bot(j) = (lvl2 - p(min(h, e2 + 40))) / lvl2; end
    end
    cc = struct('top_row', NaN, 'bottom_row', NaN, 'found_top', false, 'found_bottom', false, 'strength', [NaN NaN], ...
                'y_top_cm', NaN, 'y_bottom_cm', NaN, 'y_axis_cm', NaN, 'y0_cm', NaN, 'mid_apparent_cm', NaN);
    ok = @(v) sum(abs(v - median(v, 'omitnan')) <= 40) >= 0.6 * numel(v) && sum(~isnan(v)) >= 0.6 * numel(v);
    cc.found_top = ok(top);  cc.found_bottom = ok(bot);
    if cc.found_top, cc.top_row = median(top, 'omitnan');  cc.strength(1) = median(st_top, 'omitnan'); end
    if cc.found_bottom, cc.bottom_row = median(bot, 'omitnan');  cc.strength(2) = median(st_bot, 'omitnan'); end
    uc = mean(cols);
    xy_axis = surf_px2world(cal, cal.center);
    cc.y_axis_cm = xy_axis(2) / 10;
    cc.plate_centre_y_cm = 0;
    if cc.found_top && cc.found_bottom
        xt = surf_px2world(cal, [uc cc.top_row]);  xb = surf_px2world(cal, [uc cc.bottom_row]);
        cc.y_top_cm = xt(2) / 10;  cc.y_bottom_cm = xb(2) / 10;
        cc.mid_apparent_cm = (cc.y_top_cm + cc.y_bottom_cm) / 2;
        cc.y0_cm = cc.y_axis_cm + rr * (cc.mid_apparent_cm - cc.y_axis_cm);
    end
end
