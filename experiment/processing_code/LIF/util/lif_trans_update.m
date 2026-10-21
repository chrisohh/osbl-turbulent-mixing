function lif_trans_update(h, I, R, y_cm, z_surf, ttl, t_frame)
%LIF_TRANS_UPDATE Fill a lif_trans_figure for one frame.
%
%   LIF_TRANS_UPDATE(h, I, R, y_cm, z_surf, ttl, t_frame)
%
% h       - from lif_trans_figure
% I       - [ny x nx] raw counts for the left panel
% R       - the lif_trans_conc output struct (C_img, zeta, C_bar, z99,
%           z99_y, jc, z_pen)
% z_surf  - [1 x nx] surface height per column (cm)
% ttl     - figure title (LaTeX)
% t_frame - s since wind start; moves the vertical line on the hot-wire
%           panel, if the figure has one

    Ishow = I;
    if h.opts.raw_log, Ishow = max(Ishow, h.opts.clims_log(1)); end   % log: clip <= 0
    set(h.imR, 'CData', Ishow);

    C = R.C_img;
    if h.opts.C_log, C = max(C, h.opts.C_clims_log(1)); end
    set(h.imC, 'CData', C, 'AlphaData', ~isnan(R.C_img));             % air blank

    z99_line = z_surf(R.jc) + R.z99_y;
    for p = {{h.sfR, h.z9R, h.pnR}, {h.sfC, h.z9C, h.pnC}}
        set(p{1}{1}, 'YData', z_surf);
        set(p{1}{2}, 'XData', y_cm(R.jc), 'YData', z99_line);
        set(p{1}{3}, 'YData', R.z_pen);
    end

    set(h.prP, 'XData', R.C_bar, 'YData', R.zeta);
    if isempty(h.opts.profile_xlim)
        xmax = max(R.C_bar, [], 'omitnan');
        if xmax > 0, xlim(h.axP, [-0.02 1.05] * xmax); end
    end
    xl = xlim(h.axP);
    set(h.z9P, 'XData', xl, 'YData', [R.z99 R.z99]);
    q = prctile(R.z99_y, [5 95]);
    set(h.bandP, 'XData', xl([1 2 2 1]), 'YData', q([1 1 2 2]));
    h.ttP.String = sprintf('$z_{99} = %.2f$ cm', R.z99);
    h.tt.String  = ttl;
    if nargin >= 7 && ~isempty(h.hwLine), h.hwLine.Value = t_frame; end
end
