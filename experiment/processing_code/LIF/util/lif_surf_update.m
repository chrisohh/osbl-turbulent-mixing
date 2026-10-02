function lif_surf_update(h, I, theta_img, x_cm, z_front, z_surf, z99_x, jc, zeta, theta_bar, z99, ttl, t_frame)
%LIF_SURF_UPDATE Fill a lif_surf_figure for one frame.
%
%   LIF_SURF_UPDATE(h, I, theta_img, x_cm, z_front, z_surf, z99_x, jc, ...
%                   zeta, theta_bar, z99, ttl, t_frame)
%
% I / theta_img [ny x nx]; z_front, z_surf [1 x nx]; z99_x, jc from
% lif_z99; zeta, theta_bar the all-x profile; z99 the all-x value (cm,
% negative); ttl the figure title (LaTeX); t_frame (s since wind start)
% moves the vertical line on the hot-wire panel, if the figure has one.

    set(h.imR, 'CData', I, 'AlphaData', ~isnan(I));   % I = T_img or raw counts (h.opts.left)
    T = theta_img;
    if h.opts.theta_log, T = max(T, h.opts.theta_clims_log(1)); end   % log: clip <= 0
    set(h.imT, 'CData', T, 'AlphaData', ~isnan(theta_img));            % air blank

    z99_line = z_surf(jc) + z99_x;
    for p = {{h.sfR, h.z9R}, {h.sfT, h.z9T}}
        set(p{1}{1}, 'YData', z_surf);
        set(p{1}{2}, 'XData', x_cm(jc), 'YData', z99_line);
    end

    set(h.prP, 'XData', theta_bar, 'YData', zeta);
    set(h.z9P, 'YData', [z99 z99]);
    q = prctile(z99_x, [5 95]);
    set(h.bandP, 'YData', q([1 1 2 2]));
    h.ttP.String = sprintf('$z_{99} = %.2f$ cm', z99);
    h.tt.String  = ttl;
    if nargin >= 13 && ~isempty(h.hwLine), h.hwLine.Value = t_frame; end
end
