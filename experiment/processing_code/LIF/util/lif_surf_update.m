function lif_surf_update(h, I, theta_img, x_cm, z_front, z_surf, z99_x, jc, zeta, theta_bar, z99, ttl, t_frame, st)
%LIF_SURF_UPDATE Fill a lif_surf_figure for one frame.
%
%   LIF_SURF_UPDATE(h, I, theta_img, x_cm, z_front, z_surf, z99_x, jc, ...
%                   zeta, theta_bar, z99, ttl, t_frame)
%
% I / theta_img [ny x nx]; z_front, z_surf [1 x nx]; z99_x, jc from
% lif_z99; zeta, theta_bar the all-x profile; z99 the all-x value (cm,
% negative); ttl the figure title (LaTeX); t_frame (s since wind start)
% moves the vertical line on the hot-wire panel, if the figure has one;
% st (lif_z99's 8th output) draws the grey surface band and unlit columns.
%
% Panels the figure does not have (opts.panels 'theta', opts.show_z99
% false) come back as empty handles from lif_surf_figure and are skipped,
% so the same call fills either layout.

    if ~isempty(h.imR)
        set(h.imR, 'CData', I, 'AlphaData', ~isnan(I));   % T_img or raw counts (h.opts.left)
    end
    T = theta_img;
    if h.opts.theta_log, T = max(T, h.opts.theta_clims_log(1)); end   % log: clip <= 0
    set(h.imT, 'CData', T, 'AlphaData', ~isnan(theta_img));            % air blank

    z99_line = z_surf(jc) + z99_x;
    % While the frame z99 is shallower than the side-view floor the dye is
    % still a surface layer and z99(x) is noise: draw no line (opts
    % z99_floor_cm; CoreView 95: -2.3 cm, resolved from t ~ 42 s)
    fl = h.opts.z99_floor_cm;
    unresolved = ~isempty(fl) && ~(z99 < fl);
    if unresolved, z99_line(:) = NaN; end
    for p = {{h.sfR, h.z9R}, {h.sfT, h.z9T}}
        if ~isempty(p{1}{1}), set(p{1}{1}, 'YData', z_surf); end
        if ~isempty(p{1}{2}), set(p{1}{2}, 'XData', x_cm(jc), 'YData', z99_line); end
    end

    if ~isempty(h.prP)
        set(h.prP, 'XData', theta_bar, 'YData', zeta);
        if ~isempty(h.z9P)
            if unresolved
                set(h.z9P, 'YData', [NaN NaN]);
                set(h.bandP, 'YData', [NaN NaN NaN NaN]);
                h.ttP.String = sprintf('$z_{99}$ unresolved (above %.1f cm)', fl);
            else
                set(h.z9P, 'YData', [z99 z99]);
                q = prctile(z99_x, [5 95]);
                set(h.bandP, 'YData', q([1 1 2 2]));
                h.ttP.String = sprintf('$z_{99} = %.2f$ cm', z99);
            end
        end
    end
    h.tt.String  = ttl;
    if nargin >= 13 && ~isempty(h.hwLine), h.hwLine.Value = t_frame; end

    % Grey layers from lif_z99's stats: uniform surface band (band_cm under
    % the local surface) and the unlit columns
    if nargin >= 14 && ~isempty(st) && isfield(st, 'band_cm')
        lit = st.lit(:).';
        xl  = x_cm(lit);  zt = z_surf(lit);
        for hb = [h.bdR h.bdT]
            set(hb, 'XData', [xl fliplr(xl)], 'YData', [zt fliplr(zt - st.band_cm)]);
        end
        d  = diff([0 ~lit 0]);  s = find(d == 1);  e = find(d == -1) - 1;   % unlit runs
        dx = median(diff(x_cm)) / 2;
        for hu = [h.unR h.unT]
            yl = get(get(hu, 'Parent'), 'YLim');
            V = zeros(4*numel(s), 2);  F = zeros(numel(s), 4);
            for k = 1:numel(s)
                x0 = x_cm(s(k)) - dx;  x1 = x_cm(e(k)) + dx;
                V(4*k-3:4*k, :) = [x0 yl(1); x1 yl(1); x1 yl(2); x0 yl(2)];
                F(k, :) = 4*k-3 : 4*k;
            end
            if isempty(s), set(hu, 'Vertices', [NaN NaN], 'Faces', 1);
            else,          set(hu, 'Vertices', V, 'Faces', F); end
        end
    end
end
