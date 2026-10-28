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

    % The images are colour-mapped here and drawn as plain RGB: a
    % log-ColorScale image with AlphaData took ~50 s per video frame to
    % print with software OpenGL (the analysis itself < 1 s), and
    % decimating the pixels did not help.  The axes keep their colormap,
    % CLim and ColorScale, so the colour bars are unchanged.  NaN (air
    % above the surface in C) -> black, as the black axes behind it did.
    % Display only every show_ds-th pixel: imagesc places the CData by the
    % first and last axis values, so the coarser image lands on the same axes.
    k = 1;  if isfield(h.opts, 'show_ds'), k = h.opts.show_ds; end
    set(h.imR, 'CData', to_rgb(I(1:k:end, 1:k:end), h.axR));
    set(h.imC, 'CData', to_rgb(R.C_img(1:k:end, 1:k:end), h.axC));

    z99_line = z_surf(R.jc) + R.z99_y;
    % While the frame z99 is shallower than the run's floor the dye is still
    % a surface layer: no line (h.opts.z99_floor_cm from lif_z99_gate: Inf
    % once resolved, NaN while the floor is not known yet; [] = no gate)
    fl = [];  if isfield(h.opts, 'z99_floor_cm'), fl = h.opts.z99_floor_cm; end
    unresolved = ~isempty(fl) && ~(R.z99 < fl);
    if unresolved, z99_line(:) = NaN; end
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
    set(h.z9P, 'XData', xl);  set(h.bandP, 'XData', xl([1 2 2 1]));
    if isfield(h, 'pycP'), set(h.pycP, 'XData', xl([1 2 2 1])); end
    if unresolved
        set(h.z9P, 'YData', [NaN NaN]);  set(h.bandP, 'YData', [NaN NaN NaN NaN]);
        if isnan(fl), h.ttP.String = {'$z_{99}$ unresolved:', 'dye at the surface'};
        else,         h.ttP.String = {'$z_{99}$ unresolved:', sprintf('dye in top %.1f cm', -fl)}; end
    else
        set(h.z9P, 'YData', [R.z99 R.z99]);
        if ~isfield(h.opts, 'z99_band') || isempty(h.opts.z99_band), q = [NaN NaN];
        else, q = prctile(R.z99_y, h.opts.z99_band); end
        set(h.bandP, 'YData', q([1 1 2 2]));
        h.ttP.String = sprintf('$z_{99} = %.2f$ cm', R.z99);
    end
    h.tt.String  = ttl;
    if nargin >= 7
        for L = h.hwLine(:).', L.Value = t_frame; end   % wind and wave panels
    end
end

function rgb = to_rgb(A, ax)
    % Colour-map A with the axes' colormap, CLim and ColorScale; NaN -> black
    cm = colormap(ax);  cl = ax.CLim;  n = size(cm, 1);
    nanm = isnan(A);                                 % before max(), which drops NaN
    if strcmp(ax.ColorScale, 'log')
        A = log10(max(A, cl(1)));  cl = log10(cl);
    end
    idx = round(1 + (n - 1) * (A - cl(1)) / diff(cl));
    idx = min(max(idx, 1), n);  idx(nanm) = 1;
    rgb = reshape(cm(idx, :), [size(A) 3]);
    rgb(repmat(nanm, [1 1 3])) = 0;
end
