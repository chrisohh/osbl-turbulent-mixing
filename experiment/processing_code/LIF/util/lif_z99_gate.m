function [fl, G] = lif_z99_gate(G, t, z99)
%LIF_Z99_GATE Run-adaptive "z99 resolved" gate for the side-view videos.
%
%   G = LIF_Z99_GATE(opts)            initialise
%   [fl, G] = LIF_Z99_GATE(G, t, z99) one frame (t s since wind start, frame z99 cm)
%
% fl goes into h.opts.z99_floor_cm for lif_surf_update: the z99(x) line is
% drawn when z99 < fl.  Inf once the run has resolved (latched: the dye does
% not un-deepen, and a line that blinks on and off is worse than either);
% before that, the current threshold (NaN while the floor is not known yet).
%
% A fixed floor (-2.3 cm, from run 95) hid the start of the deepening in the
% stratified runs, whose surface layer reads much thinner (118: z99 ~ -0.3
% cm before the onset, so it was already at -2 cm when the line appeared).
% The floor is instead each run's own plateau, as in lif_depth_series:
%   floor = median z99 over opts.floor_win (frames seen so far)
%   thr   = min(floor_fac * floor, floor - max(min_margin_cm, k_sigma * sigma))
%           sigma = 1.4826 * MAD of z99 in the window (frame-to-frame noise)
%   resolved once the median of z99 over the last med_s seconds < thr.
%
% opts (defaults): floor_win [20 40], floor_fac 1.25, k_sigma 3,
% min_margin_cm 0.25, med_s 0.6, min_n 5 (frames in the window before the
% floor is used), fixed [] (a number = the old fixed floor, no adaptation).

    if nargin < 2
        if nargin == 0, G = struct(); end
        def = struct('floor_win', [20 40], 'floor_fac', 1.25, 'k_sigma', 3, ...
                     'min_margin_cm', 0.25, 'med_s', 0.6, 'min_n', 5, 'fixed', []);
        f = fieldnames(def);
        for k = 1:numel(f), if ~isfield(G, f{k}), G.(f{k}) = def.(f{k}); end, end
        G.t = [];  G.z = [];  G.latched = false;  G.t_on = NaN;  G.thr = NaN;  G.floor = NaN;
        fl = G;  return
    end

    if ~isempty(G.fixed), fl = G.fixed; return; end
    G.t(end+1) = t;  G.z(end+1) = z99;
    if G.latched, fl = Inf; return; end

    w = G.t >= G.floor_win(1) & G.t <= G.floor_win(2) & ~isnan(G.z);
    if nnz(w) < G.min_n, fl = NaN; return; end
    G.floor = median(G.z(w));
    sig     = 1.4826 * median(abs(G.z(w) - G.floor));
    G.thr   = min(G.floor_fac * G.floor, G.floor - max(G.min_margin_cm, G.k_sigma * sig));

    r = G.t > t - G.med_s & ~isnan(G.z);
    if t > G.floor_win(1) && any(r) && median(G.z(r)) < G.thr
        G.latched = true;  G.t_on = t;  fl = Inf;
    else
        fl = G.thr;
    end
end
