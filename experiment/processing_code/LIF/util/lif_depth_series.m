function D = lif_depth_series(ts_file, opts)
%LIF_DEPTH_SERIES Load a surf_theta_video time series and condition it.
%
%   D = LIF_DEPTH_SERIES(ts_file, opts)
%
% Shared by plot_depth_timeseries.m (one run) and compare_depth_runs.m.
%   - time filter: running median over opts.med_win_s; frames further than
%     opts.jump_cm from it are flagged (raw values are kept)
%   - side-view floor: median of the filtered series over opts.floor_win
%     (s since wind start), while the dye is still a surface layer; frames
%     deeper than opts.floor_fac * floor are 'resolved'
%   - deep: frames where dye reaches the bottom of the view over more than
%     half the width (z99 undefined there)
%
% opts (defaults): med_win_s 1, jump_cm 1, floor_win [20 40], floor_fac 1.25,
%   measures {'z99_all','z90_all','zrms_all'}
%
% D fields: file, t (s since wind start), run_num, R (run config, if saved),
%   and per measure m (e.g. D.z99_all): .raw, .filt, .flag, .floor,
%   .resolved, .t_on (first resolved time after floor_win); D.deep (logical),
%   D.z99_x_all (if saved).

    def = struct('med_win_s', 1, 'jump_cm', 1, 'floor_win', [20 40], 'floor_fac', 1.25, ...
                 'measures', {{'z99_all', 'z90_all', 'zrms_all'}});
    if nargin < 2, opts = struct(); end
    for f = fieldnames(def).'
        if ~isfield(opts, f{1}), opts.(f{1}) = def.(f{1}); end
    end

    S = load(ts_file);
    D.file = ts_file;
    D.t = S.t_wind(:);
    D.run_num = [];  D.R = [];
    if isfield(S, 'run_num'), D.run_num = S.run_num; end
    if isfield(S, 'R'), D.R = S.R; end
    if isfield(S, 'z99_x_all'), D.z99_x_all = S.z99_x_all; end
    D.deep = false(size(D.t));
    if isfield(S, 'deep_all'), D.deep = S.deep_all(:) > 0.5; end

    dt = median(diff(D.t));
    kw = max(1, 2*round(opts.med_win_s / dt / 2) + 1);    % odd window, frames
    in_win = D.t >= opts.floor_win(1) & D.t <= opts.floor_win(2);
    D.opts = opts;  D.kw = kw;

    for m = 1:numel(opts.measures)
        nm = opts.measures{m};
        if ~isfield(S, nm), continue; end
        M.raw  = S.(nm)(:);
        M.filt = movmedian(M.raw, kw, 'omitnan');
        M.flag = abs(M.raw - M.filt) > opts.jump_cm;
        M.floor = median(M.filt(in_win), 'omitnan');
        M.resolved = M.filt < opts.floor_fac * M.floor;
        M.t_on = D.t(find(M.resolved & D.t > opts.floor_win(2), 1, 'first'));
        if isempty(M.t_on), M.t_on = NaN; end
        D.(nm) = M;
    end
end
