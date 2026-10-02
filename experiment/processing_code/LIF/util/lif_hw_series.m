function hw = lif_hw_series(H, avg_s)
%LIF_HW_SERIES Hot-wire velocity series for the lif_surf_figure panel.
%
%   hw = LIF_HW_SERIES(H, avg_s)
%
% H     - struct loaded from a run_experiment.m hot-wire file
% avg_s - block-average length (s), e.g. 0.02 = one camera frame at 50 Hz.
%         The raw record is 4 kHz; plotting it whole is slow and unreadable.
%
% hw.t     - block centres, s since wind start (hwT_wind)
% hw.Y     - [nt x nseries] block means: U (hot-wire) and U_ref (reference
%            probe), whichever were logged
% hw.names - LaTeX legend names

    t = H.hwT_wind(:);
    Y = [];  names = {};
    if isfield(H, 'U') && ~isempty(H.U),         Y = [Y H.U(:)];     names{end+1} = '$U$ (hot-wire)'; end
    if isfield(H, 'U_ref') && ~isempty(H.U_ref), Y = [Y H.U_ref(:)]; names{end+1} = '$U_{\rm ref}$'; end
    if isempty(Y), error('lif_hw_series:empty', 'No U or U_ref in the hot-wire file.'); end

    dt = median(diff(t));
    n  = max(1, round(avg_s / dt));
    nb = floor(numel(t) / n);
    idx = 1 : nb*n;
    hw.t = mean(reshape(t(idx), n, nb), 1).';
    hw.Y = squeeze(mean(reshape(Y(idx, :), n, nb, []), 1));
    if isrow(hw.Y), hw.Y = hw.Y.'; end
    hw.names = names;
end
