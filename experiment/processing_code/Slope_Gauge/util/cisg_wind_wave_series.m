function ts = cisg_wind_wave_series(coreview_num, avg_s)
%CISG_WIND_WAVE_SERIES Wind and wave-gauge time series of a run for the CISG video panels.
%
%   ts = cisg_wind_wave_series(coreview_num)
%   ts = cisg_wind_wave_series(coreview_num, avg_s)      % block average (s), default 0.02 (one frame)
%
% Files from common/run_registry.csv; time is s since wind (fan) start, the clock of
% cisg_wind_delay, so a frame sits at R.t_wind(n).
%   wind - hot-wire / reference probe velocity (lif_hw_series) when the run has one; else the registry
%          wind_standin (another run's hot-wire, labelled "stand-in: core N"; run 101 -> core 95); else the
%          fan command converted to the nominal speed with Run_exp/fan_transfer.mat (labelled so), or
%          the fan command voltage itself where fan_transfer.mat is not on this computer
%   wave - wave-gauge eta (cm) of wave_file (or of the wind file when that is the wave gauge)
% Each is a struct t, Y, names (LaTeX), ylabel, or [] when the run has none.
% The result is cached in Slope_Gauge/calibration/windwave_CoreView_<N>_avg<avg_s>.mat, which is
% loaded when the record files are not on this computer.

    if nargin < 2 || isempty(avg_s), avg_s = 0.02; end
    sg = fileparts(fileparts(mfilename('fullpath')));
    pc = fileparts(sg);
    addpath(fullfile(pc, 'common'), fullfile(pc, 'LIF', 'util'));      % run_registry, lif_hw_series
    cache = fullfile(sg, 'calibration', sprintf('windwave_CoreView_%d_avg%g.mat', coreview_num, avg_s));
    ts = struct('wind', [], 'wave', [], 'run', coreview_num, 'avg_s', avg_s);
    try
        R = run_registry(coreview_num);
    catch ME
        warning('cisg_wind_wave_series:registry', 'Run %d not in run_registry (%s)', coreview_num, ME.message);
        return;
    end
    wave_path = R.wave_path;
    if isempty(wave_path) && contains(R.wind_kind, 'wave'), wave_path = R.wind_path; end
    have = @(p) ~isempty(p) && isfile(p);
    if ~have(R.wind_path) && ~have(wave_path)
        if isfile(cache)
            S = load(cache, 'ts');  ts = S.ts;
            fprintf('wind / wave records not on this computer: loaded %s\n', cache);
        else
            warning('cisg_wind_wave_series:none', 'Run %d: no wind or wave record found (and no cache)', coreview_num);
        end
        return;
    end

    % wind
    if have(R.wind_path)
        H = load(R.wind_path);
        hw = [];
        if ~strcmp(R.wind_kind, 'wave')
            try hw = lif_hw_series(H, avg_s); catch, end
        end
        if ~isempty(hw) && ~strcmp(hw.kind, 'eta')
            ts.wind = struct('t', hw.t, 'Y', hw.Y, 'names', {hw.names}, 'ylabel', hw.ylabel, 'symbol', '$U$', 'unit', 'm/s', 'note', ' (hot-wire)');
        elseif have(R.standin_path)
            % no velocity record: another run's hot-wire (registry wind_standin, same fan ramp), labelled with its run
            hw = lif_hw_series(load(R.standin_path), avg_s);
            [~, nm] = fileparts(R.standin_path);
            T = run_registry();  rs = T.run(strcmp(T.wind_file, R.wind_standin));
            lab = strrep(nm, '_', '\_');  if ~isempty(rs), lab = sprintf('core %d', rs(1)); end   % file name only when the run is unknown
            ts.wind = struct('t', hw.t, 'Y', hw.Y, 'ylabel', hw.ylabel, ...
                             'names', {cellfun(@(s) sprintf('%s (stand-in: %s)', s, lab), hw.names, 'UniformOutput', false)}, ...
                             'symbol', '$U$', 'unit', 'm/s', 'note', sprintf(' (%s hot-wire, stand-in)', lab), 'source', nm);
        elseif all(isfield(H, {'sentT', 'sentV', 'fanStartElapsed'}))
            % no velocity record: the fan command through the fan transfer curve (nominal speed)
            t = H.sentT(:) - H.fanStartElapsed;  V = H.sentV(:);
            tf = {fullfile(pc, 'Run_exp', 'fan_transfer.mat'), ...                  % this checkout, or the lab computer's
                  'D:\Chris\osbl-turbulent-mixing\experiment\processing_code\Run_exp\fan_transfer.mat'};
            tf = tf(cellfun(@isfile, tf));
            if ~isempty(tf)
                tf = tf{1};
                F = load(tf, 'fanTF');
                U = interp1(F.fanTF.V, F.fanTF.U, V, 'pchip', NaN);
                U(V < F.fanTF.V(1)) = 0;                                % below the transfer range: fan off
                ts.wind = struct('t', t, 'Y', U, 'names', {{'$U$ from fan command (fan\_transfer)'}}, 'ylabel', 'U (m/s)', ...
                                 'symbol', '$U$', 'unit', 'm/s', 'note', ' (from fan command)');
            else
                ts.wind = struct('t', t, 'Y', V, 'names', {{'fan command'}}, 'ylabel', 'fan (V)', 'symbol', 'fan', 'unit', 'V', 'note', ' (command)');
            end
        end
    end

    % wave gauge
    if have(wave_path)
        if ~exist('H', 'var') || ~strcmp(wave_path, R.wind_path), H = load(wave_path); end
        if all(isfield(H, {'hwT_wind', 'eta'}))
            t = H.hwT_wind(:);  e = H.eta(:);
            n = max(1, round(avg_s / median(diff(t))));  nb = floor(numel(t) / n);
            unit = 'cm';  if isfield(H, 'etaUnit'), unit = H.etaUnit; end
            ts.wave = struct('t', mean(reshape(t(1:nb*n), n, nb), 1).', 'Y', mean(reshape(e(1:nb*n), n, nb), 1).', ...
                             'names', {{'$\eta$ wave gauge'}}, 'ylabel', sprintf('\\eta (%s)', unit), ...
                             'symbol', '$\eta$', 'unit', unit, 'note', ' (wave gauge, 1 m downstream)');
        end
    end
    save(cache, 'ts');
end
