function [camDelay, hw, wave] = lif_wind_clock(R, avg_s, cal_dir)
%LIF_WIND_CLOCK Camera start after wind start + wind and wave series for a run.
%
%   [camDelay, hw, wave] = LIF_WIND_CLOCK(R, avg_s, cal_dir)
%
% R from lif_run_config.  With a wind record (hot-wire, reference probe or
% wave-gauge file) camDelay and hw come from it (lif_hw_cached).  Without
% one (registry wind_kind 'none', e.g. strat runs 116-120), camDelay is the
% registry cam_delay_s (10 s: every recorded run shows 10.02-10.10 s) and
% hw is R.standin_file -- another run's wind; hw = [] if there is none.
%
% wave - wave-gauge eta (cm) for a panel under the wind: R.wave_file when
%        the run has one (and the wind panel is not already that gauge),
%        else R.wave_standin_file (strat runs: wave_gauge_20261008_231731_r1); [] if neither.
%
% Legend names say where a record comes from, as in the CISG videos:
% "(simultaneous)" for this run's own record, "(stand-in: core N)" for
% another run's (N looked up in the registry), never the file name.

    camDelay = [];  hw = [];  wave = [];
    if ~isempty(R.hw_file)
        [~, nm] = fileparts(R.hw_file);
        [camDelay, hw] = lif_hw_cached(R.hw_file, avg_s, fullfile(cal_dir, ...
                                       sprintf('hwsync_%s_avg%g.mat', nm, avg_s)));
        hw.names = label(hw.names, '(simultaneous)');
    else
        camDelay = R.cam_delay_s;
        if isempty(camDelay) || isnan(camDelay)
            error('lif_wind_clock:none', 'Run %d has no wind record and no cam_delay_s in the registry.', R.run);
        end
        fprintf('Run %d: no wind record -- camera start assumed %.2f s after wind start.\n', R.run, camDelay);
        if ~isempty(R.standin_file)
            [~, nm] = fileparts(R.standin_file);
            [~, hw] = lif_hw_cached(R.standin_file, avg_s, fullfile(cal_dir, ...
                                    sprintf('hwsync_%s_avg%g.mat', nm, avg_s)));
            hw.names = label(hw.names, sprintf('(stand-in: %s)', source_run(nm)));
            hw.standin = true;
        end
    end

    % wave gauge
    own = isfield(R, 'wave_file') && ~isempty(R.wave_file) && ...
          ~(~isempty(hw) && strcmp(hw.kind, 'eta'));       % wind panel already shows this gauge
    if own
        wave = wave_series(R.wave_file, avg_s, cal_dir);
        if ~isempty(wave), wave.names = label(wave.names, '(simultaneous)'); end
    elseif isfield(R, 'wave_standin_file') && ~isempty(R.wave_standin_file)
        wave = wave_series(R.wave_standin_file, avg_s, cal_dir);
        if ~isempty(wave)
            [~, nm] = fileparts(R.wave_standin_file);
            wave.names = label(wave.names, sprintf('(stand-in: %s)', source_run(nm)));
            wave.standin = true;
        end
    end
end

function names = label(names, tag)
    names = strrep(names, ' (hot-wire)', ' hot-wire');
    names = strrep(names, ' (wave gauge)', ' wave gauge');
    names = cellfun(@(s) [s ' ' tag], names, 'UniformOutput', false);
end

function s = source_run(nm)
    % the registry run whose own record this file is ('core N'); the file
    % name only when no run claims it
    T = run_registry();
    r = T.run(strcmp(T.wind_file, nm) | strcmp(T.wave_file, nm));
    if isempty(r), s = strrep(nm, '_', '\_'); else, s = sprintf('core %d', r(1)); end
end

function wv = wave_series(f, avg_s, cal_dir)
    % eta (cm) vs s since wind start, block-averaged; cached for the other computer
    [~, nm] = fileparts(f);
    cache = fullfile(cal_dir, sprintf('wave_%s_avg%g.mat', nm, avg_s));
    wv = [];
    if isfile(f)
        H = load(f);
        if ~all(isfield(H, {'hwT_wind', 'eta'}))
            warning('lif_wind_clock:wave', '%s has no hwT_wind / eta.', f);  return
        end
        t = H.hwT_wind(:);  e = H.eta(:);
        n = max(1, round(avg_s / median(diff(t))));  nb = floor(numel(t) / n);
        unit = 'cm';  if isfield(H, 'etaUnit'), unit = H.etaUnit; end
        wv = struct('t', mean(reshape(t(1:nb*n), n, nb), 1).', 'Y', mean(reshape(e(1:nb*n), n, nb), 1).', ...
                    'names', {{'$\eta$ (wave gauge)'}}, 'ylabel', sprintf('\\eta (%s)', unit), 'kind', 'eta');
        if ~isfolder(cal_dir), mkdir(cal_dir); end
        save(cache, 'wv');
    elseif isfile(cache)
        S = load(cache, 'wv');  wv = S.wv;
    else
        warning('lif_wind_clock:wave', 'Neither %s nor its cache %s exists.', f, cache);
    end
end
