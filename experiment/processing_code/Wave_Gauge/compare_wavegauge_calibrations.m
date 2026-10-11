function out = compare_wavegauge_calibrations(fileA, fileB, varargin)
%COMPARE_WAVEGAUGE_CALIBRATIONS  Side-by-side check of two wave gauge calibrations.
%
%   compare_wavegauge_calibrations                     % newest two found
%   compare_wavegauge_calibrations(fileA, fileB)       % explicit files
%   out = compare_wavegauge_calibrations(...)          % also return the numbers
%   compare_wavegauge_calibrations(..., 'Plot', false) % numbers only
%
% Answers the question you actually have when two calibration files exist:
% does it matter which one I use? It separates the two ways they can differ --
%
%   GAIN   (sensitivity, V/cm): a real change in the gauge or its amplifier.
%          Scales eta, so it matters for everything including wave heights.
%   OFFSET (volts at the same physical height): electrode fouling, amplifier
%          DC drift, or the water level itself having moved between the two
%          calibrations. A pure offset shifts eta by a constant, so it is
%          invisible to anything demeaned or detrended (spectrogram, Hs, the
%          slope-gauge comparison) and only matters for absolute elevation.
%
% and reports the bottom line: the cm error you incur by applying one file's
% polynomial to the other's measured voltages.
%
% INPUTS
%   fileA, fileB - paths to wave_gauge_calibration_*.mat files (each holding a
%                  'calibration' struct as saved by run_wavegauge_calibration.m).
%                  Omit both, or pass '', to use the two newest found in the
%                  Wave_Gauge folder, Run_exp, and the current folder.
%                  By convention A is the OLDER and B the NEWER file; they are
%                  swapped automatically if given the other way round, so the
%                  sign of every difference below is "older minus newer".
%
% OPTIONS
%   'Plot'   - true (default) | false
%   'NEval'  - points used to sample the two polynomials (default 2001)
%
% OUTPUT (struct, when requested)
%   .A, .B            per-file scalars (sensitivity, residuals, span, ...)
%   .dV_at_height     table of voltage difference at each shared nominal height
%   .offset_V         mean voltage offset B->A at shared heights
%   .offset_cm        that offset in cm, via the mean sensitivity
%   .gain_ratio       B sensitivity / A sensitivity
%   .eta_diff         stats of eta_A - eta_B over the shared voltage band
%   .wrong_file_cm    rms/max cm error from using the wrong file
%
% ASSUMPTIONS
%   - Both files store the same fields run_wavegauge_calibration.m writes
%     (polynomial_coefficients, voltage_calib, rel_height_calib, clipped, ...).
%     'clipped' is optional: older files predate it and are treated as unclipped.
%   - rel_height_calib is only comparable between files when height_zero and
%     sign_flip match. If they don't, this warns and the shared-height section
%     is skipped -- the polynomials are still compared, since those map volts
%     to whatever datum each file used.
%
% NOT VERIFIED
%   - Nothing here touches hardware; it reads saved files only.

    %% --- Pick the files ---------------------------------------------------
    if nargin < 1, fileA = ''; end
    if nargin < 2, fileB = ''; end
    opt = struct('Plot', true, 'NEval', 2001);
    for i = 1:2:numel(varargin)
        name = validatestring(varargin{i}, fieldnames(opt));
        opt.(name) = varargin{i+1};
    end

    if isempty(fileA) || isempty(fileB)
        found = find_cal_files();
        if numel(found) < 2
            error('Need two wave_gauge_calibration_*.mat files; found %d.', numel(found));
        end
        % found is newest-first
        if isempty(fileA), fileA = found{2}; end
        if isempty(fileB), fileB = found{1}; end
        fprintf('Using the two newest calibrations found:\n');
    end

    A = load_cal(fileA);
    B = load_cal(fileB);

    % Keep A older than B so every difference below reads "older minus newer".
    if A.cal.date > B.cal.date
        [A, B] = deal(B, A);
    end
    fprintf('  A (older): %s   %s\n', A.name, string(A.cal.date));
    fprintf('  B (newer): %s   %s\n', B.name, string(B.cal.date));

    sameDatum = A.cal.height_zero == B.cal.height_zero && ...
                A.cal.sign_flip  == B.cal.sign_flip;
    if ~sameDatum
        warning(['height_zero or sign_flip differ (A: %g/%d, B: %g/%d) -- ' ...
                 'rel heights are on different datums, so the shared-height ' ...
                 'offset is not meaningful and is skipped.'], ...
                 A.cal.height_zero, A.cal.sign_flip, B.cal.height_zero, B.cal.sign_flip);
    end

    %% --- Per-file summary -------------------------------------------------
    sA = summarize(A);
    sB = summarize(B);
    print_summary_table(sA, sB);

    %% --- Shared voltage band: compare the two polynomials -----------------
    vLo = max(min(A.v), min(B.v));
    vHi = min(max(A.v), max(B.v));
    vq  = linspace(vLo, vHi, opt.NEval)';
    etaA = polyval(A.cal.polynomial_coefficients, vq);
    etaB = polyval(B.cal.polynomial_coefficients, vq);
    dEta = etaA - etaB;

    fprintf('\nShared voltage band %.3f .. %.3f V\n', vLo, vHi);
    fprintf('  eta_A - eta_B : mean %+.4f cm, rms %.4f cm, max abs %.4f cm\n', ...
        mean(dEta), rms(dEta), max(abs(dEta)));
    % Offset vs shape: the mean is the part a demeaned analysis ignores, what
    % is left after removing it is the part that survives detrending.
    fprintf('  of which      : %+.4f cm offset, %.4f cm rms shape difference\n', ...
        mean(dEta), rms(dEta - mean(dEta)));

    %% --- Same nominal heights: the offset, measured directly --------------
    dVtable = table();
    offsetV = NaN; offsetCm = NaN;
    if sameDatum
        [hc, ia, ib] = intersect(round(A.cal.rel_height_calib(:), 3), ...
                                 round(B.cal.rel_height_calib(:), 3));
        if isempty(hc)
            fprintf('\nNo shared nominal heights between the two ladders.\n');
        else
            vA = A.cal.voltage_calib(ia);
            vB = B.cal.voltage_calib(ib);
            dVtable = table(hc, vA(:), vB(:), vA(:) - vB(:), ...
                A.clip(ia(:)) | B.clip(ib(:)), ...
                'VariableNames', {'RelHeight_cm','V_A','V_B','dV','EitherClipped'});
            keep = ~dVtable.EitherClipped;
            offsetV  = mean(dVtable.dV(keep));
            offsetCm = offsetV / mean([sA.sens_VperCm, sB.sens_VperCm]);
            fprintf('\nVoltage at the same nominal heights (%d shared, %d clipped and excluded):\n', ...
                height(dVtable), nnz(~keep));
            disp(dVtable);
            fprintf('  mean offset A-B : %+.4f V  (spread %.4f V)\n', ...
                offsetV, max(dVtable.dV(keep)) - min(dVtable.dV(keep)));
            fprintf('  as elevation    : %+.4f cm\n', offsetCm);
            % A uniform dV across heights means offset drift; a dV that trends
            % with height means the gain changed too.
            pTrend = polyfit(dVtable.RelHeight_cm(keep), dVtable.dV(keep), 1);
            fprintf('  dV trend w/ h   : %+.4f V/cm (0 = pure offset, no gain change)\n', pTrend(1));
        end
    end

    %% --- Bottom line: the cost of using the wrong file --------------------
    eAB = polyval(A.cal.polynomial_coefficients, B.v) - B.h;   % A applied to B's points
    eBA = polyval(B.cal.polynomial_coefficients, A.v) - A.h;   % B applied to A's points
    fprintf('\nUsing the WRONG file (cm error against the heights actually measured):\n');
    fprintf('  A''s poly on B''s points : mean %+.4f, rms %.4f, max abs %.4f cm\n', ...
        mean(eAB), rms(eAB), max(abs(eAB)));
    fprintf('  B''s poly on A''s points : mean %+.4f, rms %.4f, max abs %.4f cm\n', ...
        mean(eBA), rms(eBA), max(abs(eBA)));
    fprintf('  For reference, each fit''s own rms residual is %.4f / %.4f cm.\n', ...
        sA.rms_resid_cm, sB.rms_resid_cm);
    if max(abs([mean(eAB) mean(eBA)])) > 3 * max([sA.rms_resid_cm sB.rms_resid_cm])
        fprintf(['  -> The bias is several times the fit scatter: it matters for ABSOLUTE\n' ...
                 '     elevation. It is almost all offset, so demeaned/detrended results\n' ...
                 '     (spectra, Hs, slope comparison) are unaffected.\n']);
    end

    %% --- Amplitude error budget -------------------------------------------
    % The number to quote for WAVE AMPLITUDE is not the offset above. Amplitude
    % is a difference of two elevations, so a constant offset cancels exactly,
    % and so does whatever you remove by zeroing on still water at the start of
    % a run. What survives is the gain difference and the curvature difference,
    % both of which scale with how far the wave swings -- so the honest figure
    % is a PERCENTAGE of amplitude, not a fixed number of cm.
    amp = amplitude_budget(A, B);
    fprintf('\nAmplitude error from picking the other calibration (after still-water zeroing):\n');
    fprintf('   a_true (cm)   a_measured (cm)     error (cm)    error (%%)\n');
    for i = 1:numel(amp.a_true)
        fprintf('    %7.2f        %7.4f        %+9.5f     %+7.3f\n', ...
            amp.a_true(i), amp.a_meas(i), amp.err_cm(i), amp.err_pct(i));
    end
    fprintf('  -> quote about %.1f%% of amplitude (plus the %.4f cm noise floor),\n', ...
        max(abs(amp.err_pct)), mean([sA.noise_cm sB.noise_cm]));
    fprintf('     NOT the %+.3f cm offset, which cancels in any amplitude.\n', ...
        mean(dEta));

    % A height visited twice inside one calibration measures the drift rate
    % directly -- the thing still-water zeroing is protecting you from.
    for S = [A B]
        [tw, dvw, dtw] = repeat_point_drift(S);
        if tw
            fprintf('  Repeat of the same height inside %s: %+.4f V in %.1f min', S.name, dvw, dtw);
            fprintf(' = %+.4f cm (%.4f cm/min)\n', dvw / sB.sens_VperCm, dvw / sB.sens_VperCm / dtw);
        end
    end

    %% --- Rail headroom, using the newer file ------------------------------
    % The +-10 V input is the real limit on wave height, and an offset eats
    % into it asymmetrically. Worth knowing before a big-wave run.
    rail = 10;
    etaUp   = polyval(B.cal.polynomial_coefficients,  rail);
    etaDown = polyval(B.cal.polynomial_coefficients, -rail);
    fprintf('\nRail headroom on the newer file (extrapolated to +-%g V):\n', rail);
    fprintf('  eta at +%gV = %+.2f cm, at -%gV = %+.2f cm -> usable about %.2f cm p-p\n', ...
        rail, etaUp, rail, etaDown, etaUp - etaDown);
    fprintf('  (offset makes this asymmetric; re-zero the amplifier to centre it)\n');

    %% --- Order check: is the polynomial order earning its keep? -----------
    fprintf('\nIs the cubic justified? (departure from a straight line vs fit scatter)\n');
    for s = [sA sB]
        fprintf('  %-16s order %d, curvature %.4f cm, rms residual %.4f cm%s\n', ...
            s.name, s.order, s.curvature_cm, s.rms_resid_cm, ...
            ternary(s.curvature_cm < 2*s.rms_resid_cm, '   <-- at the noise level', ''));
    end
    if numel(A.cal.polynomial_coefficients) == numel(B.cal.polynomial_coefficients)
        cA = A.cal.polynomial_coefficients(1);
        cB = B.cal.polynomial_coefficients(1);
        fprintf('  leading coeff: A %.3e, B %.3e (ratio %.1fx)\n', cA, cB, max(abs([cA cB]))/min(abs([cA cB])));
        fprintf('  A leading coefficient that swings by orders of magnitude between\n');
        fprintf('  calibrations is not determined by the data -- consider a lower order.\n');
    end

    %% --- Output ----------------------------------------------------------
    if nargout > 0
        out = struct();
        out.A = sA; out.B = sB;
        out.files = {A.path, B.path};
        out.dV_at_height = dVtable;
        out.offset_V  = offsetV;
        out.offset_cm = offsetCm;
        out.gain_ratio = sB.sens_VperCm / sA.sens_VperCm;
        out.eta_diff = struct('v', vq, 'd', dEta, 'mean_cm', mean(dEta), ...
            'rms_cm', rms(dEta), 'max_cm', max(abs(dEta)), ...
            'shape_rms_cm', rms(dEta - mean(dEta)));
        out.wrong_file_cm = struct('A_on_B_rms', rms(eAB), 'A_on_B_max', max(abs(eAB)), ...
            'B_on_A_rms', rms(eBA), 'B_on_A_max', max(abs(eBA)));
        out.amplitude = amp;
    end

    if ~opt.Plot, return; end

    %% --- Plots -----------------------------------------------------------
    figure('Name', 'Wave gauge calibration comparison');
    tiledlayout(2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

    % 1. Both curves, with the fitted points
    nexttile; hold on;
    plot(vq, etaA, '-',  'LineWidth', 1.2, 'DisplayName', ['A: ' A.name]);
    plot(vq, etaB, '--', 'LineWidth', 1.2, 'DisplayName', ['B: ' B.name]);
    plot(A.v, A.h, 'o', 'MarkerSize', 5, 'DisplayName', 'A points');
    plot(B.v, B.h, 's', 'MarkerSize', 5, 'DisplayName', 'B points');
    if any(A.clip), plot(A.cal.voltage_calib(A.clip), A.cal.rel_height_calib(A.clip), 'rx', ...
            'MarkerSize', 9, 'DisplayName', 'clipped (not fitted)'); end
    if any(B.clip), plot(B.cal.voltage_calib(B.clip), B.cal.rel_height_calib(B.clip), 'rx', ...
            'MarkerSize', 9, 'HandleVisibility', 'off'); end
    hold off; grid on;
    xlabel('E_{wg} (V)'); ylabel('\eta (cm)');
    title('Calibration curves'); legend('Location','northwest', 'Interpreter','none');

    % 2. The difference, split into offset and shape
    nexttile; hold on;
    plot(vq, dEta, 'LineWidth', 1.2, 'DisplayName', '\eta_A - \eta_B');
    yline(mean(dEta), 'k--', 'DisplayName', sprintf('mean %+.3f cm', mean(dEta)));
    hold off; grid on;
    xlabel('E_{wg} (V)'); ylabel('\Delta\eta (cm)');
    title('A minus B (flat = pure offset)'); legend('Location','best');

    % 3. Residuals of each fit, on height
    nexttile; hold on;
    plot(A.h, A.resid, 'o-', 'DisplayName', sprintf('A (rms %.4f cm)', sA.rms_resid_cm));
    plot(B.h, B.resid, 's-', 'DisplayName', sprintf('B (rms %.4f cm)', sB.rms_resid_cm));
    yline(0, 'k:', 'HandleVisibility','off');
    hold off; grid on;
    xlabel('\eta (cm)'); ylabel('residual (cm)');
    title('Fit residuals'); legend('Location','best');

    % 4. Offset at each shared height
    nexttile;
    if ~isempty(dVtable)
        keep = ~dVtable.EitherClipped;
        plot(dVtable.RelHeight_cm(keep), dVtable.dV(keep), 'o-', 'LineWidth', 1.2); hold on;
        yline(offsetV, 'k--', 'LineWidth', 1);
        hold off; grid on;
        xlabel('\eta (cm)'); ylabel('V_A - V_B (V)');
        title(sprintf('Offset at shared heights (mean %+.3f V = %+.3f cm)', offsetV, offsetCm));
    else
        axis off;
        text(0.5, 0.5, 'no shared nominal heights', 'HorizontalAlignment','center');
    end
end

%% ------------------------------------------------------------------------
function files = find_cal_files()
% Newest-first list of calibration files, same search path as the run scripts.
    dirs = unique({fileparts(mfilename('fullpath')), ...
                   'D:\Chris\osbl-turbulent-mixing\experiment\processing_code\Run_exp', ...
                   pwd});
    f = [];
    for i = 1:numel(dirs)
        f = [f; dir(fullfile(dirs{i}, 'wave_gauge_calibration_*.mat'))]; %#ok<AGROW>
    end
    files = {};
    if isempty(f), return; end
    [~, ord] = sort([f.datenum], 'descend');
    f = f(ord);
    files = arrayfun(@(x) fullfile(x.folder, x.name), f, 'UniformOutput', false);
end

function S = load_cal(path)
% Load one file and pre-compute the fitted (unclipped) subset.
    if ~isfile(path)
        error('Calibration file not found: %s', path);
    end
    L = load(path, 'calibration');
    if ~isfield(L, 'calibration')
        error('%s has no ''calibration'' struct.', path);
    end
    S.path = path;
    [~, S.name] = fileparts(path);
    S.name = erase(S.name, 'wave_gauge_calibration_');
    S.cal  = L.calibration;
    if isfield(S.cal, 'clipped')
        S.clip = logical(S.cal.clipped(:));
    else
        % Older files predate the clipped flag -- nothing was excluded then.
        S.clip = false(numel(S.cal.voltage_calib), 1);
    end
    S.v = S.cal.voltage_calib(~S.clip);
    S.h = S.cal.rel_height_calib(~S.clip);
    S.v = S.v(:); S.h = S.h(:);
    S.resid = S.h - polyval(S.cal.polynomial_coefficients, S.v);
end

function s = summarize(S)
    p1 = polyfit(S.v, S.h, 1);          % cm per V
    s.name         = S.name;
    s.date         = S.cal.date;
    s.n            = numel(S.cal.voltage_calib);
    s.n_clipped    = nnz(S.clip);
    s.order        = S.cal.polynomial_order;
    s.h_min        = min(S.h);
    s.h_max        = max(S.h);
    s.v_min        = min(S.v);
    s.v_max        = max(S.v);
    s.sens_VperCm  = 1 / p1(1);
    s.rms_resid_cm = rms(S.resid);
    s.max_resid_cm = max(abs(S.resid));
    s.noise_V      = mean(S.cal.voltage_std);
    s.noise_cm     = s.noise_V * abs(p1(1));
    s.v_at_zero    = interp_unique(S.h, S.v, 0);
    % How far the fitted polynomial departs from a straight line over its own
    % span: the curvature the extra orders are buying.
    s.curvature_cm = max(abs(polyval(S.cal.polynomial_coefficients, S.v) - polyval(p1, S.v)));
end

function print_summary_table(sA, sB)
    rows = {
        'date',                  @(s) char(string(s.date))
        'points (clipped)',      @(s) sprintf('%d (%d)', s.n, s.n_clipped)
        'eta span (cm)',         @(s) sprintf('%.2f .. %.2f', s.h_min, s.h_max)
        'voltage span (V)',      @(s) sprintf('%.3f .. %.3f', s.v_min, s.v_max)
        'sensitivity (V/cm)',    @(s) sprintf('%.4f', s.sens_VperCm)
        'V at eta = 0',          @(s) sprintf('%+.4f', s.v_at_zero)
        'rms residual (cm)',     @(s) sprintf('%.4f', s.rms_resid_cm)
        'max residual (cm)',     @(s) sprintf('%.4f', s.max_resid_cm)
        'noise (V / cm)',        @(s) sprintf('%.4f / %.4f', s.noise_V, s.noise_cm)
        'poly order',            @(s) sprintf('%d', s.order)
        'curvature (cm)',        @(s) sprintf('%.4f', s.curvature_cm) };
    w = 22;
    fprintf('\n%-*s %-26s %-26s\n', w, '', ['A: ' sA.name], ['B: ' sB.name]);
    fprintf('%s\n', repmat('-', 1, w + 54));
    for i = 1:size(rows, 1)
        fprintf('%-*s %-26s %-26s\n', w, rows{i,1}, rows{i,2}(sA), rows{i,2}(sB));
    end
    fprintf('\nsensitivity ratio B/A : %.4f (%.2f%% change)\n', ...
        sB.sens_VperCm / sA.sens_VperCm, 100*(sB.sens_VperCm/sA.sens_VperCm - 1));
end

function vq = interp_unique(h, v, hq)
% interp1 refuses duplicate sample points, and a repeated height is normal
% (run_wavegauge_calibration.m allows revisiting a point), so average them.
    [hu, ~, iu] = unique(h(:));
    vu = accumarray(iu, v(:), [], @mean);
    if numel(hu) < 2
        vq = vu(1);
    else
        vq = interp1(hu, vu, hq, 'linear', 'extrap');
    end
end

function amp = amplitude_budget(A, B)
% Synthetic sinusoid of known amplitude, built in VOLTS through the newer
% calibration (B = the one in force when the data would be taken), then read
% back with each file and demeaned -- exactly what a run does when it zeroes on
% still water. The leftover difference is the amplitude error.
    pA = A.cal.polynomial_coefficients;
    pB = B.cal.polynomial_coefficients;

    % Numerical inverse of B: eta -> volts, over B's own fitted span.
    vGrid = linspace(min(B.v), max(B.v), 20001)';
    eGrid = polyval(pB, vGrid);
    halfSpan = min(abs([min(eGrid) max(eGrid)]));   % biggest symmetric wave that fits

    fracs = [0.1 0.25 0.5 0.75 1.0];
    t = linspace(0, 1, 4096)';
    amp = struct('a_true', [], 'a_meas', [], 'err_cm', [], 'err_pct', []);
    for f = fracs
        a = f * halfSpan;
        eTrue = a * sin(2*pi*t);
        vWave = interp1(eGrid, vGrid, eTrue, 'linear', 'extrap');
        eA = polyval(pA, vWave); eA = eA - mean(eA);
        amp.a_true(end+1)  = a;
        amp.a_meas(end+1)  = (max(eA) - min(eA)) / 2;
        amp.err_cm(end+1)  = amp.a_meas(end) - a;
        amp.err_pct(end+1) = 100 * amp.err_cm(end) / a;
    end
end

function [found, dv, dtMin] = repeat_point_drift(S)
% If one nominal height was measured twice in the same calibration, the change
% in voltage between those visits is a direct measurement of drift rate.
    found = false; dv = NaN; dtMin = NaN;
    if ~isfield(S.cal, 'point_time'), return, end
    h = round(S.cal.rel_height_calib(:), 3);
    [hu, ~, iu] = unique(h);
    for k = 1:numel(hu)
        idx = find(iu == k);
        if numel(idx) >= 2
            dv    = S.cal.voltage_calib(idx(end)) - S.cal.voltage_calib(idx(1));
            dtMin = minutes(S.cal.point_time(idx(end)) - S.cal.point_time(idx(1)));
            found = abs(dtMin) > 0;
            return
        end
    end
end

function s = ternary(cond, a, b)
    if cond, s = a; else, s = b; end
end
