%% plot_midheight_stress.m
% Air-side stress and turbulence from the hotwire_midheight runs.
%
% Run section by section (Ctrl+Enter) to follow each step, or all at once (F5).
%   Step 1  find the files and group them by height
%   Step 2  load each run, keep only the ramp-up        -> check figure
%   Step 3  stats in short time windows, per run         -> check figure (one window)
%   Step 4  average the repeats at each height           -> printed table
%   Step 5  surface stress u* from the lowest heights    -> printed table
%   Fig A   stress vs time, and u*(t) vs Wagner et al.
%   Fig B   timeline at one height (U, turbulence, stress)
%   Fig C   spectra before / after inception

%% Settings
DATA_DIR      = 'D:\Chris\osbl-turbulent-mixing\experiment\processing_code\Run_exp';
AIR_CENTER_CM = 16.85;     % cm above water; file "<k>cm" = AIR_CENTER_CM - k above water
W_SIGN        = -1;        % probe W looks positive DOWN; -1 makes it positive up (confirm with the flip test)
WIN           = 2;         % s, length of each averaging window
STEP          = 0.5;       % s, spacing between windows
SURF_Z_MAX    = 3.9;       % cm, heights at or below this are averaged into the surface u*
TIMELINE_Z    = 2.85;      % cm, height used for Fig B
SPEC_Z        = 2.85;      % cm, height used for Fig C
SPEC_WINDOWS  = [14 20; 25 31; 40 46];   % s, windows compared in Fig C
T_AIR_TRANS   = 12;        % s, air goes laminar -> turbulent (measured in these runs)
T_INCEPTION   = [22.7 24.4];   % s, Langmuir inception from IR frames 637-721
% Wagner et al. (2023) water-side stress tau = ALPHA*sqrt(t - T0), converted to air:
ALPHA = 0.12e-4;           % m^2 s^-5/2 (0.12 cm^2 s^-5/2)
T0    = 6.5;               % s, their t = 0 on our clock -- ASSUMED (inception 22.7 s vs their 16 s)
RHO_A = 1.2;  RHO_W = 998;
OUT_DIR = '';              % folder to save PNGs, '' = don't save

%% Step 1: find the files and group them by height
files = dir(fullfile(DATA_DIR, 'hotwire_midheight-*cm*.mat'));
ok = ~cellfun(@isempty, regexp({files.name}, '^hotwire_midheight-[\d.]+cm(_stopT[\d.]+s)?(_r\d+)?\.mat$', 'once'));
files = files(ok);
fileZ = zeros(numel(files), 1);
for i = 1:numel(files)
    k = str2double(regexp(files(i).name, '(?<=midheight-)[\d.]+(?=cm)', 'match', 'once'));
    fileZ(i) = round(AIR_CENTER_CM - k, 4);   % round so 16.85-14 equals 2.85 exactly
end
heights = unique(fileZ);
fprintf('Step 1: %d files at %d heights\n', numel(files), numel(heights));
for h = 1:numel(heights)
    fprintf('  %5.2f cm: %d run(s)\n', heights(h), sum(fileZ == heights(h)));
end

%% Step 2: load each run, keep only the ramp-up
% Time axis: seconds since wind start (rebuilt from hwT; the saved hwT_wind in
% older files is ~0.1 s off). Data after the ramp's peak (the ramp-down) is
% dropped -- it's not comparable between runs cut at different times.
R = struct('name', {}, 'z', {}, 't', {}, 'U', {}, 'W', {}, 'tPeak', {});
for i = 1:numel(files)
    S = load(fullfile(DATA_DIR, files(i).name), 'hwT', 'hwStartElapsed', 'fanStartElapsed', 'U', 'W', 'sentT', 'sentV');
    t = S.hwT(:) + S.hwStartElapsed - S.fanStartElapsed;
    [~, iPk] = max(S.sentV);
    tPeak = S.sentT(iPk) - S.fanStartElapsed;
    keep = t >= 0 & t <= tPeak;
    R(i) = struct('name', files(i).name, 'z', fileZ(i), 't', t(keep), ...
                  'U', S.U(keep), 'W', W_SIGN * S.W(keep), 'tPeak', tPeak);
    fprintf('Step 2: %-42s z = %5.2f cm, ramp peak at %4.1f s\n', R(i).name, R(i).z, tPeak);
end

% Check figure: one run's raw signals
iEx = find([R.z] == TIMELINE_Z, 1);
figure('Name', 'Step 2 check: raw signals', 'Color', 'w');
subplot(2,1,1); plot(R(iEx).t, R(iEx).U); ylabel('U (m/s)'); grid on;
title(sprintf('Step 2 check: %s', R(iEx).name), 'Interpreter', 'none');
subplot(2,1,2); plot(R(iEx).t, R(iEx).W); ylabel('W, up positive (m/s)'); grid on;
xlabel('Time since wind start (s)');

%% Step 3: statistics in short time windows, for every run
% In each WIN-second window:
%   (a) rotate u,w so the mean vertical velocity is zero (removes probe tilt)
%   (b) remove a straight-line trend (the ramp), leaving the fluctuations u', w'
%   (c) U = mean speed, u_rms = std(u'), stress = -mean(u'w')
tc = 0:STEP:max([R.tPeak]);          % window centres, common to all runs
for i = 1:numel(R)
    [R(i).Um, R(i).urms, R(i).wrms, R(i).tau] = deal(nan(size(tc)));
    for k = 1:numel(tc)
        if tc(k) - WIN/2 < 0 || tc(k) + WIN/2 > R(i).tPeak, continue; end   % window must fit in the ramp
        m  = abs(R(i).t - tc(k)) <= WIN/2;
        tt = R(i).t(m);  u = R(i).U(m);  w = R(i).W(m);
        th = atan2(mean(w), mean(u));                                        % (a)
        ur =  u*cos(th) + w*sin(th);
        wr = -u*sin(th) + w*cos(th);
        up = ur - polyval(polyfit(tt, ur, 1), tt);                          % (b)
        wp = wr - polyval(polyfit(tt, wr, 1), tt);
        R(i).Um(k) = mean(u);  R(i).urms(k) = std(up);  R(i).wrms(k) = std(wp);   % (c)
        R(i).tau(k) = -mean(up .* wp);
    end
end
fprintf('Step 3: window stats done (%d windows of %g s per run)\n', numel(tc), WIN);

% Check figure: what one window looks like (run from Step 2, window at 25 s)
m  = abs(R(iEx).t - 25) <= WIN/2;
tt = R(iEx).t(m);  u = R(iEx).U(m);  w = R(iEx).W(m);
th = atan2(mean(w), mean(u));
ur =  u*cos(th) + w*sin(th);   wr = -u*sin(th) + w*cos(th);
up = ur - polyval(polyfit(tt, ur, 1), tt);   wp = wr - polyval(polyfit(tt, wr, 1), tt);
figure('Name', 'Step 3 check: one window', 'Color', 'w', 'Position', [100 100 900 350]);
subplot(1,3,1); plot(tt, ur, tt, polyval(polyfit(tt, ur, 1), tt), 'LineWidth', 1.5);
xlabel('t (s)'); ylabel('u (m/s)'); title('u and its trend'); grid on;
subplot(1,3,2); plot(tt, wr, tt, polyval(polyfit(tt, wr, 1), tt), 'LineWidth', 1.5);
xlabel('t (s)'); ylabel('w (m/s)'); title(sprintf('w and its trend (tilt %.2f deg)', rad2deg(th))); grid on;
subplot(1,3,3); plot(up(1:20:end), wp(1:20:end), '.', 'MarkerSize', 4);
xlabel('u'' (m/s)'); ylabel('w'' (m/s)'); axis equal; grid on;
title(sprintf('-u''w'' = %.4f m^2/s^2', -mean(up.*wp)));

%% Step 4: average the repeats at each height
% H(h).tau = mean over runs, H(h).tauStd = spread between runs, H(h).n = runs used
H = struct('z', num2cell(heights));
for h = 1:numel(heights)
    Rh = R([R.z] == heights(h));
    for f = ["Um", "urms", "wrms", "tau"]
        A = vertcat(Rh.(f));
        H(h).(f) = mean(A, 1, 'omitnan');
        H(h).(f + "Std") = std(A, 0, 1, 'omitnan');
    end
    H(h).n = sum(~isnan(vertcat(Rh.tau)), 1);
end
showT = [11 13 20 23 27 33];
fprintf('\nStep 4: -u''w'' (x1e-3 m^2/s^2), mean over runs [runs used]\n%8s', 'z (cm)');
for s = showT, fprintf('%11s', sprintf('t=%gs', s)); end
fprintf('\n');
for h = 1:numel(H)
    fprintf('%8.2f', H(h).z);
    for s = showT
        k = find(abs(tc - s) < 1e-9);
        fprintf('%7.2f [%d]', 1e3*H(h).tau(k), H(h).n(k));
    end
    fprintf('\n');
end

%% Step 5: surface stress u* from the lowest heights
% Early in the ramp the stress hardly changes between 1.85 and 3.85 cm, so the
% surface stress is taken as the average over those heights (all runs pooled).
% +/- is the standard error of that average.
low = [R.z] <= SURF_Z_MAX;
A = vertcat(R(low).tau);
tauS   = mean(A, 1, 'omitnan');
tauSE  = std(A, 0, 1, 'omitnan') ./ sqrt(sum(~isnan(A), 1));
% (not max(x,0): MATLAB's max turns NaN into 0, which would plot "no data" as zero stress)
posPart = @(x) x .* (x > 0 | isnan(x));        % negative -> 0, NaN stays NaN
ustar   = sqrt(posPart(tauS));
ustarLo = sqrt(posPart(tauS - tauSE));   ustarHi = sqrt(posPart(tauS + tauSE));
% Note: which heights contribute changes with time (1.85 cm runs end at 35 s,
% 2.85 cm at 50 s, 3.85 cm at 60 s), so small steps at those times are expected.
tauW   = (RHO_W/RHO_A) * ALPHA * sqrt(max(tc - T0, 0));   % Wagner et al., air side
fprintf('\nStep 5: surface u* from %d runs at z <= %.2f cm\n', sum(low), SURF_Z_MAX);
fprintf('%8s %10s %14s\n', 't (s)', 'u* (m/s)', 'Wagner (m/s)');
for s = showT
    k = find(abs(tc - s) < 1e-9);
    fprintf('%8g %10.3f %14.3f\n', s, ustar(k), sqrt(tauW(k)));
end

%% Shared look for the figures
ink = [0.2 0.2 0.2];  muted = [0.45 0.45 0.45];  bandGrey = [0.92 0.92 0.92];
blue = [42 120 214]/255;  orange = [235 104 52]/255;  aqua = [27 175 122]/255;
rampHex = {'86b6ef','5598e7','2a78d6','1c5cab','104281','0d366b'};   % light = high, dark = near water
ramp = cell2mat(cellfun(@(x) sscanf(x, '%2x%2x%2x', [1 3])/255, rampHex(:), 'UniformOutput', false));
hCol = interp1(linspace(0,1,size(ramp,1)), ramp, linspace(1,0,numel(heights)));   % row h = colour for heights(h)
styleAx = @(ax) set(ax, 'Box', 'off', 'TickDir', 'out', 'XColor', muted, 'YColor', muted, ...
    'GridColor', [0.85 0.85 0.85], 'GridAlpha', 1, 'FontSize', 10);
markEvents = @(ax) mark_events(ax, T_AIR_TRANS, T_INCEPTION, muted, bandGrey);   % local function at the end

%% Fig A: stress vs time at each height, and surface u* vs Wagner et al.
figA = figure('Name', 'Fig A: stress', 'Color', 'w', 'Position', [80 80 900 650]);
axA1 = subplot(2,1,1); hold on;
markEvents(axA1);
hA = gobjects(numel(H), 1);
for h = numel(H):-1:1      % highest first, so the legend reads top-down like the tank
    hA(h) = plot(tc, 1e3*H(h).tau, 'Color', hCol(h,:), 'LineWidth', 1.5, ...
        'DisplayName', sprintf('%.2f cm (%d runs)', H(h).z, max(H(h).n)));
end
hold off; grid on; styleAx(axA1);
ylabel('-u''w'' (\times10^{-3} m^2/s^2)'); xlim([0 max(tc)]);
legend(flipud(hA), 'Location', 'northwest', 'Box', 'off', 'TextColor', ink);
title('Stress at each height (grey band = Langmuir inception)', 'FontWeight', 'normal', 'Color', ink);

axA2 = subplot(2,1,2); hold on;
markEvents(axA2);
ok5 = ~isnan(ustar);
fill([tc(ok5) fliplr(tc(ok5))], [ustarLo(ok5) fliplr(ustarHi(ok5))], blue, ...
    'FaceAlpha', 0.2, 'EdgeColor', 'none', 'HandleVisibility', 'off');
h1 = plot(tc, ustar, 'Color', blue, 'LineWidth', 2, 'DisplayName', sprintf('measured (z \\leq %.1f cm)', SURF_Z_MAX));
h2 = plot(tc, sqrt(tauW), '--', 'Color', orange, 'LineWidth', 2, ...
    'DisplayName', sprintf('Wagner et al. \\alpha\\surd(t - %.1f s), air side', T0));
hold off; grid on; styleAx(axA2);
xlabel('Time since wind start (s)'); ylabel('u* (m/s)'); xlim([0 max(tc)]);
legend([h1 h2], 'Location', 'northwest', 'Box', 'off', 'TextColor', ink);
title('Surface friction velocity (band = \pm1 standard error)', 'FontWeight', 'normal', 'Color', ink);

%% Fig B: timeline at one height
hB = find(heights == TIMELINE_Z);
figB = figure('Name', 'Fig B: timeline', 'Color', 'w', 'Position', [120 60 800 700]);
qty   = {'Um', 'urms', 'tau'};
scale = [1 1 1e3];
lab   = {'U (m/s)', 'u_{rms} (m/s)', '-u''w'' (\times10^{-3} m^2/s^2)'};
for p = 1:3
    ax = subplot(3,1,p); hold on;
    markEvents(ax);
    y  = scale(p) * H(hB).(qty{p});
    sd = scale(p) * H(hB).(qty{p} + "Std");
    okB = ~isnan(y) & ~isnan(sd);
    fill([tc(okB) fliplr(tc(okB))], [y(okB)-sd(okB) fliplr(y(okB)+sd(okB))], blue, ...
        'FaceAlpha', 0.2, 'EdgeColor', 'none');
    plot(tc, y, 'Color', blue, 'LineWidth', 2);
    hold off; grid on; styleAx(ax); ylabel(lab{p}); xlim([0 max(tc)]);
    if p == 1
        title(sprintf('Timeline at %.2f cm (%d runs; band = \\pm1 std between runs)', ...
            TIMELINE_Z, max(H(hB).n)), 'FontWeight', 'normal', 'Color', ink);
    end
end
xlabel('Time since wind start (s)');

%% Fig C: spectra before / after inception
% Note: ripples at inception (~3 cm long) are invisible at these heights --
% their air motion dies off like exp(-k z), ~2% left at 1.85 cm. Changes here
% reflect the air turbulence and, later, longer waves.
Rs = R([R.z] == SPEC_Z);
winCol = [blue; orange; aqua];
figC = figure('Name', 'Fig C: spectra', 'Color', 'w', 'Position', [160 80 950 420]);
axC1 = subplot(1,2,1); axC2 = subplot(1,2,2);
hold(axC1, 'on'); hold(axC2, 'on');
fprintf('\nFig C: spectra at %.2f cm\n', SPEC_Z);
for w = 1:size(SPEC_WINDOWS, 1)
    Puu = []; Co = []; nUsed = 0;
    for i = 1:numel(Rs)
        m = Rs(i).t >= SPEC_WINDOWS(w,1) & Rs(i).t <= SPEC_WINDOWS(w,2);
        if Rs(i).tPeak < SPEC_WINDOWS(w,2) || ~any(m), continue; end   % window must be inside this run's ramp
        tt = Rs(i).t(m);  fs = 1/median(diff(tt));
        u = detrend(Rs(i).U(m));  wv = detrend(Rs(i).W(m));
        nseg = round(2*fs);                                            % 2 s segments
        [p, f] = pwelch(u, hann(nseg), nseg/2, nseg, fs);
        c = real(cpsd(u, -wv, hann(nseg), nseg/2, nseg, fs));        % -u'w' cospectrum
        Puu = [Puu, p]; Co = [Co, c]; nUsed = nUsed + 1;  %#ok<AGROW>
    end
    fprintf('  %2d-%2d s: %d run(s)\n', SPEC_WINDOWS(w,1), SPEC_WINDOWS(w,2), nUsed);
    if nUsed == 0, continue; end
    nm = sprintf('%d-%d s (%d runs)', SPEC_WINDOWS(w,1), SPEC_WINDOWS(w,2), nUsed);
    % Average into log-spaced bands (10 per decade) -- raw 0.5 Hz bins are too noisy to read
    [fb, Pb] = log_bin(f, mean(Puu, 2), 10);
    [~,  Cb] = log_bin(f, mean(Co, 2), 10);
    loglog(axC1, fb, Pb, '-o', 'Color', winCol(w,:), 'LineWidth', 1.5, 'MarkerSize', 4, ...
        'MarkerFaceColor', winCol(w,:), 'DisplayName', nm);
    semilogx(axC2, fb, fb .* Cb, '-o', 'Color', winCol(w,:), 'LineWidth', 1.5, 'MarkerSize', 4, ...
        'MarkerFaceColor', winCol(w,:), 'DisplayName', nm);
end
set(axC1, 'XScale', 'log', 'YScale', 'log'); set(axC2, 'XScale', 'log');
grid(axC1, 'on'); grid(axC2, 'on'); styleAx(axC1); styleAx(axC2);
xlabel(axC1, 'f (Hz)'); ylabel(axC1, 'S_{uu} (m^2 s^{-2} Hz^{-1})');
title(axC1, sprintf('u spectrum at %.2f cm', SPEC_Z), 'FontWeight', 'normal', 'Color', ink);
xlabel(axC2, 'f (Hz)'); ylabel(axC2, 'f \cdot Co_{-uw} (m^2/s^2)');
title(axC2, 'Stress cospectrum (area = stress)', 'FontWeight', 'normal', 'Color', ink);
legend(axC1, 'Location', 'southwest', 'Box', 'off', 'TextColor', ink);
yline(axC2, 0, 'Color', muted, 'HandleVisibility', 'off');

%% Save
if ~isempty(OUT_DIR)
    exportgraphics(figA, fullfile(OUT_DIR, 'midheight_stress_A.png'), 'Resolution', 150);
    exportgraphics(figB, fullfile(OUT_DIR, 'midheight_stress_B_timeline.png'), 'Resolution', 150);
    exportgraphics(figC, fullfile(OUT_DIR, 'midheight_stress_C_spectra.png'), 'Resolution', 150);
    fprintf('Saved figures to %s\n', OUT_DIR);
end


%% Local functions
function mark_events(ax, tAir, tInc, lineCol, bandCol)
% Dashed line where the air turns turbulent, grey band over Langmuir inception.
% Both kept out of legends.
xregion(ax, tInc(1), tInc(2), 'FaceColor', bandCol, 'FaceAlpha', 1, 'HandleVisibility', 'off');
xline(ax, tAir, '--', 'air turbulent', 'Color', lineCol, ...
    'LabelVerticalAlignment', 'bottom', 'HandleVisibility', 'off');
end

function [fb, yb] = log_bin(f, y, perDecade)
% Average y into log-spaced frequency bands (perDecade bands per factor of 10).
% Skips f = 0. fb = geometric centre of each band that has data.
keep = f > 0;  f = f(keep);  y = y(keep);
edges = 10.^(floor(log10(f(1))) : 1/perDecade : ceil(log10(f(end))));
fb = []; yb = [];
for b = 1:numel(edges)-1
    m = f >= edges(b) & f < edges(b+1);
    if any(m)
        fb(end+1, 1) = sqrt(edges(b) * edges(b+1)); %#ok<AGROW>
        yb(end+1, 1) = mean(y(m));                  %#ok<AGROW>
    end
end
end
