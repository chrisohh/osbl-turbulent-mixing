%% plot_midheight_profile.m
% Velocity profile from the hotwire_midheight-<z>cm runs: one run per probe
% height, each with the same fan ramp, so time since wind start is the common
% axis that lines the heights up.
%
% Two views, because a single time-averaged profile would blend the whole
% ramp (2 -> 10 m/s) into one meaningless number per height:
%   1) time vs height map, colour = U  -- overview of how the profile evolves
%   2) profiles U(z) at a few times     -- the shape, readable quantitatively
%
% Caveat: each height is a SEPARATE run, so the map assumes the ramp was
% repeatable run to run -- it is not a simultaneous measurement.

DATA_DIR = 'D:\Chris\osbl-turbulent-mixing\experiment\processing_code\Run_exp';
% The <z>cm in each filename is the distance BELOW the air-side centre height,
% not height above the water: height above water = AIR_CENTER_CM - <z>.
AIR_CENTER_CM = 16.85;           % cm above the water surface
% Averaging windows: 1 s was too short -- three repeat runs at 4.85 cm
% differed by up to 6.7% with 1 s means (turbulence noise ~0.13-0.20 m/s per
% 1 s mean) but agreed within ~1.5% with 3 s. Centred on a near-linear ramp,
% a wider window doesn't bias the value, it only blurs time a little.
SMOOTH_T = 3.0;                  % s, centred moving mean for the map
PROFILE_T = [5 15 25 35 45 55 65];   % s since wind start, profile snapshots
PROFILE_WIN = 3.0;               % s, centred averaging window around each snapshot
OUT_PNG = '';                    % set a path to also save the figures as PNG
% When one height has several files (a full run plus _stopT50s/_stopT60s cut
% runs): 'average' = average every run that is still ramping at each time
% (a cut run is an identical repeat up to its cut, so this uses all of it and
% lowers the noise); 'longest' / 'shortest' = use just one file.
PICK = 'average';

% Accepted names: hotwire_midheight-<k>cm[_stopT<N>s][_r<n>].mat -- full or
% cut runs, plus numbered repeats. Anything else (e.g. a _flip probe-rotation
% test) is left out on purpose so it can't get averaged in.
files = dir(fullfile(DATA_DIR, 'hotwire_midheight-*cm*.mat'));
ok = ~cellfun(@isempty, regexp({files.name}, '^hotwire_midheight-[\d.]+cm(_stopT[\d.]+s)?(_r\d+)?\.mat$', 'once'));
if any(~ok)
    fprintf('Not used (name doesn''t match the pattern): %s\n', strjoin({files(~ok).name}, ', '));
end
files = files(ok);
if isempty(files), error('No hotwire_midheight files in %s', DATA_DIR); end

% Load everything, and find where each run's ramp actually peaked from the
% commanded voltage (sentV). Data after that is the ramp-DOWN, which is not
% the same evolution as another run at the same time, so it's dropped. Works
% for full and cut runs alike, without trusting the filename's stop time.
cands = struct('name', {}, 'offset', {}, 't', {}, 'U', {}, 'Us', {}, 'tPeak', {});
for i = 1:numel(files)
    S = load(fullfile(DATA_DIR, files(i).name), 'hwT', 'hwStartElapsed', 'fanStartElapsed', 'U', 'sentT', 'sentV');
    % Wind clock rebuilt from hwT (0 = first fan write); the saved hwT_wind
    % in these files was overwritten by hwT in run_experiment.m's plot block.
    t = S.hwT(:) + S.hwStartElapsed - S.fanStartElapsed;
    [~, iPk] = max(S.sentV);
    tPeak = S.sentT(iPk) - S.fanStartElapsed;
    keep = t <= tPeak;
    fs = 1/median(diff(t));
    cands(end+1) = struct('name', files(i).name, ...
        'offset', str2double(regexp(files(i).name, '(?<=midheight-)[\d.]+(?=cm)', 'match', 'once')), ...
        't', t(keep), 'U', S.U(keep), 'Us', movmean(S.U(keep), max(round(SMOOTH_T*fs),1)), ...
        'tPeak', tPeak); %#ok<SAGROW>
end

% Group runs by height. runs{i} is an array of every run used at that height.
offsets = unique([cands.offset]);
runs = cell(numel(offsets), 1);
z = AIR_CENTER_CM - offsets(:);
fprintf('Height above water  runs used (ramp peak)\n');
for i = 1:numel(offsets)
    c = cands([cands.offset] == offsets(i));
    [~, j] = sort([c.tPeak], 'descend');
    switch PICK
        case 'average',  use = j;
        case 'longest',  use = j(1);
        case 'shortest', use = j(end);
        otherwise, error('PICK must be ''average'', ''longest'' or ''shortest''.');
    end
    runs{i} = c(use);
    used = arrayfun(@(r) sprintf('%s (%.0f s)', r.name, r.tPeak), runs{i}, 'UniformOutput', false);
    fprintf('  %6.2f cm        %s\n', z(i), strjoin(used, ', '));
    if numel(use) < numel(c)
        fprintf('                    skipped: %s\n', strjoin({c(setdiff(j, use, 'stable')).name}, ', '));
    end
end
[z, ord] = sort(z);  runs = runs(ord);

% Common time grid up to the LONGEST ramp. At each time, a height's value is
% the mean of the runs still ramping then; after every run at that height
% has been cut it's NaN (grey), instead of truncating every row.
tEnd = max(cellfun(@(r) max([r.tPeak]), runs));
tg = 0:0.1:floor(tEnd);
Umap = nan(numel(z), numel(tg));
nRunsMap = zeros(numel(z), numel(tg));   % how many runs went into each value
for i = 1:numel(z)
    Ui = nan(numel(runs{i}), numel(tg));
    for r = 1:numel(runs{i})
        Ui(r,:) = interp1(runs{i}(r).t, runs{i}(r).Us, tg, 'linear', NaN);
    end
    Umap(i,:) = mean(Ui, 1, 'omitnan');
    nRunsMap(i,:) = sum(~isnan(Ui), 1);
end
Umap(nRunsMap == 0) = NaN;

% Single-hue sequential ramp (light = slow, dark = fast), not a rainbow
rampHex = {'cde2fb','b7d3f6','9ec5f4','86b6ef','6da7ec','5598e7','3987e5', ...
           '2a78d6','256abf','1c5cab','184f95','104281','0d366b'};
ramp = cell2mat(cellfun(@(h) sscanf(h,'%2x%2x%2x',[1 3])/255, rampHex(:), 'UniformOutput', false));
cmap = interp1(linspace(0,1,size(ramp,1)), ramp, linspace(0,1,256));
ink = [0.20 0.20 0.20];  mutedInk = [0.45 0.45 0.45];

%% 1) Time vs height map
f1 = figure('Color','w','Position',[80 80 900 420]);
% pcolor with explicit cell edges, not imagesc: imagesc assumes evenly spaced
% heights, and they aren't (4.85 -> 5.85 is 1 cm, the rest 2 cm), which put
% rows off their tick labels. Each row spans halfway to its neighbours.
if numel(z) < 2, error('Need at least 2 heights for the map.'); end
zE = [z(1) - (z(2)-z(1))/2; (z(1:end-1) + z(2:end))/2; z(end) + (z(end)-z(end-1))/2];
tE = [tg - 0.05, tg(end) + 0.05];
hImg = pcolor(tE, zE, [Umap nan(numel(z),1); nan(1, numel(tg)+1)]);   % pcolor drops the last row/col
hImg.EdgeColor = 'none';
set(gca, 'Color', [0.94 0.94 0.94], 'Layer', 'top');   % NaN (after a run's cut) shows as grey = no data, not "slow"
colormap(f1, cmap);
cb = colorbar; cb.Label.String = 'U (m/s)'; cb.Color = ink;
xlabel('Time since wind start (s)'); ylabel('Height above water (cm)');
yticks(z);
title(sprintf('Streamwise velocity vs time and height (%.0f s moving mean)', SMOOTH_T), ...
    'FontWeight','normal','Color',ink);
set(gca,'TickDir','out','Box','off','XColor',mutedInk,'YColor',mutedInk,'FontSize',10);

%% 2) Profiles at selected times, with an uncertainty band
% Notation: height z, snapshot time t, runs r = 1..n(z,t) (runs still ramping),
% window length Tw = PROFILE_WIN, sampling interval dt = 1/fs.
%
% The point value is the mean over runs of each run's window mean:
%
%   Ubar_r(z,t) = (1/N) sum_j u_r(t_j),     |t_j - t| <= Tw/2,  N = Tw*fs
%   Ubar(z,t)   = (1/n) sum_r Ubar_r
%
% (a) WITHIN-RUN sampling error -- turbulence, not instrument drift.
% The samples inside a window are correlated, so the variance of the window
% mean is NOT sigma_u^2/N. For Tw >> T_int the standard result is
%
%   var(Ubar_r) = (2*T_int/Tw) * sigma_u^2  ==  sigma_u^2 / N_eff,
%   N_eff       = Tw / (2*T_int)
%
% with the fluctuation taken about the LOCAL RAMP, not about the window mean
% (otherwise the ramp's slope is counted as turbulence):
%
%   u'(t_j)   = u(t_j) - (a + b*t_j),   a,b least squares inside the window
%   sigma_u   = std(u')
%   rho(k)    = autocorrelation of u' at lag k*dt          (xcorr, 'coeff')
%   T_int     = dt * sum_{k=0}^{K-1} rho(k),  K = first k with rho(k) < 0
%
% Truncating the sum at the first zero crossing is the usual guard against the
% noisy tail of rho. Measured here: T_int ~ 10-50 ms once the air is turbulent
% (~320 ms before transition, where the signal is still nearly laminar), giving
% N_eff ~ 30-150 and
%
%   s_within(z,t) = (1/n) sum_r sigma_u,r/sqrt(N_eff,r) / sqrt(n)     [eq. A]
%
% i.e. the mean per-run standard error, reduced again by averaging n runs.
% Works out at 0.3-1.2% of U.
%
% (b) RUN-TO-RUN error -- repeatability: fan ramp timing, probe reposition,
% ambient conditions. Per point this is just the sample sd of the run means,
%
%   s_pt(z,t) = sqrt( 1/(n-1) * sum_r (Ubar_r - Ubar)^2 )              [eq. B]
%
% but with n = 2 that has 1 dof: the relative error of an sd estimate is
% ~1/sqrt(2*nu), so nu = 1 means ~70% -- it scatters between 0.3% and 8% here
% purely by luck (see the printed perPoint column). So the scatter is POOLED
% in RELATIVE form over every (z,t) with n >= 2 and t >= UNC_POOL_TMIN, which
% is legitimate because sigma scales with U rather than being constant:
%
%   c = sqrt( sum_p (n_p - 1) * (s_pt,p / Ubar_p)^2 / sum_p (n_p - 1) )  [eq. C]
%   s_run(z,t) = c * Ubar(z,t) / sqrt(n)                                 [eq. D]
%
% Measured: c = 2.22% with sum(n_p - 1) = 91 dof (so c itself is good to ~7%).
% This is 3-5x eq. A -- repeatability, not turbulence, sets the band.
%
% (c) TOTAL -- and the two do NOT add in quadrature. Each run mean Ubar_r
% already carries its own within-run sampling error, so the scatter between
% run means contains it:
%
%   s_pt^2 = sigma_repeat^2 + sigma_u^2/N_eff                          [eq. E]
%
% Eq. D is therefore already the complete standard error of the point mean,
% and the band is simply
%
%   s(z,t) = s_run(z,t) = c * Ubar / sqrt(n)                           [eq. F]
%   band   = Ubar(z,t) +- UNC_K * s(z,t)
%
% Eq. A survives as (i) a diagnostic -- inverting eq. E separates the 2.22%
% into its sampling and genuine-repeatability parts, printed below -- and
% (ii) the fallback for a point with only one run, where no scatter exists.
%
% WHAT THE BAND IS NOT: it is not the turbulent fluctuation. sigma_u itself is
% 4-12% of U here, several times the band, but that is real variation in the
% flow, not doubt about where the mean lies. Shading sigma_u would answer a
% different question ("how much does u swing?"); see plot_midheight_stress.m
% for the turbulence quantities. The band carries 1/sqrt(n), so repeats shrink
% it -- see the note at the end of this file for how many a difference needs.
UNC_K     = 1;        % band half-width in sigma (1 = +-1 sigma, 1.96 = 95%)
UNC_POOL  = true;     % true = pooled run-to-run sd, false = per-point std
UNC_POOL_TMIN = 15;   % s, exclude the pre-transition start from the pooling

nP = numel(PROFILE_T);
Uprof  = nan(numel(z), nP);   % mean U
sWin   = nan(numel(z), nP);   % within-run sampling sigma of the point mean
sRunPt = nan(numel(z), nP);   % per-point run-to-run sd (NaN if only 1 run)
nUsed  = zeros(numel(z), nP);
for k = 1:nP
    for i = 1:numel(z)
        % Mean of each run still ramping at this time (full window inside its
        % ramp), then averaged across runs. NaN if every run was cut already.
        v = nan(numel(runs{i}), 1);  se = nan(numel(runs{i}), 1);
        for r = 1:numel(runs{i})
            if PROFILE_T(k) + PROFILE_WIN/2 <= runs{i}(r).tPeak
                m = abs(runs{i}(r).t - PROFILE_T(k)) <= PROFILE_WIN/2;
                u = runs{i}(r).U(m);  tw = runs{i}(r).t(m);
                v(r) = mean(u);
                % Detrend the ramp inside the window first, else its slope
                % would be counted as turbulence and inflate sigma_u.
                up = u - polyval(polyfit(tw, u, 1), tw);
                fs = 1/median(diff(tw));
                ac = xcorr(up - mean(up), 'coeff');  ac = ac(numel(up):end);
                iz = find(ac < 0, 1);  if isempty(iz), iz = numel(ac); end
                Tint = sum(ac(1:iz-1))/fs;
                se(r) = std(up)/sqrt(max(PROFILE_WIN/(2*Tint), 1));
            end
        end
        g = ~isnan(v);
        if ~any(g), continue; end
        nUsed(i,k) = sum(g);
        Uprof(i,k) = mean(v(g));
        sWin(i,k)  = mean(se(g))/sqrt(sum(g));    % runs average down
        if sum(g) > 1, sRunPt(i,k) = std(v(g)); end
    end
end

% Pooled relative run-to-run sd, from every height x time with >= 2 runs.
w = ~isnan(sRunPt) & repmat(PROFILE_T(:)' >= UNC_POOL_TMIN, numel(z), 1);
relPool = sqrt(sum((nUsed(w)-1) .* (sRunPt(w)./Uprof(w)).^2) / sum(nUsed(w)-1));
if UNC_POOL
    sRun = relPool * Uprof ./ sqrt(nUsed);
else
    sRun = sRunPt ./ sqrt(nUsed);     % NaN where only one run -> band = within-run only
end
% NOT sqrt(sWin^2 + sRun^2): each run mean already carries its own within-run
% sampling error, so the scatter BETWEEN run means already contains it --
% s_pt^2 = sigma_repeat^2 + s_win,1run^2. Adding sWin again double-counts it.
% sRun is the complete standard error of the point mean; sWin is kept only as
% a printed diagnostic and as the fallback where no scatter estimate exists.
sig = sRun;
sig(isnan(sig)) = sWin(isnan(sig));
% How much of the pooled scatter is genuine run-to-run, once the within-run
% sampling noise it contains is removed (small here, but worth seeing).
relWin1 = mean(sWin(w).*sqrt(nUsed(w))./Uprof(w));       % single-run, relative
relRepeat = sqrt(max(relPool^2 - relWin1^2, 0));
fprintf(['Pooled run-to-run sd: %.2f%% of U (%d dof, t >= %g s)\n' ...
         '  of which within-run sampling %.2f%%, genuine repeatability %.2f%%\n'], ...
    100*relPool, sum(nUsed(w)-1), UNC_POOL_TMIN, 100*relWin1, 100*relRepeat);

f2 = figure('Color','w','Position',[100 100 560 480]);
hold on;
lineCols = interp1(linspace(0,1,size(ramp,1)), ramp, linspace(0.25,1,nP));   % skip steps too pale to see
% Bands first so every line sits on top of every band.
for k = 1:nP
    g = ~isnan(Uprof(:,k)) & ~isnan(sig(:,k));
    % One patch per contiguous run of valid heights, so a gap isn't bridged.
    d = [true; diff(find(g)) > 1];
    gi = find(g);  seg = cumsum(d);
    for s = 1:max(seg)
        q = gi(seg == s);
        if numel(q) < 2, continue; end   % a lone point gets an errorbar below
        patch([Uprof(q,k) - UNC_K*sig(q,k); flipud(Uprof(q,k) + UNC_K*sig(q,k))], ...
              [z(q); flipud(z(q))], lineCols(k,:), 'FaceAlpha', 0.28, ...
              'EdgeColor', 'none', 'HandleVisibility', 'off');
    end
end
for k = 1:nP
    plot(Uprof(:,k), z, '-o', 'Color', lineCols(k,:), 'LineWidth', 2, ...
        'MarkerSize', 7, 'MarkerFaceColor', lineCols(k,:), 'MarkerEdgeColor', 'w', ...
        'DisplayName', sprintf('t = %g s', PROFILE_T(k)));
end
hold off;
xlabel('U (m/s)'); ylabel('Height above water (cm)'); yticks(z);
title(['Velocity profiles during the ramp (shading = \pm' num2str(UNC_K) '\sigma)'], ...
    'FontWeight','normal','Color',ink);
legend('Location','eastoutside','Box','off','TextColor',ink);
grid on; set(gca,'GridColor',[0.85 0.85 0.85],'GridAlpha',1,'TickDir','out','Box','off', ...
    'XColor',mutedInk,'YColor',mutedInk,'FontSize',10);

%% Printed numbers behind the band
% nRuns, U, the band half-width, and -- for a sanity check on the pooling --
% the per-point run-to-run sd that a 1-dof estimate would have given.
fprintf('\n   z    t   nRuns      U    +-%gsig   within%%  pooled%%  perPoint%%\n', UNC_K);
for i = 1:numel(z)
    for k = 1:nP
        if isnan(Uprof(i,k)), continue; end
        fprintf('%5.2f %4g %5d %8.3f %8.3f %8.2f %8.2f %9.2f\n', z(i), PROFILE_T(k), ...
            nUsed(i,k), Uprof(i,k), UNC_K*sig(i,k), 100*sWin(i,k)/Uprof(i,k), ...
            100*relPool/sqrt(nUsed(i,k)), ...
            100*sRunPt(i,k)/Uprof(i,k));
    end
end

if ~isempty(OUT_PNG)
    exportgraphics(f1, strrep(OUT_PNG, '.png', '_map.png'), 'Resolution', 150);
    exportgraphics(f2, strrep(OUT_PNG, '.png', '_profiles.png'), 'Resolution', 150);
end

%% How many repeats a given difference needs
% To call two heights different by dU (absolute) at significance alpha and
% power 1-beta, with n runs at EACH height and sd sigma = c*U (eq. C):
%
%   n >= 2 * (z_{1-alpha/2} + z_{1-beta})^2 * sigma^2 / dU^2
%      =  15.7 * (c / delta)^2        for 95% / 80% power, delta = dU/U
%
% With c = 2.22%: resolving delta = 3.4% (the spread across 7.85-15.85 cm at
% t = 25 s) needs n ~ 7 per height; delta = 5% needs n ~ 4; delta = 8.6%
% (t = 65 s) needs n ~ 2, which is why only the late profiles separate
% cleanly in the outer region. The near-surface heights already have 4-5 runs.
