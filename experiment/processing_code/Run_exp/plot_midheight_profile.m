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

%% 2) Profiles at selected times
f2 = figure('Color','w','Position',[100 100 560 480]);
hold on;
nP = numel(PROFILE_T);
lineCols = interp1(linspace(0,1,size(ramp,1)), ramp, linspace(0.25,1,nP));   % skip steps too pale to see
Uprof = nan(numel(z), nP);
for k = 1:nP
    for i = 1:numel(z)
        % Mean of each run still ramping at this time (full window inside its
        % ramp), then averaged across runs. NaN if every run was cut already.
        v = nan(numel(runs{i}), 1);
        for r = 1:numel(runs{i})
            if PROFILE_T(k) + PROFILE_WIN/2 <= runs{i}(r).tPeak
                m = abs(runs{i}(r).t - PROFILE_T(k)) <= PROFILE_WIN/2;
                v(r) = mean(runs{i}(r).U(m));
            end
        end
        Uprof(i,k) = mean(v, 'omitnan');   % stays NaN -> point left out
    end
    plot(Uprof(:,k), z, '-o', 'Color', lineCols(k,:), 'LineWidth', 2, ...
        'MarkerSize', 7, 'MarkerFaceColor', lineCols(k,:), 'MarkerEdgeColor', 'w', ...
        'DisplayName', sprintf('t = %g s', PROFILE_T(k)));
end
hold off;
xlabel('U (m/s)'); ylabel('Height above water (cm)'); yticks(z);
title('Velocity profiles during the ramp','FontWeight','normal','Color',ink);
legend('Location','eastoutside','Box','off','TextColor',ink);
grid on; set(gca,'GridColor',[0.85 0.85 0.85],'GridAlpha',1,'TickDir','out','Box','off', ...
    'XColor',mutedInk,'YColor',mutedInk,'FontSize',10);

disp(array2table([z Uprof], 'VariableNames', ['z_above_water_cm', compose('t%gs', PROFILE_T)]));

if ~isempty(OUT_PNG)
    exportgraphics(f1, strrep(OUT_PNG, '.png', '_map.png'), 'Resolution', 150);
    exportgraphics(f2, strrep(OUT_PNG, '.png', '_profiles.png'), 'Resolution', 150);
end
