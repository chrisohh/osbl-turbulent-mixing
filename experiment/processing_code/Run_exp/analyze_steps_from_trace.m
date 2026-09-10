%% analyze_steps_from_trace.m
% Re-cut a staircase run's per-step averages at a chosen settle time, and plot
% velocity vs fan voltage from the SETTLED portion of each hold only.
%
% Works from a recovered trace (t_rec/E_rec pulled out of a figure) or from a
% run still in the workspace (hwT_wind/E_ref). The point of re-cutting offline
% is that STEP_SETTLE_T only ever controlled which samples got AVERAGED -- the
% raw trace always held the whole hold, so the settle choice can be revisited
% afterwards without touching the tunnel.

SETTLE_T = 20;    % <-- s discarded at the start of each hold before averaging (18 s measured)
SWEEP    = [0 5 10 15 20 30 45 60];   % settle times compared in the table

CTA_DIR = 'D:\Chris\osbl-turbulent-mixing\experiment\processing_code\CTA';
addpath(CTA_DIR);

%% 1. Locate the trace -------------------------------------------------------
if exist('t_rec','var') && exist('E_rec','var')
    t = t_rec(:);  y = E_rec(:);
    fprintf('Using t_rec / E_rec from the workspace (%d samples).\n', numel(t));
elseif exist('hwT_wind','var') && exist('E_ref','var') && ~isempty(E_ref)
    t = hwT_wind(:); y = E_ref(:);
    fprintf('Using hwT_wind / E_ref from the workspace (%d samples).\n', numel(t));
else
    % Fall back to the longest line in any open figure.
    h = findobj(0,'Type','line');
    if isempty(h), error('No trace found: need t_rec/E_rec, hwT_wind/E_ref, or an open figure.'); end
    [~, ix] = max(arrayfun(@(x) numel(x.XData), h));
    t = h(ix).XData(:); y = h(ix).YData(:);
    fprintf('Pulled the longest line out of an open figure (%d samples).\n', numel(t));
end

% Is this voltage or already-converted velocity? The 54T29 tops out at 4.83 V
% on the certificate, so a trace peaking well above that is m/s, not volts.
if max(y) > 5.2
    U = y;
    fprintf('Trace looks like VELOCITY (max %.2f) -- using as-is.\n', max(y));
else
    [U, cal] = convert_Eref2Uref(y);
    fprintf('Trace looks like VOLTAGE (max %.3f V) -- converted via cert %s.\n', ...
        max(y), cal.id);
end

%% 2. Rebuild the step schedule ---------------------------------------------
% Taken from the run's own parameters when they are still in the workspace,
% otherwise the values used for the 0->9.5 V run.
if exist('FAN_V_LEVELS','var'), lv = FAN_V_LEVELS(:)'; else, lv = [0 2.1 4 6 8 8.5 9 9.5]; end
if exist('PRE_RUN_ZERO_T','var'), pz = PRE_RUN_ZERO_T; else, pz = 10;  end
if exist('STEP_DWELL_T','var'),   dw = STEP_DWELL_T;   else, dw = 90;  end

stepV = [0, lv];
bnd   = [0, pz, pz + dw*(1:numel(lv))];
sS    = bnd(1:end-1);
sE    = bnd(2:end);

% A hold cut short by an interruption is still usable as long as enough
% SETTLED data survives in it -- the settle is at the START of the hold, so a
% truncated tail loses averaging time, not the part that was contaminated.
% Clip such a hold to where the data actually ends and keep it if at least
% MIN_DWELL seconds remain past the settle; drop it otherwise.
MIN_DWELL = 10;   % s of settled data required to keep a truncated hold

sE = min(sE, t(end));
usable = (sE - sS - SETTLE_T) >= MIN_DWELL;
truncated = sE < bnd(2:end) & usable;
truncated = truncated(1:numel(usable));

if any(truncated)
    fprintf('Keeping %d truncated hold(s) (levels %s V) -- %s s of settled data left.\n', ...
        sum(truncated), strjoin(compose('%.2f', stepV(truncated)), ' '), ...
        strjoin(compose('%.0f', sE(truncated)-sS(truncated)-SETTLE_T), ' '));
end
if any(~usable)
    fprintf('Dropping %d hold(s) with <%d s of settled data (levels %s V).\n', ...
        sum(~usable), MIN_DWELL, strjoin(compose('%.2f', stepV(~usable)), ' '));
end
stepV = stepV(usable); sS = sS(usable); sE = sE(usable);
n = numel(stepV);

%% 3. Settle sweep -----------------------------------------------------------
M = nan(n, numel(SWEEP));  S = nan(n, numel(SWEEP));
for j = 1:numel(SWEEP)
    for k = 1:n
        m = t >= sS(k)+SWEEP(j) & t < sE(k);
        if sum(m) > 10, M(k,j) = mean(U(m)); S(k,j) = std(U(m)); end
    end
end

fprintf('\nMean U_ref (m/s) by settle time discarded:\n%6s', 'FanV');
% Built one label at a time: compose() returns a cell array here, which
% fprintf's %s cannot consume (strjoin can, which is why the other compose
% calls in this file are fine).
for j = 1:numel(SWEEP)
    fprintf('%9s', sprintf('%ds', SWEEP(j)));
end
fprintf('%10s\n', 'drift%');
for k = 1:n
    drift = 100*(M(k,end)-M(k,1))/M(k,1);
    fprintf('%6.2f', stepV(k)); fprintf('%9.4f', M(k,:)); fprintf('%10.2f\n', drift);
end

% Settled = the point past which extending the settle changes the mean by less
% than the sample-to-sample scatter of the mean itself.
fprintf('\nSuggested settle: the column where drift stops changing materially.\n');

%% 4. Averages at the chosen SETTLE_T ---------------------------------------
Um = nan(n,1); Us = nan(n,1); Un = zeros(n,1);
for k = 1:n
    m = t >= sS(k)+SETTLE_T & t < sE(k);
    Un(k) = sum(m);
    if Un(k) > 10, Um(k) = mean(U(m)); Us(k) = std(U(m)); end
end
T = table(stepV(:), sS(:), sE(:), Un, Um, Us, Us./Um, ...
    'VariableNames', {'FanV','tStart','tEnd','nSamples','U_mean','U_std','TurbIntensity'});
fprintf('\nSettled averages (discarding %d s per hold):\n', SETTLE_T); disp(T);

%% 5. Plots ------------------------------------------------------------------
fig = figure('Position',[80 80 1150 850],'Color','w');
tiledlayout(2,2,'TileSpacing','compact','Padding','compact');

ax = nexttile([1 2]);
plot(t, U, 'LineWidth', 0.5, 'Color', [0.2 0.55 0.25]); hold on;
hb = xregion(ax, (sS+SETTLE_T)', sE', 'FaceColor',[0.1 0.3 0.7], 'FaceAlpha',0.12);
set(hb,'HandleVisibility','off'); hold off;
xlabel('Time since first fan write (s)'); ylabel('U_{ref} (m/s)');
title(sprintf('Recovered trace -- shaded = averaged (settle %d s discarded)', SETTLE_T));
grid on; xlim([t(1) t(end)]);

nexttile; hold on;
% Every hold on a common "time since step change" axis: this is the settling
% time read directly off the data rather than guessed.
cols = parula(n);
for k = 1:n
    m = t >= sS(k) & t < sE(k);
    plot(t(m)-sS(k), U(m), 'Color', [cols(k,:) 0.85], 'LineWidth', 0.6, ...
        'DisplayName', sprintf('%.1f V', stepV(k)));
end
yl = ylim; plot([SETTLE_T SETTLE_T], yl, 'r--', 'LineWidth', 1.5, 'DisplayName','settle cut');
ylim(yl); hold off;
xlabel('Time since step change (s)'); ylabel('U_{ref} (m/s)');
title('Settling, all holds overlaid'); legend('Location','eastoutside','FontSize',7); grid on;

nexttile;
errorbar(stepV, Um, Us, 'o-', 'LineWidth', 1.6, 'MarkerSize', 7, ...
    'MarkerFaceColor', 'w', 'CapSize', 5);
xlabel('Fan command (V)'); ylabel('U_{ref} (m/s)');
title(sprintf('Velocity vs fan voltage (settled, %d s discarded)', SETTLE_T));
grid on;

fprintf('\nDone. Settled points are in table T.\n');
