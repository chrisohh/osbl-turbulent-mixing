%% fit_fan_transfer.m
% Steady-state fan transfer function: fan command V -> settled U_ref.
% Saved as a monotone lookup (pchip) so run_experiment.m can INVERT it and
% build a ramp that is linear in velocity instead of linear in voltage.
%
% Input (whichever is in the workspace, first match wins):
%   T          -- from analyze_steps_from_trace.m (FanV, U_mean)
%   stepTable  -- from run_hotwire_calibration.m, ref probe logged (FanV, U_ref)
%
% Only the fan-ON part is fitted: below FAN_V_ONSET the fan does not turn, so
% V -> U is not invertible there (every V gives ~0). The ramp starts AT the
% onset, which means the flow jumps 0 -> U(onset) at wind start; that step is
% a property of the fan, not something a transfer function can remove.
%
% Does NOT load a file itself -- it works off T or stepTable already in the
% workspace, so it can be called standalone after a manual load(), or via
% build_fan_transfer_combined.m after that script builds a merged stepTable.
FAN_V_ONSET = 2.1;   % <-- V at which the fan first turns
OUT_FILE    = 'fan_transfer.mat';

%% 1. Pull the settled (V, U) points ----------------------------------------
if exist('T','var') && all(ismember({'FanV','U_mean'}, T.Properties.VariableNames))
    Vp = T.FanV;  Up = T.U_mean;
    src = 'T (analyze_steps_from_trace)';
elseif exist('stepTable','var') && ismember('U_ref', stepTable.Properties.VariableNames)
    up = stepTable.Direction == "up";      % up-leg only, see calibration-run notes
    Vp = stepTable.FanV(up);  Up = stepTable.U_ref(up);
    src = 'stepTable (run_hotwire_calibration)';
else
    error('Need T (analyze_steps_from_trace) or stepTable with U_ref (run_hotwire_calibration).');
end

keep = ~isnan(Up) & Vp >= FAN_V_ONSET;
Vp = Vp(keep);  Up = Up(keep);
[Vp, ix] = sort(Vp);  Up = Up(ix);

if numel(Vp) < 3
    error('Only %d usable points at or above %.2f V -- need at least 3.', numel(Vp), FAN_V_ONSET);
end
% interp1 in the inverse direction needs U strictly increasing with V. A
% non-monotone point means a hold that did not settle or a bad level --
% fix the data rather than papering over it.
if any(diff(Up) <= 0)
    disp(table(Vp, Up));
    error('U is not strictly increasing with V -- check the holds above before fitting.');
end

fprintf('Fan transfer from %s: %d points, %.2f..%.2f V -> %.3f..%.3f m/s\n', ...
    src, numel(Vp), Vp(1), Vp(end), Up(1), Up(end));

%% 2. Save ------------------------------------------------------------------
fanTF.V = Vp;  fanTF.U = Up;  fanTF.V_onset = FAN_V_ONSET;
fanTF.source = src;  fanTF.created = datestr(now, 'yyyy-mm-dd HH:MM:SS');
save(OUT_FILE, 'fanTF');
fprintf('Saved %s\n', OUT_FILE);

%% 3. Plot: forward map, and what a velocity-linear ramp looks like -----------
Vf = linspace(Vp(1), Vp(end), 300);
Ut = linspace(Up(1), Up(end), 121);
Vt = interp1(Up, Vp, Ut, 'pchip');

figure('Color','w','Position',[100 100 1000 400]);
tiledlayout(1,2,'TileSpacing','compact');
nexttile;
plot(Vp, Up, 'o', 'MarkerFaceColor','w', 'LineWidth',1.4); hold on;
plot(Vf, interp1(Vp, Up, Vf, 'pchip'), 'k-'); hold off;
xlabel('Fan command (V)'); ylabel('U_{ref} (m/s)'); grid on;
title('Steady-state transfer V \rightarrow U');

nexttile;
tt = linspace(0, 1, numel(Ut));
plot(tt, Vt, 'LineWidth',1.4); hold on;
plot(tt, linspace(Vp(1), Vp(end), numel(Ut)), '--'); hold off;
xlabel('Fraction of ramp-up'); ylabel('Fan command (V)'); grid on;
legend('linear in U (new)', 'linear in V (old)', 'Location','northwest');
title('Voltage schedule for a linear velocity ramp');
