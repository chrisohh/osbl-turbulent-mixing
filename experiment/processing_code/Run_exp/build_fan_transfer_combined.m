%% build_fan_transfer_combined.m
% Merges every calibration run's up-leg points at or above 2.1V into one
% stepTable, then runs fit_fan_transfer.m on it.
%
% Per-run coverage (FAN_V_LEVELS used in each):
%   hotwire_cal_20260922_234655  [2 2.1 2.2 2.4 2.7 3 3.5 4 5 6 7 8]     -- main 2-8V sweep
%   hotwire_cal_20260923_010212  [2.01 2.1 2.2 2.3 2.4]                  -- near-onset
%   hotwire_cal_20260923_012711  [2.01 2.1]                              -- near-onset
%   hotwire_cal_20260923_013607  [2.1 2.05 2.01]                         -- near-onset
%   hotwire_cal_20260923_034304  [8 8.5 9 9.5]                           -- high-V extension
% Everything below 2.1V is dropped (fan cycles on/off at 2.00V, see memory).
% That leaves real overlap at 2.1V (4 runs) and 8V (2 runs) -- averaged
% here (weighted by nSamples) rather than picking one run arbitrarily.
%
% All five files already carry the corrected (CAL_FILE/T29) U_ref -- the
% four that predated that fix were corrected in place and the "_fixed"
% suffix dropped; 034304 was run after the fix and needed no correction.

FILES = { ...
    'hotwire_cal_20260922_234655.mat', ...
    'hotwire_cal_20260923_010212.mat', ...
    'hotwire_cal_20260923_012711.mat', ...
    'hotwire_cal_20260923_013607.mat', ...
    'hotwire_cal_20260923_034304.mat' };

V_MIN = 2.1;   % fan cycles on/off at 2.00V -- exclude everything below this

rows = table();
for i = 1:numel(FILES)
    S = load(FILES{i}, 'stepTable');
    T = S.stepTable;
    keep = T.Direction == "up" & T.FanV >= V_MIN & ~isnan(T.U_ref);
    Ti = T(keep, {'FanV','U_ref','nSamples'});
    Ti.source = repmat(string(FILES{i}), height(Ti), 1);
    rows = [rows; Ti]; %#ok<AGROW>
    fprintf('%-40s %d usable point(s) >= %.2fV\n', FILES{i}, height(Ti), V_MIN);
end

%% Average duplicate voltages (nSamples-weighted), keep everything else
[uV, ~, ic] = unique(rows.FanV);
Um = nan(numel(uV), 1);
Nsrc = zeros(numel(uV), 1);
for k = 1:numel(uV)
    m = ic == k;
    w = rows.nSamples(m);
    if all(w == 0), w = ones(sum(m),1); end   % fall back to plain mean if nSamples missing
    Um(k) = sum(rows.U_ref(m) .* w) / sum(w);
    Nsrc(k) = sum(m);
end

fprintf('\n%d unique voltages after merging (%d raw points):\n', numel(uV), height(rows));
disp(table(uV, Um, Nsrc, 'VariableNames', {'FanV','U_ref','nRunsAveraged'}));

%% Hand off to fit_fan_transfer.m in its expected stepTable format
stepTable = table(uV(:), repmat("up", numel(uV), 1), Um(:), ...
    'VariableNames', {'FanV','Direction','U_ref'});
clear T Ti S rows uV Um Nsrc ic k m w i   % keep only what fit_fan_transfer.m looks for
run('fit_fan_transfer.m');
