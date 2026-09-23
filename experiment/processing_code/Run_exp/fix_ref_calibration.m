%% fix_ref_calibration.m
% Recomputes U_ref (and the per-step table) for a calibration .mat saved
% BEFORE the 2026-09-23 fix, when run_hotwire_calibration.m /
% run_fan_kick_calibration.m used convert_Eref2Uref's built-in factory
% certificate instead of CAL_FILE's T29 block (the source run_experiment.m
% actually uses). No re-run needed: E_ref, avgStart/avgEnd, stepV/stepDir
% were all saved, so U_ref is just recomputed from them.

IN_FILE  = 'hotwire_cal_20260923_012711.mat';   % <-- edit per file
OUT_FILE = strrep(IN_FILE, '.mat', '_fixed.mat');

CTA_DIR = 'D:\Chris\osbl-turbulent-mixing\experiment\processing_code\CTA';
addpath(CTA_DIR);

S = load(IN_FILE);
if ~isfield(S, 'E_ref') || isempty(S.E_ref)
    error('%s has no E_ref -- reference probe was not logged in this run.', IN_FILE);
end

calRef = parse_probe_section(S.CAL_FILE, 'T29');
U_ref  = polyval(fliplr(calRef.C(1,:)), S.E_ref);
belowRef = S.E_ref < calRef.E_floor(1);
U_ref(belowRef) = (calRef.U_floor(1) / calRef.E_floor(1)) * S.E_ref(belowRef);
aboveRef = S.E_ref > calRef.E_ceil(1);
U_ref(aboveRef) = calRef.U_ceil(1);
fprintf('Recomputed U_ref (T29, %s): %.3f..%.3f m/s (old range was %.3f..%.3f)\n', ...
    S.CAL_FILE, min(U_ref), max(U_ref), min(S.U_ref), max(S.U_ref));

%% Re-average per step using the saved windows (unchanged fan schedule)
nSteps = numel(S.stepV);
stepUrefMean = nan(nSteps, 1);
stepUrefStd  = nan(nSteps, 1);
for k = 1:nSteps
    m = S.hwT_wind >= S.avgStart(k) & S.hwT_wind < S.avgEnd(k);
    if any(m)
        stepUrefMean(k) = mean(U_ref(m));
        stepUrefStd(k)  = std(U_ref(m));
    end
end

stepTable = S.stepTable;
stepTable.U_ref     = stepUrefMean;
stepTable.U_ref_std = stepUrefStd;
disp(stepTable);

%% Save alongside the original -- nothing overwritten
U_ref_old = S.U_ref;   % kept for comparison, not used downstream
save(OUT_FILE, '-struct', 'S');
save(OUT_FILE, 'U_ref', 'stepTable', 'calRef', 'U_ref_old', '-append');
writetable(stepTable, strrep(OUT_FILE, '.mat', '_steps.csv'));
fprintf('Saved %s and its _steps.csv\n', OUT_FILE);
