%% recover_interrupted_run.m
% Salvage a run_hotwire_calibration.m run that was stopped part-way (Ctrl-C,
% or an error during the staircase).
%
% RUN THIS AS A SCRIPT, from the same MATLAB session -- it reads hw, fanD and
% the schedule variables straight out of the base workspace. A `clear` between
% the interruption and running this loses the data permanently.
%
% Why this is needed: run_hotwire_calibration.m only calls read(hw,"all")
% AFTER the staircase completes, so an interruption never reaches it. The
% samples are not lost though -- they are still sitting in the DAQ object's
% buffer, and read() pulls them out whether or not the task was stopped.
% Ctrl-C in particular does NOT trigger the script's catch block, so the fan
% is likely still driving and the task still running.
%
% Steps that never completed are dropped: stepEnd is written BEFORE the hold
% is waited out, so the last step will have a boundary even though no data
% covers it. Windows are kept only if the acquired data actually spans them.

fprintf('=== recovering interrupted calibration run ===\n');

%% 1. Fan to zero FIRST -- it may still be driving at the last commanded level
if exist('fanD','var') && isa(fanD,'daq.interfaces.DataAcquisition')
    write(fanD, 0);
    fprintf('Fan set to 0 V.\n');
else
    warning('fanD not found -- CHECK THE FAN IS OFF MANUALLY.');
end

%% 2. Stop acquisition and drain the buffer
if ~exist('hw','var') || ~isa(hw,'daq.interfaces.DataAcquisition')
    error(['hw not found in the workspace -- the acquisition object is gone ' ...
           'and the buffered samples with it. Nothing to recover.']);
end
if hw.Running
    stop(hw);
    fprintf('Acquisition stopped.\n');
else
    fprintf('Acquisition was already stopped (catch block ran).\n');
end

% READ IS DESTRUCTIVE. read(hw,"all") CONSUMES the buffer: a second call
% returns only samples that arrived since the first. Ctrl-C does not stop the
% acquisition, so after an interruption the task keeps buffering idle air --
% and a second read hands back that idle tail instead of the run. Stash any
% existing data before overwriting it, and read into a temporary first.
if exist('data','var') && istimetable(data) && height(data) > 0
    data_prev = data;
    fprintf(['NOTE: a `data` timetable (%d samples) was already in the workspace ' ...
             'and has been copied to `data_prev` before this read.\n'], height(data_prev));
end

newData = read(hw, "all");
if isempty(newData)
    error(['read(hw,"all") returned nothing -- the buffer was already drained by ' ...
           'an earlier read, or acquisition never started. If you read it before, ' ...
           'that earlier `data` IS the run; do not re-read.']);
end

% data.Time is measured from the START OF ACQUISITION, so a first timestamp
% well past zero means someone already read this buffer and we are looking at
% the leftovers, not the run.
tFirst = seconds(newData.Time(1));
if tFirst > 5
    if exist('data_prev','var')
        tail = sprintf(['The earlier read is preserved in `data_prev` (%d samples) -- ' ...
                        'THAT is your run. The tail just read is in `newData`.'], height(data_prev));
    else
        tail = ['No earlier `data` was in the workspace to preserve, so the run ' ...
                'itself is gone unless it was saved or is still held in an open figure.'];
    end
    error('recover:alreadyDrained', ...
        ['This read STARTS at t = %.1f s after acquisition began, so the first %.1f s ' ...
         'were consumed by an EARLIER read(hw,"all"). What came back is the idle tail ' ...
         'buffered after the interruption, not the run. Aborting before anything is ' ...
         'overwritten.\n%s'], tFirst, tFirst, tail);
end

data = newData;
hwStopElapsed = hwStartElapsed + seconds(data.Time(end));

hwT      = seconds(data.Time);
hwT_t0   = hwT + hwStartElapsed;
hwT_wind = hwT_t0 - fanStartElapsed;
fprintf('Recovered %d samples, %.1f s (%.2f min) of data.\n', ...
    height(data), hwT_wind(end)-hwT_wind(1), (hwT_wind(end)-hwT_wind(1))/60);

%% 3. Rebuild the channel vectors
E1 = []; E2 = []; E3 = []; U = []; V = []; W = []; E_ref = []; U_ref = [];
if hasHotwire
    E1 = data{:, strcmp(aiNames,'Probe1')};
    E2 = data{:, strcmp(aiNames,'Probe2')};
    E3 = data{:, strcmp(aiNames,'Probe3')};
    [U, V, W] = convert_E2U_fn(E1, E2, E3, CAL_FILE);
end
if hasRef
    E_ref = data{:, strcmp(aiNames,'RefProbe')};
    [U_ref, cal_ref] = convert_Eref2Uref(E_ref);
    fprintf('Reference probe (cert %s): %.4f..%.4f V -> %.3f..%.3f m/s\n', ...
        cal_ref.id, min(E_ref), max(E_ref), min(U_ref), max(U_ref));
end

%% 4. Keep only steps the data actually covers
% stepEnd(k) is assigned before wait_until() blocks, so the interrupted step
% carries a boundary with no data behind it. Require full coverage.
tEndData = hwT_wind(end);
complete = ~isnan(stepStart) & ~isnan(stepEnd) & stepEnd(:)' <= tEndData;
nDropped = sum(~isnan(stepStart) & ~complete);

avgStart = stepStart + stepSettleT;
avgEnd   = stepEnd;

keep      = find(complete);
stepV     = stepV(keep);
stepDir   = stepDir(keep);
stepStart = stepStart(keep);
stepEnd   = stepEnd(keep);
avgStart  = avgStart(keep);
avgEnd    = avgEnd(keep);
nSteps    = numel(keep);

fprintf('Kept %d complete holds', nSteps);
if nDropped > 0, fprintf(', dropped %d incomplete', nDropped); end
fprintf('.\nLevels recovered (V): %s\n', strjoin(compose('%.2f', stepV), ' '));

%% 5. Per-step averages over the surviving windows
nCh = numel(aiNames);
stepMean = nan(nSteps, nCh); stepStd = nan(nSteps, nCh); stepN = zeros(nSteps,1);
stepUrefMean = nan(nSteps,1); stepUrefStd = nan(nSteps,1); stepUmean = nan(nSteps,1);
allData = data{:,:};
for k = 1:nSteps
    m = hwT_wind >= avgStart(k) & hwT_wind < avgEnd(k);
    stepN(k) = sum(m);
    if stepN(k) == 0, continue; end
    stepMean(k,:) = mean(allData(m,:), 1);
    stepStd(k,:)  = std(allData(m,:), 0, 1);
    if hasRef,     stepUrefMean(k) = mean(U_ref(m)); stepUrefStd(k) = std(U_ref(m)); end
    if hasHotwire, stepUmean(k)    = mean(U(m)); end
end

stepTable = table(stepV(:), string(stepDir(:)), stepStart(:), stepEnd(:), ...
    avgStart(:), avgEnd(:), stepN, ...
    'VariableNames', {'FanV','Direction','tStart','tEnd','tAvgStart','tAvgEnd','nSamples'});
for i = 1:nCh
    stepTable.(['mean_' aiNames{i}]) = stepMean(:,i);
    stepTable.(['std_'  aiNames{i}]) = stepStd(:,i);
end
if hasRef,     stepTable.U_ref = stepUrefMean; stepTable.U_ref_std = stepUrefStd; end
if hasHotwire, stepTable.U_oldcal = stepUmean; end
disp(' '); disp(stepTable);

%% 6. Save
saveName = sprintf('hotwire_cal_%s_partial.mat', datestr(now,'yyyymmdd_HHMMSS'));
saveVars = {'data','hwT','hwT_t0','hwT_wind','CAL_FILE','cal_ref', ...
            'E1','E2','E3','U','V','W','E_ref','U_ref', ...
            'aiConfig','aiNames','aiGroups','sentT','sentV', ...
            'stepV','stepDir','stepStart','stepEnd','avgStart','avgEnd', ...
            'stepMean','stepStd','stepN','stepTable', ...
            'FAN_V_LEVELS','STEP_SETTLE_T','STEP_DWELL_T','HOTWIRE_FS', ...
            'fanStartElapsed','hwStartElapsed','hwStopElapsed'};
saveVars = saveVars(cellfun(@(v) exist(v,'var') == 1, saveVars));   % skip anything absent
save(saveName, saveVars{:});
writetable(stepTable, strrep(saveName,'.mat','_steps.csv'));
fprintf('\nSaved %s (+ _steps.csv)\n', saveName);
