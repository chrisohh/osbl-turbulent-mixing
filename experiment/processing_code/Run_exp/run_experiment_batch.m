function run_experiment_batch()
%% run_experiment_batch.m
% Runs run_experiment.m N_RUNS times at one probe height, with a fixed rest
% between runs, saving each as <SAVE_NAME>_r1.mat, _r2.mat, ... (continuing
% the numbering if some repeats already exist). No prompts during the batch.
%
% Edit the settings below, then press Run (F5).
%
% SAFETY: Ctrl+C at any point (during a run or a rest) sets the fan to 0 V
% automatically -- this runner is a function, so its onCleanup fires on
% Ctrl+C, unlike a plain script (where MATLAB skips catch blocks and you'd
% have to run cleanup_daq.m yourself). If anything still looks stuck after an
% interrupt, run cleanup_daq.m anyway.

%% ---- Settings ------------------------------------------------------------
SAVE_NAME = 'hotwire_midheight-12cm';   % <-- height + cut, as in the plan
N_RUNS    = 2;          % runs in this batch
REST_T    = 480;        % s of rest between the END of one run and the START
                        % of the next -- keep it the same for every run with
                        % the same cut, so each starts from the same water state
PAUSE_AFTER_FIRST = false;  % true = stop after run 1 and WAIT for you to answer
                            % (blocks the batch until you do). Leave false to
                            % run unattended: at a new height, just watch run 1
                            % and press Ctrl+C if the clearance looks unsafe --
                            % the fan is zeroed automatically.
% Settings in run_experiment.m to override for this batch only. Anything not
% listed here uses run_experiment.m's own value.
OVERRIDES = struct( ...
    'FAN_RAMP_STOP_T',     Inf, ...    % must match the stopT in SAVE_NAME
    'PRE_WIND_BASELINE_T', 5);        % 5 s zero-flow drift check each run
DEV_ID = "Dev4";        % fan device, only used to zero the fan on interrupt

%% ---- Checks ----------------------------------------------------------------
% The cut must always be stated here and must match the name, so a file's
% name can't claim a different ramp than it ran (the plot script trusts the
% data, but people read the names).
if ~isfield(OVERRIDES, 'FAN_RAMP_STOP_T')
    error('Put FAN_RAMP_STOP_T in OVERRIDES (Inf = full ramp) so it can be checked against SAVE_NAME.');
end
tok = regexp(SAVE_NAME, 'stopT([\d.]+)s', 'tokens', 'once');
stopT = OVERRIDES.FAN_RAMP_STOP_T;
if isempty(tok) && isfinite(stopT)
    error('FAN_RAMP_STOP_T = %g but SAVE_NAME has no _stopT%gs -- add it, or use Inf for a full ramp.', stopT, stopT);
elseif ~isempty(tok) && str2double(tok{1}) ~= stopT
    error('SAVE_NAME says stopT%ss but FAN_RAMP_STOP_T = %g -- make them match.', tok{1}, stopT);
end

% A manual run of run_experiment.m leaves hw/fanD/camD in the BASE workspace,
% still holding Dev4 -- the batch then can't get the hardware ("Hardware is
% reserved"). Zero the fan through them and release them first.
release_base_daq();

% Fan to 0 V however this function ends: normally, on error, or on Ctrl+C.
fanGuard = onCleanup(@() zero_fan(DEV_ID)); %#ok<NASGU>

tBatch = tic;
fprintf('\n===== Batch: %d x %s, rest %d s =====\n', N_RUNS, SAVE_NAME, REST_T);
for n = 1:N_RUNS
    fprintf('\n----- Run %d of %d (%s) -----\n', n, N_RUNS, datestr(now, 'HH:MM:SS')); %#ok<TNOW1,DATST>
    figsBefore = findall(groot, 'Type', 'figure');
    do_one_run(struct('name', SAVE_NAME, 'overrides', OVERRIDES));
    % Keep only the last run's figures open, so they don't pile up.
    if n < N_RUNS
        close(setdiff(findall(groot, 'Type', 'figure'), figsBefore));
    end

    if n == N_RUNS, break; end

    if n == 1 && PAUSE_AFTER_FIRST
        r = input('Run 1 done. Continue with the rest of the batch? [Y/n]: ', 's');
        if ~(isempty(r) || strncmpi(strtrim(r), 'y', 1))
            disp('Batch stopped after run 1.');
            return
        end
    end

    % Rest: fan is already at 0 (the run ramps down to 0 at the end).
    fprintf('Resting %d s before run %d...', REST_T, n+1);
    tRest = tic;
    while toc(tRest) < REST_T
        pause(min(1, REST_T - toc(tRest)));
        left = REST_T - toc(tRest);
        if mod(round(left), 60) == 0 && left > 1
            fprintf(' %d s', round(left));
        end
    end
    fprintf(' done.\n');
end
fprintf('\n===== Batch finished: %d runs in %.1f min =====\n', n, toc(tBatch)/60);
end


function do_one_run(BATCH_RUN) %#ok<INUSD> -- read by run_experiment.m
% Runs the script in THIS function's workspace, so its clearvars can't touch
% the batch loop's variables.
run_experiment;
end


function zero_fan(DEV_ID)
% Two tries: a fresh session; if the output is still reserved by a leftover
% object, release everything (daqreset) and try again.
for attempt = 1:2
    try
        d = daq("ni");
        addoutput(d, DEV_ID, "ao0", "Voltage");
        write(d, 0);
        disp('Fan set to 0 V.');
        return
    catch ME
        if attempt == 1
            release_base_daq();
            daqreset;
        else
            fprintf(2, 'Could not zero the fan (%s). CHECK THE FAN PHYSICALLY.\n', ME.message);
        end
    end
end
end


function release_base_daq()
% Stop and clear DAQ objects a manual run left in the base workspace, writing
% 0 V through any fan object first (it's the one that can still drive ao0).
for nm = ["fanD", "fanCleanup", "hw", "camD"]
    if ~evalin('base', sprintf('exist(''%s'', ''var'')', nm)), continue; end
    obj = evalin('base', nm);
    if startsWith(nm, "fan")
        try, write(obj, 0); catch, end
    end
    try, stop(obj); catch, end
    evalin('base', sprintf('clear %s', nm));
    fprintf('Released leftover %s from the base workspace.\n', nm);
end
end
