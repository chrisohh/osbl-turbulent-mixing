function plot_experiment_signals(mode, varargin)
% PLOT_EXPERIMENT_SIGNALS  Plot fan, hot-wire, and camera trigger signals
% in stacked subplots -- works for either planned (pre-run) signals or
% actual logged data (post-run).
%
% USAGE 1 -- Planned/expected signals (no hardware needed):
%   plot_experiment_signals('planned', ...
%       'FanVStart', 0, 'FanVEnd', 7, ...
%       'RampUpTime', 30, 'HoldTime', 5, 'RampDownTime', 5, 'Dt', 0.5, ...
%       'DelayBeforeTrigger', 30, 'TriggerWidth', 0.01, ...
%       'CamFreq', 50, 'CamDuration', 5, ...
%       'HotwireFs', 4000, 'HotwireTotalTime', 65)
%
% USAGE 2 -- Actual logged data (after a real run):
%   plot_experiment_signals('actual', ...
%       'FanT', sentT, 'FanV', sentV, ...
%       'HotwireT', hwTime, 'HotwireV', hwData, ...
%       'TrigT', trigTime, 'TrigV', trigLogic)
%
% NOTE ON CAMERA TRIGGER SUBPLOT: this shows IR/color/mono as three
% overlaid or stacked pulse trains at 50 Hz starting at DelayBeforeTrigger.
% This assumes all three cameras trigger simultaneously at the same 50 Hz
% -- if IR/color/mono have different phase offsets or the pulse durations
% differ (e.g. your 700 microsec color exposure vs IR), this plot does NOT
% capture that detail yet. Let me know if you need per-camera pulse width
% shown accurately; right now they're drawn as simple 50% duty square
% waves for visualization only.

    p = inputParser;
    % Planned-mode params
    addParameter(p, 'FanVStart', 0);
    addParameter(p, 'FanVEnd', 7);
    addParameter(p, 'RampUpTime', 30);
    addParameter(p, 'HoldTime', 5);
    addParameter(p, 'RampDownTime', 5);
    addParameter(p, 'Dt', 0.5);
    addParameter(p, 'DelayBeforeTrigger', 30);
    addParameter(p, 'TriggerWidth', 0.01);
    addParameter(p, 'CamFreq', 50);
    addParameter(p, 'CamDuration', 5);
    addParameter(p, 'HotwireFs', 4000);
    addParameter(p, 'HotwireTotalTime', 65);
    % Actual-mode params
    addParameter(p, 'FanT', []);
    addParameter(p, 'FanV', []);
    addParameter(p, 'HotwireT', []);
    addParameter(p, 'HotwireV', []);
    addParameter(p, 'TrigT', []);
    addParameter(p, 'TrigV', []);
    parse(p, varargin{:});
    r = p.Results;

    figure('Name', ['Experiment signals -- ' mode], 'Position', [100 100 900 700]);

    switch lower(mode)
        case 'planned'
            % ---- Fan voltage profile ----
            tUp   = 0 : r.Dt : r.RampUpTime;
            vUp   = linspace(r.FanVStart, r.FanVEnd, numel(tUp));
            tHold = tUp(end) + r.Dt : r.Dt : tUp(end) + r.HoldTime;
            vHold = r.FanVEnd * ones(size(tHold));
            tDown = tHold(end) + r.Dt : r.Dt : tHold(end) + r.RampDownTime;
            vDown = linspace(r.FanVEnd, r.FanVStart, numel(tDown));
            fanT = [tUp, tHold, tDown];
            fanV = [vUp, vHold, vDown];

            % ---- Hot-wire placeholder (planned = flat "acquiring" indicator) ----
            hwT = 0 : 1/r.HotwireFs : r.HotwireTotalTime;
            hwV = zeros(size(hwT)); % planned mode has no real signal, just shows acquisition window

            % ---- Camera trigger pulses (simplified 50% duty square wave) ----
            trigStart = r.DelayBeforeTrigger;
            trigEnd = trigStart + r.CamDuration;
            dtTrig = 1/(r.CamFreq*200); % fine resolution for visualization
            trigT = 0 : dtTrig : max(fanT(end), hwT(end));
            period = 1/r.CamFreq;
            inWindow = trigT >= trigStart & trigT <= trigEnd;
            trigV = double(mod(trigT, period) < (period/2)) .* inWindow;

        case 'actual'
            fanT = r.FanT; fanV = r.FanV;
            hwT  = r.HotwireT; hwV = r.HotwireV;
            trigT = r.TrigT; trigV = r.TrigV;
            if isempty(fanT) || isempty(hwT) || isempty(trigT)
                error('Actual mode requires FanT/FanV, HotwireT/HotwireV, and TrigT/TrigV to be provided.');
            end

        otherwise
            error('mode must be ''planned'' or ''actual''');
    end

    subplot(3,1,1);
    plot(fanT, fanV, 'LineWidth', 1.5, 'Color', [0.85 0.33 0.10]);
    ylabel('Fan (V)');
    title(['Fan voltage -- ' mode]);
    grid on;

    subplot(3,1,2);
    plot(hwT, hwV, 'LineWidth', 0.75, 'Color', [0.0 0.45 0.74]);
    ylabel('Hot-wire (V)');
    title(['Hot-wire signal -- ' mode]);
    grid on;

    subplot(3,1,3);
    plot(trigT, trigV, 'LineWidth', 1, 'Color', [0.47 0.67 0.19]);
    ylabel('Trigger (logic)');
    xlabel('Time (s)');
    title(['Camera trigger pulses -- ' mode]);
    ylim([-0.2 1.2]);
    grid on;

    sgtitle('Experiment signal overview');
end
