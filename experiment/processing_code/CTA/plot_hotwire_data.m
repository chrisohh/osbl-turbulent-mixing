function plot_hotwire_data(matFile)
% PLOT_HOTWIRE_DATA  Plot hot-wire data saved by run_experiment.m
%
%   plot_hotwire_data()          - prompts to select a hotwire_*.mat file
%   plot_hotwire_data(matFile)   - loads the given .mat file directly
%
% Expects the variables saved by run_experiment.m's save() call:
%   hwT_wind, E1, E2, E3, U, V, W, E_ref, aiGroups, sentT, sentV,
%   fanStartElapsed, hwStartElapsed, hwStopElapsed
% Older/partial files (e.g. hot-wire or ref-probe disabled) are handled --
% missing channels are simply skipped.

    if nargin < 1 || isempty(matFile)
        defaultDir = fullfile(fileparts(mfilename('fullpath')), '..', 'Run_exp');
        [fname, fpath] = uigetfile(fullfile(defaultDir, 'hotwire_*.mat'), ...
            'Select hot-wire data file');
        if isequal(fname, 0)
            disp('No file selected.');
            return;
        end
        matFile = fullfile(fpath, fname);
    end

    % The reference-probe calibration is the certificate table built into
    % convert_Eref2Uref.m (54T29 s/n 0202, cert T29-202) -- no constants are
    % passed in from here any more.

    S = load(matFile);
    fprintf('Loaded %s\n', matFile);

    hasHotwire = isfield(S, 'U') && ~isempty(S.U);
    hasRef     = isfield(S, 'E_ref') && ~isempty(S.E_ref);
    hasFan     = isfield(S, 'sentT') && ~isempty(S.sentT);

    U_ref = [];
    if hasRef
        [U_ref, cal_ref] = convert_Eref2Uref(S.E_ref);
        fprintf('Reference probe: %s s/n %s, certificate %s (%s)\n', ...
            cal_ref.probe, cal_ref.serial, cal_ref.id, cal_ref.date);
    end

    if isfield(S, 'hwT_wind')
        t = S.hwT_wind;
    elseif isfield(S, 'hwT')
        t = S.hwT;
    else
        error('plot_hotwire_data:noTime', 'No hwT_wind/hwT time vector found in %s.', matFile);
    end

    fanT = [];
    if hasFan
        fanT = S.sentT;
        if isfield(S, 'fanStartElapsed')
            fanT = fanT - S.fanStartElapsed;
        end
    end

    %% Fan ramp (context for interpreting the hot-wire trace)
    if hasFan
        figure('Name', 'Fan voltage', 'Position', [80 80 800 300]);
        plot(fanT, S.sentV, 'LineWidth', 1.5, 'Color', [0.85 0.33 0.10]);
        xlabel('Time since wind start (s)');
        ylabel('Fan (V)');
        title('Fan voltage ramp');
        grid on;
    end

    %% Raw voltages
    nSub = double(hasHotwire) + double(hasRef);
    if nSub > 0
        figure('Name', 'Hot-wire raw voltages', 'Position', [80 420 800 300*nSub]);
        iSub = 0;

        if hasHotwire
            iSub = iSub + 1;
            subplot(nSub, 1, iSub);
            plot(t, S.E1, t, S.E2, t, S.E3, 'LineWidth', 1);
            xlabel('Time since wind start (s)');
            ylabel('Voltage (V)');
            legend('E1', 'E2', 'E3');
            title('Hot-wire raw voltages');
            grid on;
        end

        if hasRef
            iSub = iSub + 1;
            subplot(nSub, 1, iSub);
            plot(t, S.E_ref, 'LineWidth', 1);
            xlabel('Time since wind start (s)');
            ylabel('Voltage (V)');
            legend('E_{ref}');
            title('Reference probe (raw voltage)');
            grid on;
        end
    end

    %% Reference-probe velocity (converted via the 54T29 calibration)
    if hasRef
        figure('Name', 'Reference probe velocity', 'Position', [900 420 800 300]);
        plot(t, U_ref, 'LineWidth', 1, 'Color', [0.49 0.18 0.56]);
        xlabel('Time since wind start (s)');
        ylabel('U_{ref} (m/s)');
        title('Reference probe velocity (54T29)');
        grid on;

        onIdx = t >= 0;
        fprintf('\n=== Reference probe velocity (t >= 0) ===\n');
        fprintf('Mean U_ref = %.3f m/s, RMS = %.3f m/s\n', ...
            mean(U_ref(onIdx)), std(U_ref(onIdx)));
    end

    %% Converted velocity
    if hasHotwire
        figure('Name', 'Converted velocity', 'Position', [900 80 800 600]);

        % subplot(2, 1, 1);
        plot(t, S.U, t, S.V, t, S.W, 'LineWidth', 1);
        xlabel('Time since wind start (s)');
        ylabel('Velocity (m/s)');
        legend('U', 'V', 'W');
        title('Converted velocity (probe coordinates)');
        grid on;


        % if hasRef
        %     subplot(2, 1, 2);
        %     hold on;
        %     plot(t, U_ref, 'LineWidth', 1, 'Color', [0.49 0.18 0.56]);
        %     legend('|U| (hot-wire)', 'U_{ref} (54T29)');
        %     hold off;
        % 
        % xlabel('Time since wind start (s)');
        % ylabel('Velocity (m/s)');
        % title('Velocity magnitude');
        % end
        grid on;
    end

    % %% Summary statistics over the wind-on portion (t >= 0)
    % if hasHotwire
    %     onIdx = t >= 0;
    %     Um = mean(S.U(onIdx)); Vm = mean(S.V(onIdx)); Wm = mean(S.W(onIdx));
    %     Urms = std(S.U(onIdx)); Vrms = std(S.V(onIdx)); Wrms = std(S.W(onIdx));
    % 
    %     fprintf('\n=== Wind-on velocity statistics (t >= 0) ===\n');
    %     fprintf('Mean:  U=%.3f  V=%.3f  W=%.3f m/s\n', Um, Vm, Wm);
    %     fprintf('RMS:   U=%.3f  V=%.3f  W=%.3f m/s\n', Urms, Vrms, Wrms);
    %     fprintf('Turbulence intensity: Tu_u=%.1f%%  Tu_v=%.1f%%  Tu_w=%.1f%%\n', ...
    %         100*Urms/Um, 100*Vrms/Um, 100*Wrms/Um);
    % end
    % 
    % %% Power spectral density of U (wind-on portion), if pwelch is available
    % if hasHotwire && exist('pwelch', 'file')
    %     onIdx = t >= 0;
    %     Uon = S.U(onIdx);
    %     ton = t(onIdx);
    %     fs = 1 / mean(diff(ton), 'omitnan');
    % 
    %     figure('Name', 'PSD of U', 'Position', [900 720 700 350]);
    %     [psd, freq] = pwelch(Uon - mean(Uon), [], [], [], fs);
    %     loglog(freq, psd, 'LineWidth', 1.5);
    %     xlabel('Frequency (Hz)'); ylabel('PSD (m^2/s^2/Hz)');
    %     title('Power spectral density of U -- wind-on portion');
    %     grid on;
    % end
end
