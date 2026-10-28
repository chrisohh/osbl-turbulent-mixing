function hw = lif_hw_series(H, avg_s)
%LIF_HW_SERIES Hot-wire velocity series for the lif_surf_figure panel.
%
%   hw = LIF_HW_SERIES(H, avg_s)
%
% H     - struct loaded from a run_experiment.m hot-wire file
% avg_s - block-average length (s), e.g. 0.02 = one camera frame at 50 Hz.
%         The raw record is 4 kHz; plotting it whole is slow and unreadable.
%
% hw.t     - block centres, s since wind start (hwT_wind)
% hw.Y     - [nt x nseries] block means: U (hot-wire) and U_ref (reference
%            probe), whichever were logged
% hw.names - LaTeX legend names
% hw.kind  - 'velocity' (U / U_ref, m/s) or 'eta' (wave gauge, cm: runs with
%            a wave gauge and no hot-wire, e.g. run 102)
% hw.ylabel
%
% The sample rate is taken from the timestamps, so 1 kHz records (run 96's
% hot-wire, the run-89 reference probe, wave gauges) need nothing special.

    t = H.hwT_wind(:);
    Y = [];  names = {};
    if isfield(H, 'U') && ~isempty(H.U),         Y = [Y H.U(:)];     names{end+1} = '$U$ (hot-wire)'; end
    if isfield(H, 'U_ref') && ~isempty(H.U_ref), Y = [Y H.U_ref(:)]; names{end+1} = '$U_{\rm ref}$'; end
    hw.kind = 'velocity';  hw.ylabel = 'velocity (m/s)';
    if isempty(Y) && isfield(H, 'E_ref') && ~isempty(H.E_ref)
        % Older runs (e.g. 89, 9/10) saved the reference probe only as the
        % voltage E_ref.  Convert as run_experiment.m does now: T29 block of
        % the run's CAL_FILE, origin-to-floor line below the calibrated
        % range, clamped above it.
        Uref = [];
        try
            addpath(fullfile(fileparts(mfilename('fullpath')), '..', '..', 'CTA'));   % parse_probe_section
            cf = H.CAL_FILE;
            if ~isfile(cf)          % saved from an older checkout (D:\Chris\osbl-turbulent-mixing\...)
                k = strfind(cf, 'osbl-turbulent-mixing');
                repo = fullfile(fileparts(mfilename('fullpath')), '..', '..', '..', '..');
                if ~isempty(k), cf = fullfile(repo, cf(k(1) + numel('osbl-turbulent-mixing') + 1 : end)); end
            end
            E = H.E_ref(:);
            try
                c = parse_probe_section(cf, 'T29');
                Uref = polyval(fliplr(c.C(1,:)), E);
                lo = E < c.E_floor(1);  Uref(lo) = (c.U_floor(1) / c.E_floor(1)) * E(lo);
                Uref(E > c.E_ceil(1)) = c.U_ceil(1);
            catch
                % CAL_FILE without a T29 block (older runs, e.g. 89: probe3.txt
                % holds only the 55R95 hot-wire): the probe's certificate table
                Uref = convert_Eref2Uref(E);
            end
        catch ME
            warning('lif_hw_series:Eref', 'E_ref could not be converted (%s); plotting the voltage.', ME.message);
        end
        if ~isempty(Uref)
            Y = Uref;  names = {'$U_{\rm ref}$ (from $E_{\rm ref}$)'};
        else
            Y = H.E_ref(:);  names = {'$E_{\rm ref}$ (V)'};
            hw.kind = 'voltage';  hw.ylabel = 'E_{ref} (V)';
        end
    end
    if isempty(Y) && isfield(H, 'eta') && ~isempty(H.eta)
        Y = H.eta(:);  names = {'$\eta$ (wave gauge)'};
        hw.kind = 'eta';  hw.ylabel = '\eta (cm)';
    end
    if isempty(Y), error('lif_hw_series:empty', 'No U, U_ref or eta in the file.'); end

    dt = median(diff(t));
    n  = max(1, round(avg_s / dt));
    nb = floor(numel(t) / n);
    idx = 1 : nb*n;
    hw.t = mean(reshape(t(idx), n, nb), 1).';
    hw.Y = squeeze(mean(reshape(Y(idx, :), n, nb, []), 1));
    if isrow(hw.Y), hw.Y = hw.Y.'; end
    hw.names = names;
end
