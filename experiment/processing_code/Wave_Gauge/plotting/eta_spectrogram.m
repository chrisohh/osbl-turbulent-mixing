function S = eta_spectrogram(data, win_dur, t_range)
% ETA_SPECTROGRAM  Sliding-window PSD of the wave-gauge elevation eta(t).
% Shared by plot_f_t_eta and plot_k_t_eta -- the gauge counterpart of the
% windowed spectra in Slope_Gauge/plotting/plot_f_t_slopes_csig.m.
%
%   S = eta_spectrogram(data)                 % 4 s windows, whole record
%   S = eta_spectrogram(data, win_dur)
%   S = eta_spectrogram(data, win_dur, t_range)
%
% data    : struct loaded from a wave_ramp_*.mat (run_wavegauge_fan_calibration.m),
%           needs hwT_wind, eta, etaUnit (camOnT/camOffT used if present)
% win_dur : window length in seconds (default 4, as in Leckler et al.)
% t_range : [t1 t2] on the wind clock, or 'cam' = the camera window only
%           (default [] = whole record)
%
% Returns S with fields
%   t     window-centre times (s, wind clock)       [1 x Nt]
%   f     frequency (Hz), 0..F_MAX                  [Nf x 1]
%   P     one-sided PSD of eta, in m^2/Hz (or V^2/Hz if the gauge was not
%         calibrated)                               [Nf x Nt]
%   unit  'm' or 'V';  camOn/camOff  camera window (NaN if none)

if nargin < 2 || isempty(win_dur), win_dur = 4; end
if nargin < 3, t_range = []; end

STEP  = 0.25;   % s between window centres
F_MAX = 10;     % Hz, the wire gauge does not resolve much above this

t   = data.hwT_wind(:);
eta = double(data.eta(:));
switch data.etaUnit
    case 'cm', eta = eta / 100; S.unit = 'm';
    case 'm',  S.unit = 'm';
    otherwise, S.unit = data.etaUnit;
end

S.camOn = NaN; S.camOff = NaN;
if isfield(data, 'camOnT') && ~isempty(data.camOnT)
    S.camOn = data.camOnT(1); S.camOff = data.camOffT(1);
end

if ischar(t_range) && strcmpi(t_range, 'cam')
    if isnan(S.camOn), error('No camera window in this file.'); end
    t_range = [S.camOn S.camOff];
end
if ~isempty(t_range)
    m = t >= t_range(1) & t <= t_range(2);
    t = t(m); eta = eta(m);
end

fs   = 1 / median(diff(t));
nWin = round(win_dur * fs);
if numel(eta) < nWin
    error('Record (%.1f s) shorter than one %.1f s window.', numel(eta)/fs, win_dur);
end

% Remove the slow set-up/drift (the ramp raises the mean level) so it does
% not leak into the lowest bins -- same role as the row-mean removal in the
% CSIG spatial FFT.
eta = eta - movmean(eta, nWin);

nStep = max(round(STEP * fs), 1);
nfft  = max(2^nextpow2(nWin), nWin);
[~, f, tc, P] = spectrogram(eta, hann(nWin), nWin - nStep, nfft, fs);

keep  = f <= F_MAX;
S.f   = f(keep);
S.P   = P(keep, :);
S.t   = tc(:).' + t(1);   % spectrogram times are from the first sample
S.win_dur = win_dur;
S.fs  = fs;
end
