function C = plot_komega_gauge_csig(gauge, csig, win_dur, step)
% PLOT_KOMEGA_GAUGE_CSIG  Peak wavenumber from the CSIG slope field vs peak
% frequency from the wave gauge, on the (k, omega) plane with the linear
% dispersion relation overlaid.
%
% k and omega come from DIFFERENT instruments here (k: spatial FFT of Sx
% over the image, as plot_kx_t_slopes_csig.m; omega: temporal spectrum of
% the gauge, as plot_f_t_eta.m), so their agreement with the dispersion
% curve is a real test -- unlike plot_k_t_eta.m, whose k is computed FROM
% the gauge's f and lies on the curve by construction. No Doppler shift from
% the wind drift is included yet; a systematic offset above the curve that
% grows with wind is the signature of that drift.
%
% Figures:
%   1. peak pairs (CSIG k_p, gauge omega_p) on the (k, omega) plane, coloured by t
%   2. full |Sx(kx, omega)|^2 of the CSIG record (as plot_kx_omega_Sx_csig.m,
%      omega axis) with dispersion curves and the peak pairs overlaid, and the
%      gauge S_eta(omega) beside it on the same omega axis
%   3. omega_p vs t: gauge measured vs CSIG k_p pushed through dispersion
%
%   gauge = load('wave_ramp_20260928_120000.mat');     % run_wavegauge_fan_calibration.m
%   csig  = struct('Sx', Sx, 'time', time_subset, 'dx', dx);   % from make_all_plots_csig.m
%   plot_komega_gauge_csig(gauge, csig)           % 4 s windows, one point per second
%   plot_komega_gauge_csig(gauge, csig, 4, 0.5)
%
% csig.time must be camera time as in make_all_plots_csig.m, i.e.
% (frame_number - 1) * dt, so frame n <-> gauge.frameT(n). If the gauge file
% has no frameT, the software camera start (camOnT) is used instead.
%
% Resolution limits (printed on run): omega is quantised to 2*pi/win_dur by
% the gauge window, k to 2*pi/L by the image width L.

if nargin < 3 || isempty(win_dur), win_dur = 4; end
if nargin < 4 || isempty(step),    step = 1;    end

F_LO     = 0.5;   % Hz, gauge peak search floor (as plot_f_t_eta.m)
K_LO_BIN = 3;     % skip DC and the first bin (domain-scale leakage) in the CSIG peak search
g = 9.81; gamma = 7.4e-5;   % as Slope_Gauge/util/dispersion_curves.m

%% CSIG time -> wind clock (the gauge's time base)
camFreq = 50;
if isfield(gauge, 'camConfig') && ~isempty(gauge.camConfig), camFreq = gauge.camConfig(1).freq; end
tCam = csig.time(:).';
if isfield(gauge, 'frameT') && ~isempty(gauge.frameT)
    n = round(tCam * camFreq) + 1;
    if any(n < 1 | n > numel(gauge.frameT))
        error('CSIG frames %d..%d fall outside the %d frames logged in the gauge file.', ...
            min(n), max(n), numel(gauge.frameT));
    end
    tW = gauge.frameT(n).';
    fprintf('CSIG -> wind clock via gauge frameT (%s).\n', gauge.frameSource);
else
    tW = gauge.camOnT(1) + tCam;
    fprintf('CSIG -> wind clock via software camera start (+-~50-70 ms).\n');
end

%% CSIG: y-averaged kx power spectrum of every frame (as plot_kx_t_slopes_csig.m)
[~, Nx, Nt] = size(csig.Sx);
Nk    = floor(Nx/2) + 1;
kx    = 2*pi * (0:Nk-1) / (Nx * csig.dx);
win_x = hann(Nx).';
Pk = zeros(Nk, Nt);
for j = 1:Nt
    fx = double(csig.Sx(:,:,j));
    fx = fx - mean(fx, 2);
    F  = fft(fx .* win_x, [], 2);
    Pk(:, j) = mean(abs(F(:, 1:Nk)).^2, 1).';
end

%% Window centres common to both instruments
tc = (tW(1) + win_dur/2) : step : (tW(end) - win_dur/2);
if isempty(tc), error('CSIG record shorter than one %.1f s window.', win_dur); end
nc = numel(tc);

kp = nan(1, nc);
for j = 1:nc
    in = abs(tW - tc(j)) <= win_dur/2;
    P  = mean(Pk(:, in), 2);
    [~, i] = max(P(K_LO_BIN:end));
    kp(j) = kx(i + K_LO_BIN - 1);
end

%% Gauge: peak frequency on the same windows
S = eta_spectrogram(gauge, win_dur, [tW(1) tW(end)]);
iLo = find(S.f >= F_LO, 1);
[~, iPk] = max(S.P(iLo:end, :), [], 1);
fpG = S.f(iPk + iLo - 1).';
omegaG = 2*pi * interp1(S.t, fpG, tc, 'nearest', 'extrap');   % S.t steps 0.25 s

omegaFromK   = sqrt(g*kp + gamma*kp.^3);   % what the gauge SHOULD see for the CSIG k
omegaFromK_g = sqrt(g*kp);

fprintf('Resolution: d(omega) = %.2f rad/s (gauge, %g s window), dk = %.1f rad/m (CSIG, L = %.1f cm)\n', ...
    2*pi/win_dur, win_dur, kx(2), 2*pi/kx(2)*100);

C = struct('t', tc, 'kp', kp, 'omega_gauge', omegaG, ...
           'omega_disp', omegaFromK, 'omega_disp_gravity', omegaFromK_g, 'win_dur', win_dur);

%% (k, omega) plane
kmax = max([kp, sqrt(max(omegaG)^2/g)]) * 1.4;
kc   = linspace(0, kmax, 500);

figure('Name', 'k-omega: CSIG k_p vs gauge omega_p', 'Position', [100 100 800 650], 'Color', 'w');
hold on;
plot(kc, sqrt(g*kc + gamma*kc.^3), 'k-',  'LineWidth', 1.5, 'DisplayName', '$\omega^2 = gk + \gamma k^3$');
plot(kc, sqrt(g*kc),               'k--', 'LineWidth', 1.2, 'DisplayName', '$\omega^2 = gk$');
scatter(kp, omegaG, 40, tc, 'filled', 'DisplayName', 'CSIG $k_p$, gauge $\omega_p$');
cb = colorbar; cb.Label.String = '$t$ (s)'; cb.Label.Interpreter = 'latex'; cb.Label.FontSize = 14;
colormap(gca, pick_cmap());
xlabel('$k_x$ (rad m$^{-1}$)', 'Interpreter', 'latex');
ylabel('$\omega$ (rad s$^{-1}$)', 'Interpreter', 'latex');
title(sprintf('Peak $(k, \\omega)$, %g s windows', win_dur), 'Interpreter', 'latex');
legend('Interpreter', 'latex', 'Location', 'northwest');
xlim([0 kmax]); grid on;
set(gca, 'fontsize', 14, 'fontname', 'times');

%% Full (kx, omega) spectrum of Sx, with the gauge omega-spectrum alongside
% Background: |Sx(kx, omega)|^2 over the whole CSIG record, same '2d' method
% as Slope_Gauge/plotting/plot_kx_omega_Sx_csig.m (y-averaged, Hann in x and
% t), but on an omega axis. Right panel: the gauge's S_eta(omega) over the
% same time span, sharing the omega axis, so the gauge peak can be read
% straight across onto the (k, omega) map. The peak pairs from the first
% figure are overlaid as dots.
Fs = 1 / median(diff(csig.time));
Sx_xt = reshape(mean(double(csig.Sx), 1), Nx, Nt);          % [Nx x Nt]
Sx_xt = (Sx_xt - mean(Sx_xt(:))) .* (hann(Nx) * hann(Nt).');
P2    = fftshift(abs(fft2(Sx_xt)).^2);
kx_full = 2*pi * ((0:Nx-1) - floor(Nx/2)) / (Nx * csig.dx);  % matches fftshift for odd N too
f_full  = ((0:Nt-1) - floor(Nt/2)) * Fs / Nt;
ip      = f_full > 0;
omegaI  = 2*pi * f_full(ip);
PkwI    = P2(:, ip);

Sg     = mean(S.P, 2) / (2*pi);   % S_eta(omega) = S_eta(f)/(2 pi), time-averaged over the record
omegaS = 2*pi * S.f;

omegaMax = min(max(omegaI), max(omegaS));
kLim     = min(max(abs(kx_full)), 1.2 * fzero(@(k) sqrt(g*k + gamma*k^3) - omegaMax, [1e-3 1e5]));
kd       = linspace(0, kLim, 600);

figure('Name', 'k-omega spectrum: CSIG Sx with gauge', 'Position', [80 80 1200 650], 'Color', 'w');
tl = tiledlayout(1, 4, 'TileSpacing', 'compact', 'Padding', 'compact');

ax1 = nexttile(tl, [1 3]);
imagesc(kx_full, omegaI, log10(PkwI.' + eps));
set(gca, 'YDir', 'normal'); colormap(ax1, pick_cmap());
cb = colorbar('westoutside');
cb.Label.String = '$\log_{10}|S_x(k_x,\omega)|^2$'; cb.Label.Interpreter = 'latex'; cb.Label.FontSize = 14;
hold on;
h1 = plot( kd, sqrt(g*kd + gamma*kd.^3), 'w-',  'LineWidth', 1.5);
     plot(-kd, sqrt(g*kd + gamma*kd.^3), 'w-',  'LineWidth', 1.5);
h2 = plot( kd, sqrt(g*kd),               'w--', 'LineWidth', 1.2);
     plot(-kd, sqrt(g*kd),               'w--', 'LineWidth', 1.2);
h3 = scatter(kp, omegaG, 30, 'c', 'filled', 'MarkerEdgeColor', 'k');
legend([h1 h2 h3], {'$\omega^2 = gk + \gamma k^3$', '$\omega^2 = gk$', 'CSIG $k_p$, gauge $\omega_p$'}, ...
       'Interpreter', 'latex', 'Location', 'northwest', 'TextColor', 'w', 'Color', 'none');
xlim([-kLim kLim]); ylim([0 omegaMax]);
xlabel('$k_x$ (rad m$^{-1}$)', 'Interpreter', 'latex');
ylabel('$\omega$ (rad s$^{-1}$)', 'Interpreter', 'latex');
title(sprintf('$(k_x,\\omega)$ spectrum of $S_x$, $t$ = %.0f--%.0f s', tW(1), tW(end)), 'Interpreter', 'latex');
set(gca, 'fontsize', 14, 'fontname', 'times');

ax2 = nexttile(tl);
semilogx(Sg, omegaS, 'b-', 'LineWidth', 1.3); hold on;
yline(median(omegaG), 'c-', 'LineWidth', 1.5);
ylim([0 omegaMax]); grid on;
xlabel(sprintf('$S_\\eta(\\omega)$ (%s$^2$ s rad$^{-1}$)', S.unit), 'Interpreter', 'latex');
title('Wave gauge', 'Interpreter', 'latex');
set(gca, 'fontsize', 14, 'fontname', 'times', 'YTickLabel', []);
linkaxes([ax1 ax2], 'y');

%% Same comparison through time
figure('Name', 'omega_p vs time: gauge vs CSIG via dispersion', 'Position', [100 100 1000 450], 'Color', 'w');
hold on;
plot(tc, omegaG,       'b-o', 'LineWidth', 1.3, 'MarkerSize', 4, 'DisplayName', 'gauge $\omega_p$ (measured)');
plot(tc, omegaFromK,   'r-s', 'LineWidth', 1.3, 'MarkerSize', 4, 'DisplayName', 'CSIG $k_p \to \omega$ ($gk + \gamma k^3$)');
plot(tc, omegaFromK_g, 'r--', 'LineWidth', 1.0, 'DisplayName', 'CSIG $k_p \to \omega$ ($gk$)');
xlabel('$t$ (s since wind start)', 'Interpreter', 'latex');
ylabel('$\omega_p$ (rad s$^{-1}$)', 'Interpreter', 'latex');
legend('Interpreter', 'latex', 'Location', 'best'); grid on;
set(gca, 'fontsize', 14, 'fontname', 'times');
end


function cm = pick_cmap()
    if exist('inferno', 'file'), cm = inferno; else, cm = parula; end
end
