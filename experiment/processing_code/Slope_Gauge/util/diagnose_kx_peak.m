function diagnose_kx_peak(data, cache_dir, rest_frames, kx_suspect, dt)
% DIAGNOSE_KX_PEAK  Diagnose a bright vertical stripe at fixed kx in the
% (kx, f) spectrogram of Sx.
%
% Usage:
%   diagnose_kx_peak(data, cache_dir, rest_frames, kx_suspect, dt)
%
% Inputs:
%   data         : struct with Sx, dx (from make_all_plots_csig)
%   cache_dir    : path to per-frame cache (for Sx_offset check)
%   rest_frames  : frame indices known to be still water (e.g. 1:50)
%   kx_suspect   : wavenumber of bright stripe (rad/m, e.g. 200)
%   dt           : frame interval in seconds (1/Fs)
%
% Five panels:
%   1. kx power spectrum: rest mean (Sx_offset) vs analysis cube row-avg
%   2. kx power spectrum of single rest frame vs. (rest frame - Sx_offset)
%   3. Spatial map of one Sx row with kx_suspect component overlaid
%   4. Phase map of Sx at kx_suspect -> constant phase = fixed pattern
%   5. Time series of |FFT at kx_suspect| to see if artifact is static

dx  = data.dx;
Sx  = data.Sx;
[Ny, Nx, Nt] = size(Sx);

kx_axis = 2*pi * (0:floor(Nx/2)) / (Nx * dx);   % one-sided, rad/m
[~, ikx] = min(abs(kx_axis - kx_suspect));
fprintf('kx_suspect = %.0f rad/m  -> lambda = %.2f cm\n', ...
        kx_suspect, 2*pi/kx_suspect * 100);
fprintf('Nearest bin: kx = %.1f rad/m  (index %d)\n', kx_axis(ikx), ikx);

%% 1. Load Sx_offset (rest mean) and compare kx spectra
fprintf('Loading rest cube (%d frames) ...\n', numel(rest_frames));
slope_offset_file = strcat([cache_dir '\slope_offset.mat']);
    load(slope_offset_file)

% Row-averaged kx power spectra
win_x = hann(Nx).';
P_offset   = mean(abs(fft(bsxfun(@times, Sx_offset - mean(Sx_offset,2), win_x), [], 2)).^2, 1);
P_analysis = mean(mean(abs(fft(bsxfun(@times, ...
    double(Sx(:,:,1:min(50,Nt))) - mean(double(Sx(:,:,1:min(50,Nt))),2), win_x), [], 2)).^2, 1), 3);

% Single rest frame before and after offset subtraction
Sx_r1 = double(Sx_rest(:,:,1));
Sx_r1_corr = Sx_r1 - Sx_offset;
P_r1      = mean(abs(fft(bsxfun(@times, Sx_r1      - mean(Sx_r1,2),      win_x), [],2)).^2,1);
P_r1_corr = mean(abs(fft(bsxfun(@times, Sx_r1_corr - mean(Sx_r1_corr,2), win_x), [],2)).^2,1);

%% 2. Time series of |FFT| at kx_suspect across all frames
P_t = zeros(1, Nt);
for t = 1:Nt
    row_fft  = abs(fft(double(Sx(:,:,t)) .* win_x, [], 2));
    P_t(t)   = mean(row_fft(:, ikx).^2);
end
time_vec = (0:Nt-1) * dt;

%% 3. Spatial map of kx_suspect component (phase)
% Use first analysis frame
frame_ex = double(Sx(:,:,1));
F_ex     = fft(frame_ex .* win_x, [], 2);
phase_map = angle(F_ex(:, ikx));      % [Ny x 1] phase at each row

%% 4. One row of Sx with filtered kx_suspect overlaid
row_ex   = double(Sx(round(Ny/2), :, 1));
F_row    = fft(row_ex .* hann(Nx).');
F_bp     = zeros(size(F_row));
bw       = max(1, round(Nx * dx / (2*pi/kx_suspect) * 0.1));  % ±10% bandwidth
F_bp(ikx-bw:ikx+bw) = F_row(ikx-bw:ikx+bw);
F_bp(end-ikx-bw:end-ikx+bw) = F_row(end-ikx-bw:end-ikx+bw);
row_bp   = real(ifft(F_bp));
x_cm     = (0:Nx-1) * dx * 100;

%% Plotting
fig = figure('Name','kx stripe diagnostic','Position',[50 50 1400 900],'Color','w');

% --- Panel 1: kx spectra comparison ---
subplot(2,3,1);
semilogy(kx_axis, P_offset(1:numel(kx_axis)), 'b', 'LineWidth',1.5); hold on;
semilogy(kx_axis, P_analysis(1:numel(kx_axis)), 'r', 'LineWidth',1.5);
xline(kx_suspect, 'k--', sprintf('k=%d',round(kx_suspect)), 'LabelVerticalAlignment','bottom');
xlabel('k_x (rad/m)'); ylabel('Power'); grid on;
title('kx spectrum: offset (blue) vs analysis (red)');
legend('Sx\_offset','Analysis frames','Location','best');
set(gca,'fontsize',11,'fontname','times');

% --- Panel 2: rest frame before/after correction ---
subplot(2,3,2);
semilogy(kx_axis, P_r1(1:numel(kx_axis)), 'Color',[.6 .6 .6],'LineWidth',1.5); hold on;
semilogy(kx_axis, P_r1_corr(1:numel(kx_axis)), 'b','LineWidth',1.5);
xline(kx_suspect, 'k--');
xlabel('k_x (rad/m)'); ylabel('Power'); grid on;
title('Rest frame: raw (grey) vs offset-subtracted (blue)');
legend('Raw','Corrected','Location','best');
set(gca,'fontsize',11,'fontname','times');

% --- Panel 3: spatial row + bandpassed kx_suspect component ---
subplot(2,3,3);
plot(x_cm, row_ex, 'Color',[.7 .7 .7],'LineWidth',0.8); hold on;
plot(x_cm, row_bp, 'r','LineWidth',1.8);
xlabel('x (cm)'); ylabel('S_x');
title(sprintf('Mid-row S_x: grey=full, red=%.0f rad/m component (\\lambda=%.1f cm)', ...
      kx_suspect, 2*pi/kx_suspect*100));
grid on;
set(gca,'fontsize',11,'fontname','times');

% --- Panel 4: phase at kx_suspect across y rows ---
subplot(2,3,4);
plot(phase_map, 1:Ny, 'b.','MarkerSize',4);
xlabel('Phase (rad)'); ylabel('y row');
title(sprintf('Phase of k_x=%.0f component across y\n(constant = fixed pattern, random = noise)', kx_suspect));
xlim([-pi pi]); grid on;
set(gca,'fontsize',11,'fontname','times');

% --- Panel 5: time evolution of power at kx_suspect ---
subplot(2,3,5);
plot(time_vec, P_t, 'b','LineWidth',1);
xlabel('t (s)'); ylabel(sprintf('|FFT|^2 at k_x=%.0f',kx_suspect));
title('Power at k\_x suspect vs time');
yline(mean(P_t(1:min(numel(rest_frames),Nt))), 'r--', 'Rest mean');
grid on;
set(gca,'fontsize',11,'fontname','times');

% --- Panel 6: 2D spatial map of Sx with kx_suspect highlighted ---
subplot(2,3,6);
x_cm_ax = (0:Nx-1)*dx*100;
y_cm_ax = (0:Ny-1)*dx*100;   % approximate
imagesc(x_cm_ax, y_cm_ax, Sx_offset);
axis equal tight; set(gca,'YDir','normal');
colormap(gca, brewermap([],'Spectral')); colorbar;
title('Sx\_offset spatial map (fixed pattern)');
xlabel('x (cm)'); ylabel('y (cm)');
set(gca,'fontsize',11,'fontname','times');

sgtitle(sprintf('Diagnostic: kx = %.0f rad/m  (\\lambda = %.2f cm)', ...
        kx_suspect, 2*pi/kx_suspect*100), 'FontSize',13,'FontName','times');

%% Summary printout
fprintf('\n--- DIAGNOSTIC SUMMARY ---\n');
fprintf('Power at kx=%.0f in Sx_offset      : %.2e\n', kx_suspect, P_offset(ikx));
fprintf('Power at kx=%.0f in analysis frames : %.2e\n', kx_suspect, P_analysis(ikx));
fprintf('Power at kx=%.0f in raw rest frame  : %.2e\n', kx_suspect, P_r1(ikx));
fprintf('Power at kx=%.0f after correction   : %.2e\n', kx_suspect, P_r1_corr(ikx));
ratio = P_r1_corr(ikx) / P_r1(ikx);
fprintf('Residual after correction: %.1f%%\n', ratio*100);
if ratio < 0.1
    fprintf('  -> Artifact is STATIC. Rest subtraction removes it.\n');
elseif ratio < 0.5
    fprintf('  -> Artifact is PARTIALLY static. Some time-variation remains.\n');
else
    fprintf('  -> Artifact is TIME-VARYING. Rest subtraction does not help.\n');
end
phase_std = std(phase_map);
fprintf('Phase std across y rows: %.2f rad  ', phase_std);
if phase_std < 0.3
    fprintf('(coherent across y -> fixed spatial pattern)\n');
else
    fprintf('(incoherent across y -> noise or real waves)\n');
end
end
