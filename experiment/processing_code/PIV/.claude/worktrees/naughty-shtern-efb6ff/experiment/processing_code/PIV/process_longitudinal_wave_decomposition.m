clear; clc;

%% =========================================================================
%% PARAMETERS
%% =========================================================================
root    = get_server_root();
addpath([root 'codes\202204 - IR PIV to compute surface velocity\']);
addpath([root 'codes\202204 - Transverse structures PIV\']);

Fs_PIV  = 7.2;          % Hz
delay   = 73 + 5;       % s  (imaging starts delay s after trigger)
DX      = 6.7861e-05;   % m/pixel
GS      = 4;            % PIV window size (pixels); physical grid = DX*GS
t_wind_onset = 40;      % s  trigger→wind onset (for JFM time axis)
t_snap_target = 55;     % s  (JFM) snapshot to visualise

%% =========================================================================
%% SECTION 1 — LOAD / REPROCESS LONGITUDINAL PIV CUBE
%% =========================================================================
cube_file = [root 'codes\202204 - Transverse structures PIV\R2_EXP1_Longitudinal_PIV_Cube.mat'];

if exist(cube_file, 'file')
    fprintf('Loading saved longitudinal cube: %s\n', cube_file);
    load(cube_file);                     % loads STAT_LG_R2_EXP1
else
    fprintf('Cube file not found — reprocessing from PIVMat...\n');
    dirin = [root 'Longitudinal\PIV\ExpLCL_2_01\PIVMat\'];
    STAT_LG_R2_EXP1 = process_one_directory_transverse_PIV(dirin);

    % Time, X (along-wind), Z (depth) axes
    STAT_LG_R2_EXP1.time = delay + STAT_LG_R2_EXP1.filenumber / Fs_PIV;
    STAT_LG_R2_EXP1.X    = (1:length(STAT_LG_R2_EXP1.u(1,1,:))) * DX * GS;  % along-wind (m)
    STAT_LG_R2_EXP1.Y    = (1:length(STAT_LG_R2_EXP1.u(1,:,1))) * DX * GS;  % depth      (m)

    save(cube_file, 'STAT_LG_R2_EXP1', '-v7.3');
    fprintf('Saved to %s\n', cube_file);
end

STAT = STAT_LG_R2_EXP1;

% --- Convenience variables ---
% process_one_directory_transverse_PIV convention (confirmed via imagesc call):
%   squeeze(STAT.u(t,:,:)) is [N_2nd x N_3rd]
%   imagesc(STAT.X, STAT.Y, squeeze(STAT.u(t,:,:))) puts STAT.X on x-axis (cols)
%   and STAT.Y on y-axis (rows), so:
%     STAT.Y → 2nd dim → horizontal (along-wind x)
%     STAT.X → 3rd dim → vertical   (depth z)
%
% For pcolor(x_piv, z_piv, data) we need data = [Nz x Nx], so squeeze(...)' is required.

t_piv = STAT.time - t_wind_onset;  % JFM time [Nt x 1]
x_piv = STAT.Y(:);                 % along-wind [Nx x 1] m   (2nd dim)
z_piv = STAT.X(:);                 % depth      [Nz x 1] m   (3rd dim), positive downward
%   z_piv(1) = DX*GS ≈ first grid cell below image top.
%   Image top ≈ mean water surface by camera design.
%   → z_piv is already depth below mean water level (check: fprintf below).

% Cube dimensions: u(time, x_along_wind, z_depth)
t_piv = STAT.time - t_wind_onset;  % JFM time  [Nt x 1]
x_piv = STAT.X(:);                 % along-wind [Nx x 1]  m
z_piv = STAT.Y(:);                % depth, negative downward [Nz x 1]  m

[Nt, Nx, Nz] = size(STAT.u);
fprintf('Longitudinal cube: Nt=%d  Nx=%d  Nz=%d\n', Nt, Nx, Nz);
fprintf('  x: %.1f–%.1f mm,  z: %.1f–%.1f mm\n', ...
    1e3*x_piv(1), 1e3*x_piv(end), 1e3*z_piv(1), 1e3*z_piv(end));

%% =========================================================================
%% SECTION 2 — LOAD ETA
%% =========================================================================
dirin_eta = [root 'Longitudinal\PIV\ExpLCL_2_01\PIVRaw\EXTRACTED_SURFACES\'];
load([dirin_eta 'A.mat']);                           % A1 [Nx_eta x Nt_eta]
load([dirin_eta 'ExpLCL_2_01_Surface_000.mat']);     % surface_x1

x_eta = surface_x1(:);
t_eta = delay + (0:size(A1,2)-1)/Fs_PIV - t_wind_onset;
eta   = A1;   % [Nx_eta x Nt_eta]  m
clear A1 surface_x1;
fprintf('eta: Nx=%d  Nt=%d  t=%.1f–%.1f s\n', length(x_eta), length(t_eta), t_eta(1), t_eta(end));

%% =========================================================================
%% SECTION 2b — VERIFICATION
%%   Run this section after loading STAT and eta.
%%   Each check prints PASS / WARN / FAIL with a diagnosis.
%% =========================================================================
fprintf('\n========== VERIFICATION ==========\n');
ok = true;

% ---- (1) STAT dimension ordering ----------------------------------------
% Expected: size(STAT.u) = [Nt, Nx_along_wind, Nz_depth]
% From user's imagesc call: STAT.Y is 2nd dim (along-wind), STAT.X is 3rd dim (depth)
assert_dim = (length(STAT.Y) == Nx) && (length(STAT.X) == Nz);
if assert_dim
    fprintf('  [PASS] STAT dimensions: Nt=%d  Nx(along-wind)=%d  Nz(depth)=%d\n', Nt, Nx, Nz);
else
    fprintf('  [FAIL] STAT.Y length (%d) != Nx (%d) or STAT.X length (%d) != Nz (%d)\n', ...
        length(STAT.Y), Nx, length(STAT.X), Nz);
    fprintf('         Check whether X/Y assignment is swapped.\n');
    ok = false;
end

% ---- (2) Grid spacing matches DX*GS -------------------------------------
dx_x = mean(diff(x_piv));
dx_z = mean(diff(z_piv));
expected_gs = DX * GS;
if abs(dx_x - expected_gs)/expected_gs < 0.01 && abs(dx_z - expected_gs)/expected_gs < 0.01
    fprintf('  [PASS] Grid spacing: dx=%.4f mm  dz=%.4f mm  (expected %.4f mm)\n', ...
        dx_x*1e3, dx_z*1e3, expected_gs*1e3);
else
    fprintf('  [WARN] Grid spacing mismatch: dx=%.4f mm, dz=%.4f mm, expected=%.4f mm\n', ...
        dx_x*1e3, dx_z*1e3, expected_gs*1e3);
    fprintf('         DX or GS parameter may be wrong.\n');
end

% ---- (3) eta units: metres vs pixels ------------------------------------
eta_rms_mm = rms(eta(:)) * 1e3;
eta_mean_mm = mean(eta(:)) * 1e3;
if eta_rms_mm < 20
    fprintf('  [PASS] eta RMS = %.3f mm  (capillary-gravity range ~1–5 mm) → units = metres\n', eta_rms_mm);
else
    fprintf('  [FAIL] eta RMS = %.1f mm — likely still in PIXELS (expect ~1000–3000)\n', eta_rms_mm);
    fprintf('         Fix:  eta = (MWL - eta_pixels) * DX\n');
    fprintf('         MWL for Ramp2: 2468 px  (from Fab1_Main_PIVCalc MeanWL vector)\n');
    ok = false;
end

% ---- (4) eta mean ≈ 0 (MWL already subtracted) -------------------------
if abs(eta_mean_mm) < 1.0
    fprintf('  [PASS] eta mean = %.4f mm ≈ 0  (MWL already subtracted)\n', eta_mean_mm);
else
    fprintf('  [WARN] eta mean = %.3f mm  (non-zero; subtracting before use)\n', eta_mean_mm);
    fprintf('         This is fine — code uses eta_fluct = eta_snap - eta_mean.\n');
end

% ---- (5) eta spatial coverage of PIV FOV --------------------------------
x_piv_min = x_piv(1);   x_piv_max = x_piv(end);
x_eta_min = x_eta(1);   x_eta_max = x_eta(end);
if x_eta_min <= x_piv_min && x_eta_max >= x_piv_max
    fprintf('  [PASS] eta x-range [%.1f, %.1f] mm covers PIV FOV [%.1f, %.1f] mm\n', ...
        x_eta_min*1e3, x_eta_max*1e3, x_piv_min*1e3, x_piv_max*1e3);
else
    fprintf('  [WARN] eta x-range [%.1f, %.1f] mm does NOT fully cover PIV FOV [%.1f, %.1f] mm\n', ...
        x_eta_min*1e3, x_eta_max*1e3, x_piv_min*1e3, x_piv_max*1e3);
    fprintf('         Extrapolation will be used at edges.\n');
end

% ---- (6) Time overlap between eta and PIV -------------------------------
t_overlap_start = max(t_eta(1),  t_piv(1));
t_overlap_end   = min(t_eta(end), t_piv(end));
overlap_frac = (t_overlap_end - t_overlap_start) / (t_piv(end) - t_piv(1));
if overlap_frac > 0.95
    fprintf('  [PASS] Time overlap: %.1f–%.1f s  (%.0f%% of PIV record)\n', ...
        t_overlap_start, t_overlap_end, 100*overlap_frac);
elseif overlap_frac > 0
    fprintf('  [WARN] Partial time overlap: %.1f–%.1f s  (%.0f%% of PIV record)\n', ...
        t_overlap_start, t_overlap_end, 100*overlap_frac);
    fprintf('         Check delay and t_wind_onset parameters.\n');
else
    fprintf('  [FAIL] No time overlap between eta (%.1f–%.1f s) and PIV (%.1f–%.1f s)\n', ...
        t_eta(1), t_eta(end), t_piv(1), t_piv(end));
    fprintf('         Check delay (%g s) and t_wind_onset (%g s).\n', delay, t_wind_onset);
    ok = false;
end

% ---- (7) Typical velocity magnitude -------------------------------------
u_sample = STAT.u(round(Nt/2), :, :);
u_rms = rms(u_sample(:), 'omitnan');
if u_rms > 0 && u_rms < 5
    fprintf('  [PASS] u RMS (mid-record sample) = %.4f m/s  (physically reasonable)\n', u_rms);
elseif u_rms == 0 || isnan(u_rms)
    fprintf('  [FAIL] u is all NaN or zero — data may not have loaded correctly.\n');
    ok = false;
else
    fprintf('  [WARN] u RMS = %.2f m/s — unusually large. Check DX/DT.\n', u_rms);
end

% ---- (8) FOV depth vs expected (~14 cm) ---------------------------------
fov_depth_mm = z_piv(end) * 1e3;
if fov_depth_mm > 50 && fov_depth_mm < 250
    fprintf('  [PASS] PIV depth FOV = %.1f mm\n', fov_depth_mm);
else
    fprintf('  [WARN] PIV depth FOV = %.1f mm — expected ~100–150 mm.\n', fov_depth_mm);
    fprintf('         If ~0.27 mm: GS may be wrong (grid spacing = DX not DX*GS).\n');
    fprintf('         If ~%d px:   DX not applied — check STAT.X units.\n', round(fov_depth_mm));
end

% ---- (9) z_piv(1) ≈ one grid cell (surface is at image top) -------------
expected_z1 = DX * GS;
if abs(z_piv(1) - expected_z1)/expected_z1 < 0.05
    fprintf('  [PASS] z_piv(1) = %.4f mm ≈ DX*GS (first vector one cell below surface)\n', ...
        z_piv(1)*1e3);
else
    fprintf('  [WARN] z_piv(1) = %.4f mm  (expected %.4f mm = DX*GS)\n', ...
        z_piv(1)*1e3, expected_z1*1e3);
    fprintf('         STAT.X may start at 0 instead of DX*GS — z_piv(1)=0 means surface IS in the grid.\n');
end

% ---- (10) Masking: which side of the surface is data on? -----------------
%   Fabio's PIVMask.m:  Mask(1:floor(PF_Surface(i)), i) = 1  → keeps AIR (rows above surface)
%   Your data should be WATER side: valid data BELOW the surface, NaN above.
%
%   Diagnostic: for the middle frame, check whether data is mostly above or
%   below the expected surface row.
mid_frame = squeeze(STAT.u(round(Nt/2), :, :));  % [Nx x Nz]
nan_frac_top    = mean(isnan(mid_frame(:, 1:min(5,Nz))),     'all');  % top (near surface)
nan_frac_bottom = mean(isnan(mid_frame(:, max(1,Nz-4):Nz)),  'all');  % bottom (deep)

if nan_frac_top > nan_frac_bottom
    fprintf('  [PASS] Masking: NaN concentrated near surface (top=%.0f%% NaN, bottom=%.0f%%) → WATER side data\n', ...
        nan_frac_top*100, nan_frac_bottom*100);
elseif nan_frac_bottom > nan_frac_top
    fprintf('  [WARN] Masking: NaN concentrated at depth (top=%.0f%%, bottom=%.0f%%) → might be AIR side data\n', ...
        nan_frac_top*100, nan_frac_bottom*100);
    fprintf('         If this is wrong, your data may need re-masking (see check 11).\n');
else
    fprintf('  [INFO] Masking: NaN distribution uniform (top=%.0f%%, bottom=%.0f%%) → no mask applied?\n', ...
        nan_frac_top*100, nan_frac_bottom*100);
end

% ---- (11) Apply water-side mask if needed --------------------------------
%   If data has NO mask (both sides valid) or wrong mask (air side),
%   apply water-side mask using eta to NaN everything above the surface.
%
%   Fabio (air side):   Mask(1:floor(surface), col) = 1     → keep above surface
%   Water side (ours):  Mask(ceil(surface):end, col) = 1    → keep below surface
%
%   Set apply_water_mask = true to force masking.

apply_water_mask = false;   % ← change to true if check (10) says air-side or no mask

if apply_water_mask
    fprintf('  [ACTION] Applying water-side mask to STAT.u and STAT.v...\n');
    for it = 1:Nt
        [~, ite] = min(abs(t_eta - t_piv(it)));
        eta_t = interp1(x_eta, eta(:,ite), x_piv, 'linear', 'extrap');  % [Nx x 1]

        for ix = 1:Nx
            % Surface depth in z_piv coords: z_surface = -eta_fluct
            % NaN everything shallower than the surface (above it = air)
            z_surface = -eta_t(ix) + mean(eta(:));
            air_rows = z_piv < z_surface;
            STAT.u(it, ix, air_rows) = NaN;
            STAT.v(it, ix, air_rows) = NaN;
        end
    end
    fprintf('  [DONE]  Water-side mask applied (%d frames).\n', Nt);
end

% ---- Diagnostic plot: raw velocity with NaN pattern ---------------------
figure('Name', 'Sec 2b | Masking diagnostic', 'Position', [50 600 1000 350]);
it_diag = round(Nt/2);
u_diag = squeeze(STAT.u(it_diag, :, :));  % [Nz x Nx]

subplot(1,2,1);
imagesc(x_piv*1e3, z_piv*1e3, u_diag);
set(gca, 'YDir', 'reverse');
colorbar; colormap(brewermap([],'Spectral')); caxis_sym(gca);
xlabel('x (mm)'); ylabel('z depth (mm)');
title(sprintf('u  (frame %d)', it_diag));

hold on;plot(ETA_R2_EXP1.x*1e3,-(ETA_R2_EXP1.ETA(:,it_diag)*1e3-103.42)+25.243);

subplot(1,2,2);
imagesc(x_piv*1e3, z_piv*1e3, double(isnan(u_diag)));
set(gca, 'YDir', 'reverse');
colormap(gca, [1 1 1; 0.8 0.2 0.2]);  % white=valid, red=NaN
xlabel('x (mm)'); ylabel('z depth (mm)');
title('NaN mask (red = NaN = masked out)');

sgtitle(sprintf('Section 2b — Masking diagnostic  (frame %d of %d)', it_diag, Nt));
drawnow;

% ---- Summary ------------------------------------------------------------
if ok
    fprintf('\n  ALL CRITICAL CHECKS PASSED — safe to continue.\n');
else
    fprintf('\n  ONE OR MORE CRITICAL CHECKS FAILED — fix above before continuing.\n');
end
fprintf('===================================\n\n');

%% =========================================================================
%% SECTION 3 — SELECT SNAPSHOT & PLOT CARTESIAN VELOCITY FIELD
%% =========================================================================
[~, it_snap] = min(abs(t_piv - t_snap_target));
fprintf('\nSnapshot: t = %.2f s  (frame %d of %d)\n', t_piv(it_snap), it_snap, Nt);

% squeeze(STAT.u(t,:,:)) = [Nx x Nz];  transpose → [Nz x Nx] for pcolor(x_piv, z_piv, ...)
u_cart = squeeze(STAT.u(it_snap, :, :))';   % [Nz x Nx]
w_cart = squeeze(STAT.v(it_snap, :, :))';   % [Nz x Nx]

% u(t, x, z)  →  [Nz x Nx] for plotting (z on y-axis, x on x-axis)
u_cart = squeeze(STAT.u(it_snap, :, :));   % [Nz x Nx]
w_cart = squeeze(STAT.v(it_snap, :, :));   % [Nz x Nx]

% --- Surface elevation for this frame ---
[~, it_e] = min(abs(t_eta - t_piv(it_snap)));
eta_snap = interp1(x_eta, eta(:, it_e), x_piv, 'linear', 'extrap')';  % [1 x Nx] m

% eta reference: subtract the long-term mean so that z=0 is the mean water level.
% Fabio's equivalent: DC bin of fft(surface_fab) carries the mean; zPIV is already
% measured from the mean surface. We do the same explicitly here.
eta_mean   = mean(eta(:));           % scalar, long-term spatial+temporal mean
eta_fluct  = eta_snap - eta_mean;    % [1 x Nx] wave fluctuation, mean-zero

% In the depth axis (z positive downward, z=0 at mean surface):
%   surface is at  z = -eta_fluct
%     crest (eta_fluct > 0) → z < 0 → surface is above mean (at top of/above plot)
%     trough (eta_fluct < 0) → z > 0 → surface dips below mean (visible in plot)
% Check: fprintf prints the RMS so you can verify eta is in metres, not pixels.
fprintf('  eta mean=%.4f mm,  RMS=%.4f mm  (should be ~1-3 mm for capillary-gravity)\n', ...
    eta_mean*1e3, rms(eta_fluct)*1e3);

z_surf_line = -eta_fluct;   % [1 x Nx]  depth coordinate of instantaneous surface

figure('Name', sprintf('Sec 3 | Cartesian snapshot  t=%.1f s', t_piv(it_snap)), ...
       'Position', [50 550 1200 380]);
subplot(1,2,1);
pcolor(x_piv*1e3, z_piv*1e3, u_cart); shading flat;
hold on; plot(x_piv*1e3, z_surf_line*1e3, 'k', 'LineWidth', 1.5);
colorbar; colormap_bwr(); caxis_sym(gca);

xlabel('x (mm)'); ylabel('z (mm,  +down)');
% YDir reverse: z=0 (mean surface) at top, depth increases downward — water below
set(gca, 'YDir', 'reverse');

xlabel('x (mm)'); ylabel('z (mm)'); set(gca,'YDir','reverse');
title('u_{cart}  [m/s]');

subplot(1,2,2);
pcolor(x_piv*1e3, z_piv*1e3, w_cart); shading flat;
hold on; plot(x_piv*1e3, z_surf_line*1e3, 'k', 'LineWidth', 1.5);
colorbar; colormap_bwr(); caxis_sym(gca);

xlabel('x (mm)'); ylabel('z (mm,  +down)');
set(gca, 'YDir', 'reverse');

xlabel('x (mm)'); ylabel('z (mm)'); set(gca,'YDir','reverse');
title('w_{cart}  [m/s]');
sgtitle(sprintf('Section 3 — Cartesian snapshot  t = %.1f s', t_piv(it_snap)));
drawnow;

%% =========================================================================
%% SECTION 4 — WAVE-FOLLOWING COORDINATE TRANSFORM (WATER SIDE)
%%   Air side (Fabio): z = +zeta + sum(a exp(ikx) exp(-k*zeta)),  zeta > 0 up
%%   Water side (here): z = -zeta + sum(a exp(ikx) exp(-k*zeta)), zeta > 0 down
%%   exp(-k*zeta) same: coordinate lines follow wave at surface, flatten at depth
%% =========================================================================
% Phase-average subtraction (remove_wave_phase_avg.m)
% 
% Assign wave phase φ(t) from longitudinal Hilbert at x=12m
% Bin all transverse frames by phase, average within each bin → <v,w>(y, z, φ)
% Subtract phase average from each frame
% 
% 
% dz      = abs(z_piv(2) - z_piv(1));
% zeta_vec = (0 : dz : abs(z_piv(end)))';   % [Nz x 1]  depth below surface

% Use mean-zero eta for the wave-following transform.
% The DC offset (eta_mean) shifts z=0 but doesn't create a wavy surface to flatten.
% Fabio's equivalent: fft(surface_fab) DC term sets the absolute level of h,
% but the wave undulation comes from the AC components only.
transfo = generate_transfo_water(eta_fluct, x_piv, zeta_vec);

u_wf = transform_field_to_wavefollowing(u_cart, z_piv, transfo);
w_wf = transform_field_to_wavefollowing(w_cart, z_piv, transfo);

% ---- Plot 4a: wave-following velocity ----
figure('Name', sprintf('Sec 4a | Wave-following snapshot  t=%.1f s', t_piv(it_snap)), ...
       'Position', [50 450 1200 380]);
subplot(1,2,1);
pcolor(x_piv*1e3, -zeta_vec*1e3, u_wf); shading flat;
colorbar; colormap_bwr(); caxis_sym(gca);
xlabel('x (mm)'); ylabel('-\zeta (mm)  [0=surface]'); set(gca,'YDir','normal');
title('u(x,\zeta)  [m/s]');

subplot(1,2,2);
pcolor(x_piv*1e3, -zeta_vec*1e3, w_wf); shading flat;
colorbar; colormap_bwr(); caxis_sym(gca);
xlabel('x (mm)'); ylabel('-\zeta (mm)'); set(gca,'YDir','normal');
title('w(x,\zeta)  [m/s]');
sgtitle(sprintf('Section 4a — Wave-following coordinates  t = %.1f s  (top = flattened surface)', t_piv(it_snap)));
drawnow;

% ---- Plot 4b: constant-zeta grid lines on Cartesian ----
figure('Name', 'Sec 4b | Coordinate grid on Cartesian', 'Position', [50 350 800 450]);
pcolor(x_piv*1e3, z_piv*1e3, w_cart); shading flat;
colormap_bwr(); caxis_sym(gca); hold on;
zeta_show  = [0, 0.001, 0.003, 0.007, 0.015, 0.04, 0.10];
clr = winter(length(zeta_show));
for k = 1:length(zeta_show)
    [~, iz] = min(abs(zeta_vec - zeta_show(k)));
    plot(x_piv*1e3, transfo.Z_grid(iz,:)*1e3, '-', 'Color', clr(k,:), 'LineWidth', 1.2);
end
plot(x_piv*1e3, z_surf_line*1e3, 'k', 'LineWidth', 2);  % surface at z = -eta_fluct
xlabel('x (mm)'); ylabel('z (mm,  +down)'); set(gca,'YDir','reverse');
colorbar;
legend([{'w_{cart}', '\eta(x)'}; ...
    arrayfun(@(z) sprintf('\\zeta=%.0fmm', z*1e3), zeta_show', 'Uni',0)], ...
    'Location','southwest','FontSize',7);
title('Section 4b — Constant-\zeta lines on Cartesian field');
drawnow;

%% =========================================================================
%% SECTION 5 — WAVE PHASE vs X
%% =========================================================================
figure('Name', 'Sec 5 | Wave phase vs x', 'Position', [50 280 900 300]);
yyaxis left;
plot(x_piv*1e3, transfo.phase, 'b-', 'LineWidth', 1.5);
ylabel('Phase (rad)');
ylim([-pi pi]);
set(gca, 'YTick', [-pi -pi/2 0 pi/2 pi], ...
         'YTickLabel', {'-\pi','-\pi/2','0','\pi/2','\pi'});

yyaxis right;
plot(x_piv*1e3, eta_snap*1e3, 'r-', 'LineWidth', 1);
ylabel('\eta (mm)');

xlabel('x (mm)'); grid on;
legend('\phi(x)', '\eta(x)', 'Location','best');
title('Section 5 — Wave phase from Hilbert transform');
drawnow;

%% =========================================================================
%% SECTION 6 — PHASE & ENSEMBLE AVERAGES  (loop over all frames)
%% =========================================================================
Nbins       = 180;
phase_edges = linspace(-pi, pi, Nbins+1);
phase_ctr   = 0.5*(phase_edges(1:end-1)+phase_edges(2:end));
Nz_wf       = length(zeta_vec);

u_phSum = zeros(Nz_wf, Nbins);   w_phSum = zeros(Nz_wf, Nbins);
phCount = zeros(Nz_wf, Nbins);
u_enSum = zeros(Nz_wf, 1);       w_enSum = zeros(Nz_wf, 1);
enCount = zeros(Nz_wf, 1);

fprintf('Computing phase & ensemble averages (%d frames)...\n', Nt);
for it = 1:Nt
    u_t = squeeze(STAT.u(it, :, :))';   % [Nz x Nx]
    w_t = squeeze(STAT.v(it, :, :))';

    [~, ite] = min(abs(t_eta - t_piv(it)));
    eta_t = interp1(x_eta, eta(:,ite), x_piv, 'linear','extrap')';

    tf_t  = generate_transfo_water(eta_t, x_piv, zeta_vec);
    u_t_wf = transform_field_to_wavefollowing(u_t, z_piv, tf_t);
    w_t_wf = transform_field_to_wavefollowing(w_t, z_piv, tf_t);

    bin_t = discretize(tf_t.phase, phase_edges);   % [1 x Nx]

    for ib = 1:Nbins
        cols = find(bin_t == ib);
        if isempty(cols), continue; end
        u_phSum(:,ib) = nansum([u_phSum(:,ib), nansum(u_t_wf(:,cols),2)], 2);
        w_phSum(:,ib) = nansum([w_phSum(:,ib), nansum(w_t_wf(:,cols),2)], 2);
        phCount(:,ib) = phCount(:,ib) + sum(isfinite(u_t_wf(:,cols)), 2);
    end

    u_enSum = u_enSum + nansum(u_t_wf, 2);
    w_enSum = w_enSum + nansum(w_t_wf, 2);
    enCount = enCount + sum(isfinite(u_t_wf), 2);

    if mod(it, max(1,round(Nt/10))) == 0
        fprintf('  %d%%\n', round(100*it/Nt));
    end
end

u_phAvg = u_phSum ./ phCount;         % [Nz x Nbins]
w_phAvg = w_phSum ./ phCount;
u_enAvg = u_enSum ./ enCount;         % [Nz x 1]
w_enAvg = w_enSum ./ enCount;

% ---- Plot 6a: phase-averaged fields ----
figure('Name', 'Sec 6a | Phase average', 'Position', [50 210 1200 380]);
subplot(1,2,1);
pcolor(phase_ctr, -zeta_vec*1e3, u_phAvg); shading flat;
colorbar; colormap_bwr(); caxis_sym(gca);
set(gca, 'YDir','normal', 'XTick',[-pi -pi/2 0 pi/2 pi], ...
         'XTickLabel',{'-\pi','-\pi/2','0','\pi/2','\pi'});
xlabel('Phase (rad)'); ylabel('-\zeta (mm)');
title('<u>(\zeta,\phi)  [m/s]');

subplot(1,2,2);
pcolor(phase_ctr, -zeta_vec*1e3, w_phAvg); shading flat;
colorbar; colormap_bwr(); caxis_sym(gca);
set(gca, 'YDir','normal', 'XTick',[-pi -pi/2 0 pi/2 pi], ...
         'XTickLabel',{'-\pi','-\pi/2','0','\pi/2','\pi'});
xlabel('Phase (rad)'); ylabel('-\zeta (mm)');
title('<w>(\zeta,\phi)  [m/s]');
sgtitle('Section 6a — Phase-averaged velocity (wave-following coords)');
drawnow;

% ---- Plot 6b: ensemble-averaged profiles ----
figure('Name', 'Sec 6b | Ensemble average', 'Position', [700 210 550 420]);
subplot(1,2,1);
plot(u_enAvg, -zeta_vec*1e3, 'b-', 'LineWidth', 1.5);
xlabel('<u>_{ens} (m/s)'); ylabel('-\zeta (mm)');
title('Ensemble-averaged u'); grid on; set(gca,'YDir','normal');

subplot(1,2,2);
plot(w_enAvg, -zeta_vec*1e3, 'r-', 'LineWidth', 1.5);
xlabel('<w>_{ens} (m/s)'); ylabel('-\zeta (mm)');
title('Ensemble-averaged w'); grid on; set(gca,'YDir','normal');
sgtitle('Section 6b — Ensemble average (all frames, all phases)');
drawnow;

%% =========================================================================
%% SECTION 7 — TRIPLE DECOMPOSITION FOR THE SELECTED SNAPSHOT
%%   Turbulent:     u' = u(x,zeta) - <u>_phi(zeta, phi(x))
%%   Wave-coherent: u~ = <u>_phi(zeta, phi(x)) - <u>_ens(zeta)
%% =========================================================================
bin_snap = discretize(transfo.phase, phase_edges);

u_turb = NaN(size(u_wf));   w_turb = NaN(size(w_wf));
u_wave = NaN(size(u_wf));   w_wave = NaN(size(w_wf));

for ix = 1:Nx
    ib = bin_snap(ix);
    if isnan(ib), continue; end
    u_turb(:,ix) = u_wf(:,ix) - u_phAvg(:,ib);
    w_turb(:,ix) = w_wf(:,ix) - w_phAvg(:,ib);
    u_wave(:,ix) = u_phAvg(:,ib) - u_enAvg;
    w_wave(:,ix) = w_phAvg(:,ib) - w_enAvg;
end

% ---- Plot 7: 2×2 decomposition ----
figure('Name', sprintf('Sec 7 | Decomposition  t=%.1f s', t_piv(it_snap)), ...
       'Position', [50 50 1400 700]);

subplot(2,2,1);
pcolor(x_piv*1e3, -zeta_vec*1e3, u_turb); shading flat;
colorbar; colormap_bwr(); caxis_sym(gca);
xlabel('x (mm)'); ylabel('-\zeta (mm)'); set(gca,'YDir','normal');
title("u' = u - <u>_\phi  (turbulent)");

subplot(2,2,2);
pcolor(x_piv*1e3, -zeta_vec*1e3, w_turb); shading flat;
colorbar; colormap_bwr(); caxis_sym(gca);
xlabel('x (mm)'); ylabel('-\zeta (mm)'); set(gca,'YDir','normal');
title("w' = w - <w>_\phi  (turbulent)");

subplot(2,2,3);
pcolor(x_piv*1e3, -zeta_vec*1e3, u_wave); shading flat;
colorbar; colormap_bwr(); caxis_sym(gca);
xlabel('x (mm)'); ylabel('-\zeta (mm)'); set(gca,'YDir','normal');
title('u~ = <u>_\phi - <u>_{ens}  (wave-coherent)');

subplot(2,2,4);
pcolor(x_piv*1e3, -zeta_vec*1e3, w_wave); shading flat;
colorbar; colormap_bwr(); caxis_sym(gca);
xlabel('x (mm)'); ylabel('-\zeta (mm)'); set(gca,'YDir','normal');
title('w~ = <w>_\phi - <w>_{ens}  (wave-coherent)');

sgtitle(sprintf('Section 7 — Triple decomposition  t = %.1f s', t_piv(it_snap)));
drawnow;

%% =========================================================================
%% SAVE
%% =========================================================================
save('longitudinal_wave_decomposition.mat', ...
    'u_cart','w_cart','u_wf','w_wf', ...
    'u_phAvg','w_phAvg','u_enAvg','w_enAvg', ...
    'u_turb','w_turb','u_wave','w_wave', ...
    'transfo','zeta_vec','z_piv','x_piv','t_piv', ...
    'it_snap','Nbins','phase_edges','phCount');
disp('Saved longitudinal_wave_decomposition.mat');

%% =========================================================================
%% LOCAL HELPERS
%% =========================================================================
function colormap_bwr()
% Blue–white–red diverging colormap applied to current axes.
n = 256; h = floor(n/2);
r = [linspace(0,1,h), ones(1,n-h)];
g = [linspace(0,1,h), linspace(1,0,n-h)];
b = [ones(1,h),       linspace(1,0,n-h)];
colormap(gca, [r' g' b']);
end

function caxis_sym(ax)
% Symmetric colour limits around zero.
cl = caxis(ax);
mx = max(abs(cl));
if mx > 0, caxis(ax, [-mx mx]); end
end
