% check_piv_particles.m
% Diagnostic: overlay PIV image pair A/B to see particle displacement,
% and draw interrogation window boxes (largest and smallest) to verify
% there are enough particles (~7+) in the smallest window.
%
% Red channel = image A, Cyan channel = image B.
% Particles that moved appear as red/cyan pairs; stationary = white.

clear; clc;

%% =========================================================================
%% PARAMETERS
%% =========================================================================
LONG = 'D:\DelawareDataBackup\Longitudinal\PIV\';
addpath('D:\Scripps\GC-Wave-Gen\M-Files_FabMarcNovDec2014\');
addpath('D:\Scripps\GC-Wave-Gen\M-Files_FabMarcNovDec2014\FabriceScripts\');
addpath('D:\Scripps\GC-Wave-Gen\M-Files_FabMarcNovDec2014\CrapperOptimizedFindSurface\');

ii = 4;
image_pair_number = 216;

IntrWndw = [256 128 64 24 16 8];
GrdSpc   = [128  64 32 12  8 4];

DX = 1/17697.69;  % m per pixel

%% =========================================================================
%% LOAD IMAGE PAIR
%% =========================================================================
DIRS = dir(LONG); DIRS = DIRS(3:end);
exp_name = DIRS(ii).name;
num_of_digits = 3;
load_path = [LONG exp_name];
pair_str = sprintf(['%0' num2str(num_of_digits) 'd'], image_pair_number);

load([load_path '/PIVRaw/PIV/' exp_name '_Piv_' pair_str '_a.mat']);
IM_a = double(imgPiv);
load([load_path '/PIVRaw/PIV/' exp_name '_Piv_' pair_str '_b.mat']);
IM_b = double(imgPiv);
[h, w] = size(IM_a);


%% =========================================================================
%% FIGURE 1 — FULL IMAGE: RED/CYAN OVERLAY + INTERROGATION WINDOWS
%% =========================================================================
% % Normalize to [0, 1] for RGB composite
% IM_a_n = IM_a / max(IM_a(:));
% IM_b_n = IM_b / max(IM_b(:));

% For the RGB composite (Figures 1 & 2), adjust per-channel before building composite:
IM_a_n = IM_a / prctile(IM_a(:), 99);  % instead of max — clips top 1% to white
IM_b_n = IM_b / prctile(IM_b(:), 99);
IM_a_n = min(IM_a_n, 1);  % clamp to [0, 1]
IM_b_n = min(IM_b_n, 1);

% Red/cyan composite: A in red, B in cyan (green+blue)
composite = cat(3, IM_a_n, IM_b_n, IM_b_n);

% Pick a center point for the window boxes (middle of image, below surface)
cx = round(w/2);
cz = round(h * 0.6);  % 60% down — should be well below surface

figure('Name', 'PIV Particle Check — Full View', 'Position', [50 100 1400 800]);
imagesc(composite); axis equal;axis tight
% For the grayscale image (Figure 3):
% imagesc(...); clim([0, prctile(IM_a(:), 95)])
hold on

% Draw largest and smallest interrogation windows
win_large = IntrWndw(1);
win_small = IntrWndw(end);

% Largest window — yellow dashed
rectangle('Position', [cx - win_large/2, cz - win_large/2, win_large, win_large], ...
    'EdgeColor', 'y', 'LineWidth', 2, 'LineStyle', '--')
text(cx - win_large/2, cz - win_large/2 - 10, ...
    sprintf('Largest: %d x %d px', win_large, win_large), ...
    'Color', 'y', 'FontSize', 12, 'FontWeight', 'bold')

% Smallest window — green solid
rectangle('Position', [cx - win_small/2, cz - win_small/2, win_small, win_small], ...
    'EdgeColor', 'g', 'LineWidth', 2)
text(cx + win_small/2 + 5, cz, ...
    sprintf('Smallest: %d x %d px', win_small, win_small), ...
    'Color', 'g', 'FontSize', 12, 'FontWeight', 'bold')

title(sprintf('Red = A, Cyan = B — %s pair %s  (center: [%d, %d] px)', ...
    exp_name, pair_str, cx, cz), 'Interpreter', 'none')
xlabel('x (px)'); ylabel('z (px)')
daspect([1, 1, 1])
drawnow

%% =========================================================================
%% FIGURE 2 — ZOOMED: SMALLEST INTERROGATION WINDOW
%%   Zoom to ~4x the smallest window so you can count particles
%% =========================================================================
pad = win_small * 2;  % show 4x the smallest window area
x_range = max(1, cx-pad) : min(w, cx+pad);
z_range = max(1, cz-pad) : min(h, cz+pad);

figure('Name', 'PIV Particle Check — Zoomed to Smallest Window', ...
    'Position', [100 100 900 900]);
imagesc(x_range, z_range, composite(z_range, x_range, :))
hold on

% Smallest window — green
rectangle('Position', [cx - win_small/2, cz - win_small/2, win_small, win_small], ...
    'EdgeColor', 'g', 'LineWidth', 2)
text(cx + win_small/2 + 2, cz - win_small/2, ...
    sprintf('%d x %d px', win_small, win_small), ...
    'Color', 'g', 'FontSize', 14, 'FontWeight', 'bold')

% Second-smallest window — cyan dashed
win_2nd = IntrWndw(end-1);
rectangle('Position', [cx - win_2nd/2, cz - win_2nd/2, win_2nd, win_2nd], ...
    'EdgeColor', 'c', 'LineWidth', 1.5, 'LineStyle', '--')
text(cx + win_2nd/2 + 2, cz - win_2nd/2, ...
    sprintf('%d x %d px', win_2nd, win_2nd), ...
    'Color', 'c', 'FontSize', 11)

% Draw pixel grid lines inside smallest window to help counting
for px = (cx - win_small/2) : (cx + win_small/2)
    plot([px px], [cz - win_small/2, cz + win_small/2], '-', ...
        'Color', [0.3 1 0.3 0.2], 'LineWidth', 0.5)
end
for pz = (cz - win_small/2) : (cz + win_small/2)
    plot([cx - win_small/2, cx + win_small/2], [pz pz], '-', ...
        'Color', [0.3 1 0.3 0.2], 'LineWidth', 0.5)
end

daspect([1, 1, 1])
title(sprintf('Zoomed — count particles in %dx%d green box (need ~7+)', ...
    win_small, win_small), 'Interpreter', 'none')
xlabel('x (px)'); ylabel('z (px)')
drawnow

%% =========================================================================
%% FIGURE 3 — OVERLAP VISUALIZATION
%%   Show how the grid spacing tiles across a region
%% =========================================================================
figure('Name', 'PIV Grid Overlap Pattern', 'Position', [200 100 900 900]);

% Zoom to a region around center, size = 2x largest window
pad_large = win_large;
x_range2 = max(1, cx-pad_large) : min(w, cx+pad_large);
z_range2 = max(1, cz-pad_large) : min(h, cz+pad_large);

imagesc(x_range2, z_range2, IM_a_n(z_range2, x_range2))
colormap gray; hold on

% Draw tiled smallest windows with current grid spacing (50% overlap)
gs_small = GrdSpc(end);
x_starts = (cx - pad_large) : gs_small : (cx + pad_large);
z_starts = (cz - pad_large) : gs_small : (cz + pad_large);

for xs = x_starts
    for zs = z_starts
        rectangle('Position', [xs, zs, win_small, win_small], ...
            'EdgeColor', [0 1 0 0.3], 'LineWidth', 0.5)
    end
end

% Highlight one window
rectangle('Position', [cx, cz, win_small, win_small], ...
    'EdgeColor', 'g', 'LineWidth', 2)

title(sprintf('Smallest window %dx%d, grid spacing %d px (%.0f%% overlap)', ...
    win_small, win_small, gs_small, (1 - gs_small/win_small)*100), ...
    'Interpreter', 'none')
xlabel('x (px)'); ylabel('z (px)')
daspect([1, 1, 1])
drawnow

%% =========================================================================
%% PRINT SUMMARY
%% =========================================================================
fprintf('\n=== PIV Interrogation Window Summary ===\n');
fprintf('  Pass   IntrWndw   GrdSpc   Overlap\n');
for k = 1:length(IntrWndw)
    overlap_pct = (1 - GrdSpc(k)/IntrWndw(k)) * 100;
    fprintf('   %d      %3d       %3d      %.0f%%\n', k, IntrWndw(k), GrdSpc(k), overlap_pct);
end
fprintf('\nSmallest window: %d x %d px = %.1f x %.1f um\n', ...
    win_small, win_small, win_small*DX*1e6, win_small*DX*1e6);
fprintf('Need ~7+ bright particles in each %dx%d box for reliable PIV.\n', ...
    win_small, win_small);
fprintf('Inspect Figure 2 to verify.\n');
