%% cisg_synthetic_amplitude_test.m
% Does the slope AMPLITUDE come out right at realistic wave displacements?
% obs(x) = ref(x + delta(x)) is synthesised from the run-95 reference frame
% (uint8, 8-bit quantisation included), then cisg_calculate_slopes is run with
% colour_map = 'linear' and 'local' and compared with the true slopes
%   S = delta_px * (pattern_scale * cm_per_pixel) / ((n_w - n_a) * L_eff).
% Cases: A sin in x, B sin in y, C broadband random-phase field (k^-3 power).
% Writes calibration/cisg_synthetic_amplitude_check.png.
clear; clc;
here = fileparts(mfilename('fullpath'));
sg = fileparts(here);
addpath(here, fullfile(sg, 'util'), sg, fullfile(fileparts(sg), 'IR Camera', 'calibration'));

%% ---- config ----
coreview_num = 95;
folder = fullfile('D:\HLAB_2026\Core', sprintf('CoreView_%d', coreview_num), 'Flare 12M125 CCL C1112A00');
ref = cisg_load_coreview(fullfile(folder, sprintf('CoreView_%d_Flare 12M125 CCL C1112A00_0001.raw', coreview_num)));
[H, W, ~] = size(ref);
S = load(fullfile(here, 'cisg_cal_CoreView_97.mat'));  cal = S.cal;
depth = cisg_apparent_depth(ref, cal);
setup = struct('n_water', 1.333, 'n_air', 1, 'use_apparent_depth', true, 'water_depth', 48, ...
               'colour_norm', 'none', 'colour_sigma', 12, 'max_disp_px', 600, 'colour_n_iter', 15, ...
               'L_eff_cm', depth.L_eff_mm / 10, 'pattern_scale', depth.r, 'target_cols', depth.edge_cols);
setup.roi = struct('row_min', round(0.2*H), 'row_max', round(0.8*H), 'col_min', round(0.2*W), 'col_max', round(0.8*W));
r = setup.roi;
[setup.cm_per_pixel_x, setup.cm_per_pixel_y] = cisg_local_scale(cal, [mean([r.col_min r.col_max]) mean([r.row_min r.row_max])]);
px2S_x = setup.pattern_scale * setup.cm_per_pixel_x / ((setup.n_water - setup.n_air) * setup.L_eff_cm);   % slope per px
px2S_y = setup.pattern_scale * setup.cm_per_pixel_y / ((setup.n_water - setup.n_air) * setup.L_eff_cm);
fprintf('slope per px of pattern shift: x %.3e, y %.3e -> rms delta 215/120 px = rms S %.3f/%.3f\n', ...
        px2S_x, px2S_y, 215*px2S_x, 120*px2S_y);
rg = @(A) A(r.row_min:r.row_max, r.col_min:r.col_max, :);

mg = 80;                                                   % margin of the synthetic block around the ROI
rows = r.row_min - mg : r.row_max + mg;  cols = r.col_min - mg : r.col_max + mg;
[X, Y] = meshgrid(single(cols), single(rows));
refd = double(ref);

%% ---- cases ----
cases = struct('name', {}, 'dx', {}, 'dy', {}, 'noise', {}, 'group', {}, 'amp', {});
gradlim = 0.65;
Aamp = [25 50 100 150 200 300 400 500];
for A = Aamp
    lam = max(1500, 2*pi*A / gradlim);
    if A == 300, lam = 3000; end
    if A == 500, lam = 4500; end
    cases(end+1) = struct('name', sprintf('A x: %d px, lam %d', A, round(lam)), 'dx', A * sin(2*pi*X / lam), ...
                          'dy', zeros(size(X), 'single'), 'noise', 0, 'group', 'A', 'amp', A); %#ok<SAGROW>
end
for k = 1:3
    A = [50 150 300];  lam = [1500 1500 3000];
    cases(end+1) = struct('name', sprintf('B y: %d px, lam %d', A(k), lam(k)), 'dx', zeros(size(X), 'single'), ...
                          'dy', A(k) * sin(2*pi*Y / lam(k)), 'noise', 0, 'group', 'B', 'amp', A(k)); %#ok<SAGROW>
end
% C: broadband random-phase, k^-3 power spectrum (smooth), on a coarse grid
q = 8;  rng(3);
[nyc, nxc] = deal(ceil(numel(rows) / q), ceil(numel(cols) / q));
[kx, ky] = meshgrid(([0:floor(nxc/2)-1, -ceil(nxc/2):-1]) / nxc, ([0:floor(nyc/2)-1, -ceil(nyc/2):-1]) / nyc);
kk = hypot(kx, ky);  amp = kk.^-1.5;  amp(kk < 1/ min(nxc, nyc) * 1.5) = 0;  amp(kk == 0) = 0;
amp(kk > 1 / (200 / q)) = 0;                                                          % no structure below ~200 px
mkfield = @(seed) real(ifft2(fft2(randn(RandStream('mt19937ar', 'Seed', seed), nyc, nxc)) .* amp));
fc = {mkfield(11), mkfield(12)};
for i = 1:2
    f = imresize(fc{i}, [numel(rows) numel(cols)], 'bicubic');
    fc{i} = single(f / std(f(:)));
end
for sc = [1 0.5 0.25]
    dxc = 215 * sc * fc{1};  dyc = 120 * sc * fc{2};
    cases(end+1) = struct('name', sprintf('C broadband x%.2f', sc), 'dx', dxc, 'dy', dyc, 'noise', 0, 'group', 'C', 'amp', sc); %#ok<SAGROW>
end
cases(end+1) = struct('name', 'C broadband x0.25 + noise', 'dx', 0.25*215*fc{1}, 'dy', 0.25*120*fc{2}, 'noise', 1.5, 'group', 'C', 'amp', 0.25);

%% ---- run ----
nc = numel(cases);  mnames = {'linear', 'local'};  short = {'lin', 'loc'};
res = struct('gain', nan(nc, 2, 2), 'corr', nan(nc, 2, 2), 'rmse', nan(nc, 2, 2), 'bias', nan(nc, 2, 2), 'mask', nan(nc, 2), 'time', nan(nc, 2));
cmp_sc = cell(nc, 2, 2);                                   % subsampled [true rec] for scatter (C cases)
for ic = 1:nc
    c = cases(ic);
    [gxx, dxy] = gradient(double(c.dx));  [dyx, gyy] = gradient(double(c.dy));   % d/dcol, d/drow
    det_j = (1 + gxx) .* (1 + gyy) - dxy .* dyx;               % x -> x + delta must not fold
    ok_inv = min(det_j(:));
    obs = ref;
    for k = 1:3
        v = interp2(refd(:,:,k), double(X + c.dx), double(Y + c.dy), 'linear');
        if c.noise > 0, v = v + c.noise * randn(size(v)); end
        obs(rows, cols, k) = uint8(round(v));
    end
    tx = rg_block(c.dx, mg) * px2S_x;  ty = rg_block(c.dy, mg) * px2S_y;
    for im = 1:2
        s = setup;  s.colour_map = mnames{im};
        t1 = tic;  [Sx, Sy] = cisg_calculate_slopes(ref, obs, s);  tm = toc(t1);
        for iv = 1:2
            rec = double(rg(Sx));  tru = tx;
            if iv == 2, rec = double(rg(Sy));  tru = ty; end
            ok = isfinite(rec);
            if std(tru(:)) < 1e-9, continue; end
            p = polyfit(tru(ok), rec(ok), 1);
            res.gain(ic, im, iv) = p(1);
            res.corr(ic, im, iv) = corr(tru(ok), rec(ok));
            res.rmse(ic, im, iv) = rms(rec(ok) - tru(ok));
            res.bias(ic, im, iv) = mean(rec(ok) - tru(ok));
            if strcmp(c.group, 'C') && c.noise == 0
                cmp_sc{ic, im, iv} = [tru(ok(1:25:end)) rec(ok(1:25:end))];
            end
        end
        res.mask(ic, im) = 100 * mean(~isfinite(rg(Sx)), 'all');
        res.time(ic, im) = tm;
    end
    fprintf('%-28s min det(1+grad d) %.2f | linear: gain x %.3f y %.3f corr %.3f %.3f | local: gain x %.3f y %.3f corr %.3f %.3f rmse %.4f %.4f mask %.1f%% (%.0f s)\n', ...
            c.name, ok_inv, res.gain(ic,1,1), res.gain(ic,1,2), res.corr(ic,1,1), res.corr(ic,1,2), ...
            res.gain(ic,2,1), res.gain(ic,2,2), res.corr(ic,2,1), res.corr(ic,2,2), res.rmse(ic,2,1), res.rmse(ic,2,2), ...
            res.mask(ic,2), res.time(ic,2));
end
save(fullfile(here, 'cisg_synthetic_amplitude_results.mat'), 'res', 'cases', 'setup');
for k = 1:nc, cases(k).dx = []; cases(k).dy = []; end

%% ---- table ----
fprintf('\ncase | method | gain x | gain y | corr x | corr y | rms err x | rms err y | bias x | bias y | masked %%\n');
for ic = 1:nc
    for im = 1:2
        fprintf('%-28s %-6s %7.3f %7.3f %7.3f %7.3f %8.4f %8.4f %+8.4f %+8.4f %6.1f\n', cases(ic).name, short{im}, ...
                res.gain(ic,im,1), res.gain(ic,im,2), res.corr(ic,im,1), res.corr(ic,im,2), res.rmse(ic,im,1), ...
                res.rmse(ic,im,2), res.bias(ic,im,1), res.bias(ic,im,2), res.mask(ic,im));
    end
end

%% ---- figure ----
fig = figure('Visible', 'off', 'Position', [100 100 1500 500]);
subplot(1,3,1);  ia = find(strcmp({cases.group}, 'A'));
plot([cases(ia).amp], res.gain(ia,1,1), 'ko-', [cases(ia).amp], res.gain(ia,2,1), 'rs-'); yline(1, 'b--');
xlabel('A (px)'); ylabel('gain (recovered / true S_x)'); legend('linear', 'local', 'ideal'); title('case A: gain vs amplitude'); grid on;
ic = find(strcmp({cases.name}, 'C broadband x1.00'), 1);  if isempty(ic), ic = find(strcmp({cases.group}, 'C'), 1); end
for im = 1:2
    subplot(1,3,1+im);  d = cmp_sc{ic, im, 1};
    if ~isempty(d), plot(d(:,1), d(:,2), '.', 'MarkerSize', 3); hold on; lim = max(abs(d(:))); plot([-lim lim], [-lim lim], 'k-'); end
    xlabel('true S_x'); ylabel('recovered S_x'); title(sprintf('case C (%s): %s', cases(ic).name, mnames{im})); axis equal; grid on;
end
exportgraphics(fig, fullfile(here, 'cisg_synthetic_amplitude_check.png'), 'Resolution', 100);
close(fig);
fprintf('Saved calibration/cisg_synthetic_amplitude_check.png\n');

function A = rg_block(F, mg)
    % ROI part of a field defined on the ROI + margin block
    A = double(F(mg + 1 : end - mg, mg + 1 : end - mg));
end
