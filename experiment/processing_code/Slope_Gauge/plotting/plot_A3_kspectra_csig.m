function plot_A3_kspectra_csig(data)
% A.3  Wavenumber spectra E(k_x), E(k_y), and saturation B(k) = k^3 E(k).

eta = data.eta;
roi = data.setup.roi;
dx  = data.dx;  dy = data.dy;

eta_roi = double(eta(roi.row_min:roi.row_max, roi.col_min:roi.col_max, :));
[Ny, Nx, Nt] = size(eta_roi);

% --- E(k_x): pwelch along x for each y row, averaged over y and t ---
nx_win = min(128, 2^floor(log2(Nx)));
[Ex, kx] = avg_pwelch_dim(eta_roi, 2, nx_win, dx);
% --- E(k_y): pwelch along y for each x col, averaged over x and t ---
ny_win = min(128, 2^floor(log2(Ny)));
[Ey, ky] = avg_pwelch_dim(eta_roi, 1, ny_win, dy);

% Skip DC bin
kx = kx(2:end);  Ex = Ex(2:end);
ky = ky(2:end);  Ey = Ey(2:end);

Bx = kx.^3 .* Ex;
By = ky.^3 .* Ey;

fig = figure('Name','A3 wavenumber spectra (CSIG)','Position',[100 100 1500 450],'Color','w');
tiledlayout(1,3,'TileSpacing','compact','Padding','compact');

nexttile;
loglog(kx, Ex, 'b-', 'LineWidth', 1.5);
xlabel('k_x (rad/m)'); ylabel('E(k_x) (m^3)');
title('E(k_x)'); grid on;

nexttile;
loglog(ky, Ey, 'r-', 'LineWidth', 1.5);
xlabel('k_y (rad/m)'); ylabel('E(k_y) (m^3)');
title('E(k_y)'); grid on;

nexttile;
loglog(kx, Bx, 'b-', 'LineWidth', 1.5); hold on;
loglog(ky, By, 'r-', 'LineWidth', 1.5);
% k^{-3} reference (flat in B = k^3 E if E ~ k^{-3})
xlabel('k (rad/m)'); ylabel('B(k) = k^3 E(k)');
legend('B_x','B_y','Location','best');
title('Saturation form'); grid on;

sgtitle('A.3  1D wavenumber spectra');

% --- Verification print ---
var_eta = var(eta_roi(:), 'omitnan');
int_Ex  = trapz(kx, Ex);
fprintf('A.3 verification: var(eta) = %.3e m^2,  integral E(k_x) dk = %.3e m^2\n', var_eta, int_Ex);

save_figure(fig, 'A3_kspectra_csig');
end


function [E, k] = avg_pwelch_dim(field, dim, win, dgrid)
% Average pwelch along dim, then over the other spatial dim and time.
[Ny, Nx, Nt] = size(field);
Fs = 1/dgrid;
if dim == 2
    other = 1;
    other_idx = round(linspace(1, Ny, min(20, Ny)));
    accE = []; nL = 0;
    for ii = 1:length(other_idx)
        for it = 1:Nt
            s = squeeze(field(other_idx(ii), :, it)).';
            s = s - mean(s);
            [p, k] = pwelch(s, hann(win), floor(win/2), [], Fs);
            if isempty(accE), accE = zeros(length(k),1); end
            accE = accE + p; nL = nL + 1;
        end
    end
elseif dim == 1
    other_idx = round(linspace(1, Nx, min(20, Nx)));
    accE = []; nL = 0;
    for ii = 1:length(other_idx)
        for it = 1:Nt
            s = squeeze(field(:, other_idx(ii), it));
            s = s - mean(s);
            [p, k] = pwelch(s, hann(win), floor(win/2), [], Fs);
            if isempty(accE), accE = zeros(length(k),1); end
            accE = accE + p; nL = nL + 1;
        end
    end
end
E = accE / nL;
k = 2*pi * k;             % cyclic Hz -> rad/m
E = E / (2*pi);           % preserve variance: int E(k) dk over rad/m
end
