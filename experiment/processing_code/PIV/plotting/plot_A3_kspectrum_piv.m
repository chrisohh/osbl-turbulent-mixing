function plot_A3_kspectrum_piv(data)
% A.3 (PIV part)  1D wavenumber spectrum E(k_x) and saturation B(k_x).

eta   = data.eta;
x_eta = data.x_eta;
[Nx, Nt] = size(eta);

dx = mean(diff(x_eta));
nx_win = min(128, 2^floor(log2(Nx)));

acc = []; nL = 0;
for it = 1:Nt
    s = eta(:, it); s = s - mean(s);
    [p, fk] = pwelch(s, hann(nx_win), floor(nx_win/2), [], 1/dx);
    if isempty(acc), acc = zeros(length(fk),1); end
    acc = acc + p; nL = nL + 1;
end
E_k = acc / nL / (2*pi);
k = 2*pi * fk;
k = k(2:end);  E_k = E_k(2:end);
B_k = k.^3 .* E_k;

% Noise floor: median over the upper 20% of k
ix_floor = round(0.8 * length(k)):length(k);
noise = median(E_k(ix_floor));
k_nyq = pi/dx;

fig = figure('Name','A3 wavenumber spectrum (PIV)','Position',[100 100 1200 500],'Color','w');
tiledlayout(1,2,'TileSpacing','compact','Padding','compact');

nexttile;
loglog(k, E_k, 'b-', 'LineWidth', 1.5); hold on;
yline(noise, 'k--', 'noise floor');
xline(k_nyq, 'r:',  'Nyquist');
xlabel('k_x (rad/m)'); ylabel('E(k_x) (m^3)');
title('E(k_x)'); grid on;

nexttile;
loglog(k, B_k, 'r-', 'LineWidth', 1.5); hold on;
xline(k_nyq, 'r:', 'Nyquist');
xlabel('k_x (rad/m)'); ylabel('B(k_x) = k_x^3 E(k_x)');
title('Saturation form'); grid on;

sgtitle('A.3  PIV 1D wavenumber spectrum');

var_eta = var(eta(:), 'omitnan');
fprintf('A.3 (PIV) verification: var(eta)=%.3e m^2,  int E dk = %.3e m^2\n', ...
        var_eta, trapz(k, E_k));

save_figure(fig, 'A3_kspectrum_piv');
end
