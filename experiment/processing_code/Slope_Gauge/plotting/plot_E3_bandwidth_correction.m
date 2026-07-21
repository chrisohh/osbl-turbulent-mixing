function plot_E3_bandwidth_correction(data)
% E.3  Wave-bandwidth correction to u^S(z): spectral integral vs monochromatic peak.

eta = data.eta;
roi = data.setup.roi;
dx  = data.dx;

eta_roi = double(eta(roi.row_min:roi.row_max, roi.col_min:roi.col_max, :));
[Ny, Nx, Nt] = size(eta_roi);

% Time-averaged 1D E(k_x): pwelch along x for each y row, averaged
nx_win = min(128, 2^floor(log2(Nx)));
acc = []; nL = 0;
y_idx = round(linspace(1, Ny, min(20, Ny)));
for ii = 1:length(y_idx)
    for it = 1:Nt
        s = eta_roi(y_idx(ii), :, it).';
        s = s - mean(s);
        [p, fk] = pwelch(s, hann(nx_win), floor(nx_win/2), [], 1/dx);
        if isempty(acc), acc = zeros(length(fk),1); end
        acc = acc + p; nL = nL + 1;
    end
end
E_k = acc / nL / (2*pi);
k = 2*pi * fk;
k = k(2:end);  E_k = E_k(2:end);

% omega(k) deep-water
g     = 9.81;
gamma = 7.4e-5;
om_k  = sqrt(g .* k + gamma .* k.^3);

% Stokes spectral integral
z = (0:-1:-100)' * 1e-3;     % m
uS_int = zeros(length(z), 1);
for iz = 1:length(z)
    integrand = 2 .* k .* om_k .* exp(2 .* k .* z(iz)) .* E_k;
    uS_int(iz) = trapz(k, integrand);
end

% Monochromatic peak estimate
[~, ip] = max(k.^3 .* E_k);            % saturation peak
k_p   = k(ip);
om_p  = om_k(ip);
% peak amplitude from variance contribution near peak
dk = mean(diff(k));
a_p2 = max(2 * E_k(ip) * dk, eps);     % variance ~ a^2/2 -> a^2 ~ 2*E*dk
uS_mono = a_p2 .* om_p .* k_p .* exp(2 .* k_p .* z);

fig = figure('Name','E3 bandwidth correction','Position',[100 100 1100 600],'Color','w');
tiledlayout(1,2,'TileSpacing','compact','Padding','compact');

nexttile;
semilogx(uS_int,  z*1000, 'b-', 'LineWidth', 1.6); hold on;
semilogx(uS_mono, z*1000, 'r--', 'LineWidth', 1.6);
xlabel('u^S(z)  (m/s)'); ylabel('z (mm)');
legend('spectral integral','monochromatic peak','Location','best');
title(sprintf('Stokes drift profile (k_p = %.0f rad/m)', k_p));
grid on;

nexttile;
ratio = uS_int ./ max(uS_mono, eps);
plot(ratio, z*1000, 'k-', 'LineWidth', 1.6);
xline(1, 'k:');
xlabel('u^S_{int} / u^S_{mono}'); ylabel('z (mm)');
title('Bandwidth correction factor');
grid on;

sgtitle('E.3  Wave-bandwidth correction to u^S(z)');
save_figure(fig, 'E3_bandwidth_correction');
end
