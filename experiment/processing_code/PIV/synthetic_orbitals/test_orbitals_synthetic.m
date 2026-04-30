% test_orbitals_synthetic.m
% Verify that compute_orbitals_linear and generateTransfo_LC_noLFV_2023
% produce correct ORBX by testing against an analytical monochromatic wave.
%
% Analytical solution (deep-water linear theory):
%   eta(x)   = A * cos(k0*x)
%   u_orb    = A * omega0 * exp(-k0*zeta) * cos(k0*x)   [m/s, +x direction]
%   w_orb    = A * omega0 * exp(-k0*zeta) * sin(k0*x)   [m/s, +z upward]
%
% What to look for regarding current k/omega range:
%   dx      = GS * DX          % grid spacing in metres
%   Nx      = number of PIV columns
%   k_vec   = (0:Nx-1)/Nx * (2*pi/dx)   % rad/m
%   omega_v = sqrt(k_vec * g)            % rad/s  (deep-water, no cap.)
%   generateTransfo truncates at i_max = (Nx/2+1)/4  -> k_max_used ~ k_vec(i_max)

clear; clc;

%% --- Path to Fabrice's functions (adjust if needed) ---
addpath('D:\Scripps\GC-Wave-Gen\M-Files_FabMarcNovDec2014\');
addpath('D:\Scripps\GC-Wave-Gen\M-Files_FabMarcNovDec2014\FabriceScripts\');
addpath(fileparts(mfilename('fullpath')));   % this folder
addpath(fullfile(fileparts(mfilename('fullpath')), '..'));   % PIV folder

%% =========================================================================
%% 1. Grid — matches Delaware experiment
%% =========================================================================
DX  = 1/17697.69;   % m/pixel
DT  = 10e-3;        % s/frame
GS  = 4;            % pixels (grid spacing used in compVel)
dx  = GS * DX;      % m  (~2.26e-4 m)
Nx  = 256;          % PIV columns  (adjust to match real data if desired)
Nz  = 40;           % PIV rows

x        = (0:Nx-1) * dx;           % m  [1 x Nx]
zeta_vec = (1:Nz)'  * dx;           % m  [Nz x 1], positive downward

%% --- Print k/omega range so you know what the code actually sees ---
k_vec   = (0:Nx-1) / Nx * (2*pi/dx);
g       = 9.81;
gamma   = 0.073/998;   % surface tension / density  [m^3/s^2]
omega_v = sqrt(k_vec .* (g + gamma .* k_vec.^2));

i_max_transfo = floor((Nx/2+1)/4);
fprintf('--- k/omega range ---\n');
fprintf('  dx            = %.4e m\n', dx);
fprintf('  Nx            = %d\n', Nx);
fprintf('  k_max (Nyq)   = %.1f rad/m   lambda_min = %.4f m\n', k_vec(end), 2*pi/k_vec(end));
fprintf('  i_max transfo = %d  ->  k_max used = %.1f rad/m\n', i_max_transfo, k_vec(i_max_transfo+1));
fprintf('  omega at k_max_used = %.1f rad/s  (f = %.2f Hz)\n\n', ...
    omega_v(i_max_transfo+1), omega_v(i_max_transfo+1)/(2*pi));

%% =========================================================================
%% 2. Synthetic single-frequency wave
%%    Choose lambda = 5 cm  (k ~ 126 rad/m), well inside transfo range
%% =========================================================================
lambda0 = 0.05;                             % m
k0      = 2*pi / lambda0;                   % rad/m
omega0  = sqrt(k0*(g + gamma*k0^2));        % rad/s
A       = 2e-3;                             % m  (2 mm amplitude)

fprintf('--- Synthetic wave ---\n');
fprintf('  lambda = %.3f m,  k0 = %.1f rad/m,  omega0 = %.2f rad/s  (f = %.2f Hz)\n\n', ...
    lambda0, k0, omega0, omega0/(2*pi));

eta = A * cos(k0 * x);   % [1 x Nx]  in metres, mean-zero

%% =========================================================================
%% 3. Analytical solution  (ground truth)
%% =========================================================================
u_analytical = zeros(Nz, Nx);
w_analytical = zeros(Nz, Nx);
for iz = 1:Nz
    decay = exp(-k0 * zeta_vec(iz));
    u_analytical(iz,:) =  A * omega0 * decay * cos(k0*x);
    w_analytical(iz,:) =  A * omega0 * decay * sin(k0*x);
end

%% =========================================================================
%% 4. compute_orbitals_linear  (SI, clean function)
%% =========================================================================
params.g     = g;
params.gamma = gamma;
[u_orb, w_orb] = compute_orbitals_linear(eta, x, zeta_vec, params);

err_u_orb = max(abs(u_orb(:) - u_analytical(:)));
err_w_orb = max(abs(w_orb(:) - w_analytical(:)));
fprintf('--- compute_orbitals_linear ---\n');
fprintf('  max |u err| = %.2e m/s\n', err_u_orb);
fprintf('  max |w err| = %.2e m/s\n\n', err_w_orb);

if err_u_orb < 1e-8
    fprintf('  PASS: compute_orbitals_linear u matches analytical\n\n');
else
    fprintf('  FAIL: compute_orbitals_linear u mismatch — check FFT one-sided factor\n\n');
end

%% =========================================================================
%% 5. generateTransfo_LC_noLFV_2023  (pixel-space, converts to m/s via DX/DT)
%% =========================================================================
% Build minimal compVel / pivRes structs with synthetic eta
eta_pix = eta / DX;   % metres -> pixels

compVel.xPIV  = (1:Nx) * GS;
compVel.zPIV  = (1:Nz) * GS;
compVel.GS    = GS;
compVel.DX    = DX;
compVel.DT    = DT;
compVel.Mask  = ones(Nz, Nx);

pivRes        = compVel;
pivRes.mask   = compVel.Mask;

Surface_PIV   = eta_pix;   % [1 x Nx] pixel row coordinates of surface

try
    transfo = generateTransfo_LC_noLFV_2023(compVel, Surface_PIV, pivRes);

    ORBX = transfo.ORBX(2:end, :);   % remove zeta=0 row (surface itself)
    ORBZ = transfo.ORBZ(2:end, :);

    % Convert pixel/frame -> m/s
    ORBX_ms = ORBX * DX / DT;
    ORBZ_ms = ORBZ * DX / DT;

    % Trim to same Nz as analytical (transfo may have different row count)
    Nz_cmp = min(Nz, size(ORBX_ms,1));
    err_orbx = max(abs(ORBX_ms(1:Nz_cmp,:) - u_analytical(1:Nz_cmp,:)), [], 'all');
    err_orbz = max(abs(ORBZ_ms(1:Nz_cmp,:) - w_analytical(1:Nz_cmp,:)), [], 'all');

    fprintf('--- generateTransfo_LC_noLFV_2023 ---\n');
    fprintf('  max |ORBX err| = %.2e m/s\n', err_orbx);
    fprintf('  max |ORBZ err| = %.2e m/s\n\n', err_orbz);

    if err_orbx < 1e-4
        fprintf('  PASS: ORBX matches analytical\n\n');
    else
        fprintf('  FAIL: ORBX mismatch\n');
        fprintf('  Check: is k0 = %.1f rad/m within transfo k range (0 to %.1f rad/m)?\n\n', ...
            k0, k_vec(i_max_transfo+1));
    end

    has_transfo = true;
catch ME
    fprintf('  generateTransfo not found or errored: %s\n  Skipping ORBX test.\n\n', ME.message);
    has_transfo = false;
    ORBX_ms = []; ORBZ_ms = [];
end

%% =========================================================================
%% 6. Diagnostic plots
%% =========================================================================
figure('Name','Synthetic wave surface', 'Color','white');
plot(x*1e3, eta*1e3, 'b', 'LineWidth',1.5);
xlabel('x (mm)'); ylabel('\eta (mm)');
title(sprintf('Synthetic surface  \\lambda=%.0f mm  A=%.1f mm', lambda0*1e3, A*1e3));
grid on;

figure('Name','u orbital: surface row comparison', 'Color','white');
hold on;
plot(x*1e3, u_analytical(1,:)*1e3, 'k--', 'LineWidth',2,  'DisplayName','Analytical');
plot(x*1e3, u_orb(1,:)*1e3,        'b-',  'LineWidth',1.5, 'DisplayName','compute\_orbitals\_linear');
if has_transfo && ~isempty(ORBX_ms)
    plot(x*1e3, ORBX_ms(1,:)*1e3,  'r:',  'LineWidth',2,  'DisplayName','generateTransfo ORBX');
end
xlabel('x (mm)'); ylabel('u_{orb} (mm/s)');
title('u orbital at \zeta_1  (surface row)');
legend('Location','best'); grid on;

figure('Name','u orbital depth profiles', 'Color','white');
hold on;
plot(u_analytical(:, 1)*1e3, zeta_vec*1e3, 'k--', 'LineWidth',2,  'DisplayName','Analytical');
plot(u_orb(:, 1)*1e3,        zeta_vec*1e3, 'b-',  'LineWidth',1.5, 'DisplayName','compute\_orbitals\_linear');
if has_transfo && ~isempty(ORBX_ms)
    Nz_cmp = min(Nz, size(ORBX_ms,1));
    plot(ORBX_ms(1:Nz_cmp, 1)*1e3, zeta_vec(1:Nz_cmp)*1e3, 'r:','LineWidth',2,'DisplayName','generateTransfo ORBX');
end
set(gca,'YDir','reverse');
xlabel('u_{orb} (mm/s)'); ylabel('\zeta (mm, + downward)');
title('Depth decay of u orbital at x=0');
legend('Location','best'); grid on;
