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
root_path='C:\Users\airsealab\Documents\GitHub\GC-Wave-Gen';
% root_path='D:\Scripps';
addpath(strcat(root_path, '\GC-Wave-Gen\M-Files_FabMarcNovDec2014\'));
addpath(strcat(root_path,'\GC-Wave-Gen\M-Files_FabMarcNovDec2014\FabriceScripts\'));
addpath(fileparts(mfilename('fullpath')));   % this folder
addpath(fullfile(fileparts(mfilename('fullpath')), '..'));   % PIV folder

%% =========================================================================
%% 1. Grid — matches Delaware experiment
%% =========================================================================
DX  = 1/17697.69;   % m/pixel
DT  = 10e-3;        % s/frame
GS  = 4;            % pixels (grid spacing used in compVel)
dx  = GS * DX;      % m  (~2.26e-4 m)
% NOTE: Nx*dx must equal Nfull*DX (where Nfull=2048 = surface FFT length in
% generateTransfo) so both FFTs share the same bin spacing. With GS=4 and
% Nfull=2048 this requires Nx = Nfull/GS = 512 (and Nz to match).
Nx  = 511;          % PIV columns
Nz  = 511;          % PIV rows
Nfull = 2048;       % raw image columns

compVel.xPIV  = (1:Nx) * GS;
compVel.zPIV  = (1:Nz) * GS;
compVel.GS    = GS;
compVel.DX    = DX;
compVel.DT    = DT;
compVel.Mask  = ones(Nz, Nx);

x = (0:Nx-1) * dx;           % m  [1 x Nx]
xfull = (0:Nfull-1) * DX; % m [1 x N_full]
% depth coordinate: cell centres below surface, matching run_decomposition
zeta_vec = (compVel.zPIV - compVel.GS/2)' * DX;   % m  [Nz x 1]
z_ax     = zeta_vec;                                % m  [Nz x 1]

k_vec   = (0:Nx-1) / Nx * (2*pi/dx);
g       = 9.81;
gamma   = 0;%0.073/998;   % surface tension / density  [m^3/s^2]
omega_v = sqrt(k_vec .* (g + gamma .* k_vec.^3)); %capillary


%% =========================================================================
%% 2. Synthetic single-frequency wave
%%    Choose lambda = 5 cm  (k ~ 126 rad/m), well inside transfo range
%% =========================================================================
% Pick k0 to land exactly on an FFT bin: k0 = m * 2*pi/(Nfull*DX).
% Otherwise spectral leakage at high k gets amplified by huge capillary
% omega(k) ~ sqrt(gamma)*k^(3/2) and dominates the error.

m_modes = [2, 12];     % choose FFT bin indices
A_modes = [2e-3,0.5e-3];   % amplitudes
phi_modes = [0 0];% 0]%, pi/4, pi/2];      % phases (optional)

% FFT bin index (1 = fundamental)
% k0      = m_mode * 2*pi / (Nx * dx);     % rad/m
% lambda0 = 2*pi / k0;                        % m  (~5.79 cm for m=2)
% omega0  = sqrt(k0*(g + gamma*k0^3));        % rad/s
% A       = 2e-3;                             % m  (2 mm amplitude)

% fprintf('--- Synthetic wave ---\n');
% fprintf('  lambda = %.3f m,  k0 = %.1f rad/m,  omega0 = %.2f rad/s  (f = %.2f Hz)\n\n', ...
%     lambda0, k0, omega0, omega0/(2*pi));
%single 
% eta = A * cos(k0 * x);   % [1 x Nx]  in metres, mean-zero
%multiple wave modes
eta = zeros(1, Nfull); % eta in image pixel 

k_modes = zeros(size(m_modes));
omega_modes = zeros(size(m_modes));

figure;hold on
for i = 1:length(m_modes)
    k_modes(i) = m_modes(i) * 2*pi / (Nx * dx);
    omega_modes(i) = sqrt(k_modes(i) * (g + gamma*k_modes(i)^3));
    
    eta_i=A_modes(i) * cos(k_modes(i).*xfull + phi_modes(i));
    eta = eta + eta_i;
    plot(xfull,eta_i,'linewidth',1)
end
plot(xfull,eta,'linewidth',1)
xlabel('x (m)')
ylabel('\eta (m)')
legend('')
%% =========================================================================
%% 3. Analytical solution  (ground truth)
%% =========================================================================
u_analytical = zeros(Nz, Nx);
w_analytical = zeros(Nz, Nx);
% for iz = 1:Nz
%     decay = exp(-k0 * zeta_vec(iz));
%     u_analytical(iz,:) =  A * omega0 * decay * cos(k0*x);
%     w_analytical(iz,:) =  A * omega0 * decay * sin(k0*x);
% end
% Sign matches generateTransfo line 78: transfo.ORBX = -real(ORB_E)
% Delaware setup assumes waves move in -x
for iz = 1:Nz
    for i = 1:length(k_modes)
        decay = exp(-k_modes(i) * zeta_vec(iz));

        u_analytical(iz,:) = u_analytical(iz,:) - ...
            A_modes(i) * omega_modes(i) * decay .* cos(k_modes(i)*x + phi_modes(i));

        w_analytical(iz,:) = w_analytical(iz,:) - ...
            A_modes(i) * omega_modes(i) * decay .* sin(k_modes(i)*x + phi_modes(i));
    end
end
%% =========================================================================
% %% 4. compute_orbitals_linear  (SI, clean function)
% %% =========================================================================
params.g     = g;
params.gamma = gamma;
[u_orb, w_orb] = compute_orbitals_linear(eta, xfull, zeta_vec, params);
% 
% err_u_orb = max(abs(u_orb(:) - u_analytical(:)));
% err_w_orb = max(abs(w_orb(:) - w_analytical(:)));
% fprintf('--- compute_orbitals_linear ---\n');
% fprintf('  max |u err| = %.2e m/s\n', err_u_orb);
% fprintf('  max |w err| = %.2e m/s\n\n', err_w_orb);
% 
% if err_u_orb < 1e-8
%     fprintf('  PASS: compute_orbitals_linear u matches analytical\n\n');
% else
%     fprintf('  FAIL: compute_orbitals_linear u mismatch — check FFT one-sided factor\n\n');
% end

%% =========================================================================
%% 5. generateTransfo_LC_noLFV_2023  (pixel-space, converts to m/s via DX/DT)
%% =========================================================================
% Build minimal compVel / pivRes structs with synthetic eta
% eta_pix = eta / DX;   % metres -> pixels

pivRes        = compVel;
pivRes.mask   = compVel.Mask;

% generateTransfo hardcodes X = GS:GS:2048-GS (511 pts) and requires even-length
% Surface_PIV (for the Nyquist sym_vec index). Provide surface over full 2048-px width.
Surface_PIV = eta / DX;   % metres -> pixels (generateTransfo is pixel-space)

    transfo = generateTransfo_LC_noLFV_2023(compVel, Surface_PIV, pivRes);

    ORBX = transfo.ORBX(2:end, :);   % remove zeta=0 row (surface itself)
    ORBZ = transfo.ORBZ(2:end, :);

    % Convert pixel/frame -> m/s
    ORBX_ms = ORBX * DX / DT;
    ORBZ_ms = ORBZ * DX / DT;


%% =========================================================================
%% 6. Diagnostic plots
%% =========================================================================
figure;
subplot(1,2,1)
imagesc(x,z_ax,u_analytical)
colorbar; colormap(gca, brewermap([],'Spectral'))
xlabel('x (m)'); ylabel('\zeta (m)')
title('analytical u orbital (m/s)')
axis equal;axis tight
% clim([-0.04,0.04])
ylim([z_ax(1), z_ax(end)])
xlim([compVel.xPIV(1), compVel.xPIV(end)] * DX)

subplot(1,2,2)
imagesc(x,z_ax,ORBX_ms)
colorbar; colormap(gca, brewermap([],'Spectral'))
xlabel('x (m)'); ylabel('\zeta (m)')
title('computed u orbital (m/s)')
axis equal;axis tight
% clim([-0.04,0.04])
ylim([z_ax(1), z_ax(end)])
xlim([compVel.xPIV(1), compVel.xPIV(end)] * DX)
